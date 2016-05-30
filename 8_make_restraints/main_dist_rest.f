c
c
c
c
c***********************************************************
        implicit real*8 (a-h,o-z)
        parameter (nmx=100000)

        character aat*4,ars*4,achn*1
        common/ref1/aat(nmx),ars(nmx),achn(nmx)
        common/ref2/nr(nmx)
        common/ref3/xp(nmx),yp(nmx),zp(nmx)
        common/ref4/natp

        character awk*3
        character abc*3
c****
        common/rescor1/io_res_coresp
        dimension ntab_orig(nmx),ntab_new(nmx)
c***********************************************************
        print*,' '
        print*,' ##################################################'
        print*,' # Genarate distance restraint file for phygene.  #'
        print*,' # This program is, at the moment, only for       #'
        print*,' # the INTRA-CHAIN distance restraint file.       #'
        print*,' ##################################################'
        print*,' '
        print*,' ===================================== '
        print*,'  (1) Input a pdb file to generate distance restraint.'
        print*,' '

        icou=0
210     continue

        read(10,100,end=900) aaa
100     format(a6)
        if(aaa.eq."ATOM  ") then
          backspace 10
          icou=icou+1

          read(10,200,end=900) idum,aat(icou),ars(icou),achn(icou),
     *               nr(icou),xp(icou),yp(icou),zp(icou)
cc        write(6,200) idum,aat(icou),ars(icou),achn(icou),nr(icou),
cc   *               xp(icou),yp(icou),zp(icou)

200       format(6x,i5,1x,a4,1x,a4,a1,i4,4x,3f8.3)
        endif

        goto 210

900     continue

        natm=icou
        print*,'       N of all atoms in reference = ',natm
c********************
c  Residue renumbering.

        print*,' ===================================== '
        print*,'  (2) Input residue correspondence table '
        print*,'      when mismatch the residue numbering '
        print*,'      between MD model and the original pdb model. '
        print*,' '

        read(5,*) io_res_coresp
        if(io_res_coresp.eq.1) then
          print*,' '
          print*,' Input residue corresp. table '
          print*,' '
        endif
        if(io_res_coresp.ne.1) then
          print*,' '
          print*,' Do not use residue-correspondence table. '
          print*,' '
        endif
c**
        if(io_res_coresp.eq.1) then
          read(5,*) npair
          print*,' N of residue pairs: ',npair
          do ii=1,npair
            read(5,*) idum,ntab_orig(ii),ntab_new(ii)
            write(6,112) idum,ntab_orig(ii),ntab_new(ii)
112         format(i4,2x," orig: ",i4," & new: ",i4)
          enddo

          do ii=1,natm
            do jj=1,npair
              if(nr(ii).eq.ntab_orig(jj)) then
                nr(ii)=ntab_new(jj)
              endif
            enddo
          enddo
        endif
c**
c____
cc      do ii=1,natm
cc        write(6,200) ii,aat(ii),ars(ii),achn(ii),nr(ii),
cc   *                 xp(ii),yp(ii),zp(ii)
cc      enddo
c____
c********************
c  (3) Generate data.

        write(20,433)
433     format("RDDSTC> LIST")

        write(30,611)
611     format("LINE")
c**
c  Intra-chain procedure: nc = 1

        rlim=5.0d0

        nc=1
        er=1.0d0
        awk="YES"

        dfluctd=1.5d0
        dfluctu=1.5d0
c**
        do ii=2,natm

          do jj=1,ii-1
c CA
            if(aat(ii).eq." CA " .and. aat(jj).eq." CA ") then
              dx=xp(ii)-xp(jj)
              dy=yp(ii)-yp(jj)
              dz=zp(ii)-zp(jj)
              dd=dx**2 + dy**2 + dz**2
              dd=sqrt(dd)
              if(dd.gt.rlim) goto 881

              d1=dd-dfluctd
              d2=dd+dfluctu

              write(20,411) nc,nr(ii),ars(ii),aat(ii),
     *                      nc,nr(jj),ars(jj),aat(jj),
     *                      er,er,d1,d2,awk
881           continue
            endif
c CB
            if(aat(ii).eq." CB " .and. aat(jj).eq." CB ") then
              dx=xp(ii)-xp(jj)
              dy=yp(ii)-yp(jj)
              dz=zp(ii)-zp(jj)
              dd=dx**2 + dy**2 + dz**2
              dd=sqrt(dd)
              if(dd.gt.rlim) goto 882

              d1=dd-dfluctd
              d2=dd+dfluctu

              write(20,411) nc,nr(ii),ars(ii),aat(ii),
     *                      nc,nr(jj),ars(jj),aat(jj),
     *                      er,er,d1,d2,awk
882           continue
            endif
c C 
            if(aat(ii).eq." C  " .and. aat(jj).eq." C  ") then
              dx=xp(ii)-xp(jj)
              dy=yp(ii)-yp(jj)
              dz=zp(ii)-zp(jj)
              dd=dx**2 + dy**2 + dz**2
              dd=sqrt(dd)
              if(dd.gt.rlim) goto 883

              d1=dd-dfluctd
              d2=dd+dfluctu

              write(20,411) nc,nr(ii),ars(ii),aat(ii),
     *                      nc,nr(jj),ars(jj),aat(jj),
     *                      er,er,d1,d2,awk
883           continue
            endif
c O 
            if(aat(ii).eq." O  " .and. aat(jj).eq." O  ") then
              dx=xp(ii)-xp(jj)
              dy=yp(ii)-yp(jj)
              dz=zp(ii)-zp(jj)
              dd=dx**2 + dy**2 + dz**2
              dd=sqrt(dd)
              if(dd.gt.rlim) goto 884

              d1=dd-dfluctd
              d2=dd+dfluctu

              write(20,411) nc,nr(ii),ars(ii),aat(ii),
     *                      nc,nr(jj),ars(jj),aat(jj),
     *                      er,er,d1,d2,awk
884           continue
            endif


          enddo
        enddo

411     format(2(i2,1x,i4,2x,a4,2x,a4,4x),
     *         4(f8.3,1x),a3)
c***********************************************************
        stop
        end
