from MDAnalysis import Universe
import numpy as np
import math
import sys

def DetermineRangeOfDistanceFluctuation(distance, LowerLimit=1.5, UpperLimit=3.5):
  if distance < LowerLimit:
    print "##########################################################################"
    print "The lower bound is higher than a distance, which must be lower to be positive."
    print "The distance =", distance
    print "##########################################################################"
    sys.exit()
  distanceUpperLim = distance + UpperLimit
  distanceLowerLim = distance - LowerLimit
  ListLimit = [distanceLowerLim, distanceUpperLim]
  return ListLimit



def write_disresfile(molsele_type, list_selected_atoms, cutoff, frange, fout_name):
   fout = open(fout_name, "w")
   fout.write("RDDSTC> LIST\n")

   if molsele_type == "inter":
     group1, group2 = list_selected_atoms[0], list_selected_atoms[1]
     for i in xrange(group1.n_atoms):
       for j in xrange(group2.n_atoms):
         iatom, jatom    = group1[i], group2[j]
         diff            = iatom.position - jatom.position
         distance        = np.linalg.norm(diff)
   #      print iatom.name, jatom.name, distance 
         distance_mergin = DetermineRangeOfDistanceFluctuation(distance, frange[0], frange[1])
         low_lim, up_lim = distance_mergin[0], distance_mergin[1]
         if distance < cutoff: 
           fout.write("{0:>2}{1:>5}  {2:<5}{3:>4} ".format(
                       iatom.segid, iatom.resid, iatom.resname, iatom.name))
           fout.write("{0:>2}{1:>5}  {2:<5}{3:>4} ".format(
                       jatom.segid, jatom.resid, jatom.resname, jatom.name) )
           fout.write(" {0:.1f}  {1:.1f}  ".format(1.0,1.0) )
           fout.write("%10.3f %10.3f" %(low_lim,up_lim))
           fout.write("  {0:>3}".format("YES\n"))

   elif molsele_type == "intra":
     group1 = list_selected_atoms[0]
     for i in xrange(group1.n_atoms):
       for j in xrange(i+1, group1.n_atoms):
         iatom, jatom    = group1[i], group1[j]
         diff            = iatom.position - jatom.position
         distance        = np.linalg.norm(diff)
    #     print iatom.name, jatom.name, distance 
         distance_mergin = DetermineRangeOfDistanceFluctuation(distance, frange[0], frange[1])
         low_lim, up_lim = distance_mergin[0], distance_mergin[1]
         if distance < cutoff: 
           fout.write("{0:>2}{1:>5}  {2:<5}{3:>4} ".format(
                       iatom.segid, iatom.resid, iatom.resname, iatom.name))
           fout.write("{0:>2}{1:>5}  {2:<5}{3:>4} ".format(
                       jatom.segid, jatom.resid, jatom.resname, jatom.name) )
           fout.write(" {0:.1f}  {1:.1f}  ".format(1.0,1.0) )
           fout.write("%10.3f %10.3f" %(low_lim,up_lim))
           fout.write("  {0:>3}".format("YES\n"))

   fout.write("RDDSTC> STOP\n")
   fout.close()
 
   
if __name__ == "__main__":
   import sys
   import argparse

   parser = argparse.ArgumentParser()
   parser.add_argument("-t", "--intra_inter_type", required = True)
   parser.add_argument("-i", "--reference", required = True)
   parser.add_argument("-o", "--output", required = True)
   parser.add_argument("-c", "--cutoff",type=float, required = True)
   parser.add_argument("-fr","--fluctuation_range",type=float, nargs=2)
   parser.add_argument("-s", "--selection", nargs="*", required = True)
   args = parser.parse_args()

   #***Parse here
   molsele_type, ref,outfile,cutoff,frange,selection = args.intra_inter_type, args.reference,args.output,args.cutoff,args.fluctuation_range,args.selection
   print "intra inter type: ", molsele_type
   print "Reference: ", ref
   print "Cutoff   : ", cutoff
   print "fluctuation range:", frange
   print "Selection: ", selection

   if molsele_type != "intra" and molsele_type != "inter": sys.exit("\n specify intra/inter \n")
   if molsele_type == "inter" and len(selection) <  2: sys.exit("\n if you select 'inter', two selections are required.\n")
   if molsele_type == "intra" and len(selection) != 1: sys.exit("\n if you select 'intra', one selection is required.\n")

   u = Universe(ref)
   list_selected_atoms = []
   for isele in selection:
     list_selected_atoms.append(u.select_atoms(isele))
   print list_selected_atoms

   for selected_atoms in list_selected_atoms:
#Error treatment
     if len(selected_atoms.segids) == 0: sys.exit("\n NOTE) STOP. There is no chain ID in the input PDB.\n")
     for iatom in selected_atoms:
       try:
         val = int(iatom.segid)
       except ValueError:
         print "\n NOTE) STOP. There is a Chain ID specified by char (=%s)"%iatom.segid
         print "       You must specify chain IDs on the basis of a topology file of presto.\n"
         sys.exit()

#***write a distance restrain file
   write_disresfile(molsele_type, list_selected_atoms, cutoff, frange, args.output)
   #write_disresfile(selected_atoms, cutoff, frange, args.output)
   sys.exit()

