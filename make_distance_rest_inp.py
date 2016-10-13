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

def write_disresfile(selected_atoms, cutoff, frange, fout_name):
   fout = open(fout_name, "w")
   fout.write("RDDSTC> LIST\n")

   for i in xrange(selected_atoms.n_atoms):
     for j in xrange(i + 1, selected_atoms.n_atoms):
       iatom, jatom    = selected_atoms[i], selected_atoms[j]
       diff            = iatom.position - jatom.position
       distance        = np.linalg.norm(diff)
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
   parser.add_argument("-i", "--reference", required = True)
   parser.add_argument("-o", "--output", required = True)
   parser.add_argument("-c", "--cutoff",type=float, required = True)
   parser.add_argument("-fr","--fluctuation_range",type=float, nargs=2)
   parser.add_argument("-s", "--selection", required = True)
   args = parser.parse_args()

   #***Parse here
   ref,outfile,cutoff,frange,selection = args.reference,args.output,args.cutoff,args.fluctuation_range,args.selection
   print "Reference: ", args.reference
   print "Cutoff   : ", cutoff
   print "fluctuation range:", frange
   print "Selection: ", selection

   u = Universe(ref)
   selected_atoms = u.select_atoms(selection)

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
   write_disresfile(selected_atoms, cutoff, frange, args.output)

