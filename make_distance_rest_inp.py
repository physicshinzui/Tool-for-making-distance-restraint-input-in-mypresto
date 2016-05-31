from Bio.PDB import PDBParser
from Bio.Data import IUPACData #.atom_weights
import numpy as np
import math
import sys

def GetInfoSpecifiedAtoms(structure, SpecifiedAtomTypes=["CA"]):
   info_list = []
   for model in structure:
     for chain in model:
       ChainId = chain.id
       #if ChainId ==" ":
       #  ChainId == "1" 
       for residue in chain:
	 Resname = residue.resname
         ResNo   = residue.id[1]
         for atom in residue:
           if atom.get_name() in SpecifiedAtomTypes:
             info = [[ChainId,ResNo,Resname,atom.get_name()], atom.get_coord()]
             info_list.append(info)
	   else:
	    # print "Ignored atom:", atom.get_name()
	     pass
   return info_list

def GetCoordinates(info_list):
    NoOfAtoms   = len(info_list)
    coordinates = [InfoSpecifiedAtoms[i][1] for i in xrange(NoOfAtoms)]
    return coordinates
  

def calc_distance(coordA, coordB):
  deltaVector = coordA-coordB
  dist = math.sqrt(np.dot(deltaVector,deltaVector))
  return dist

def GetInfoDistances(info_list, cutoff):
  distances = []
  natoms = len(info_list)
  for i in xrange(natoms):
    for j in xrange(i+1,natoms):
      dist = calc_distance(info_list[i][1],info_list[j][1] )
      if dist < cutoff:
        distances.append([info_list[i][0],info_list[j][0],dist])
  return distances

def ReplacePositiveChargedResNameForTplgeneFormat(Resname):
   if Resname == "ARG":
      Resname =  "ARG+"
   if Resname == "LYS":
      Resname =  "LYS+"
   if Resname == "ASP":
      Resname =  "ASP-"
   if Resname == "GLU":
      Resname =  "GLU-"
   return Resname
  
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

def MakingInpDistanceRest(DistancesInfo, FluctuationRange):
   fout = open("restreint.inp", "w")
   fout.write("RDDSTC> LIST" + "\n")
   for line in DistancesInfo:
     line[0][2] = ReplacePositiveChargedResNameForTplgeneFormat(line[0][2])
     line[1][2] = ReplacePositiveChargedResNameForTplgeneFormat(line[1][2])
     distance   = line[2]
     fluctuation_limit_list = DetermineRangeOfDistanceFluctuation(distance,FluctuationRange[0], FluctuationRange[1]) 
     LowerLim = fluctuation_limit_list[0]
     UpperLim = fluctuation_limit_list[1]

     group1 = "  ".join(map(str,line[0]))
     group2 = "  ".join(map(str,line[1]))

     if LowerLim < 0:
       print group1, group2, distance
       sys.exit("???Lower limit is negative.???")

     fout.write(" " +
                group1   + " " +
                group2   + " " +
                str(1.0) + " " +
                str(1.0) + " " +
                str(LowerLim) + " " +
                str(UpperLim) + " " +
                "YES\n")
   fout.write("RDDSTC> STOP\n")
   fout.close()

if __name__ == "__main__":
   import sys
   import argparse
   
   parser = argparse.ArgumentParser()
   parser.add_argument("-i", "--reference", required = True) 
   parser.add_argument("-c", "--cutoff", required = True) 
   parser.add_argument("-fr","--fluctuation_range",   nargs=2) 
   parser.add_argument("-at","--specified_atom_type", nargs="*", required = True) 
   args = parser.parse_args()

   ref      = args.reference
   cutoff   = float(args.cutoff) 
   frange   = map(float,args.fluctuation_range)
   atomtype = args.specified_atom_type
   print "Reference: ", args.reference
   print "Cutoff   : ", cutoff
   print "fluctuation range:", frange

   pParser = PDBParser()
   st      = pParser.get_structure("ref_system", ref)

   InfoSpecifiedAtoms  = GetInfoSpecifiedAtoms(st, atomtype)

   distancesInfo = GetInfoDistances(InfoSpecifiedAtoms, cutoff)
   MakingInpDistanceRest(distancesInfo, frange)
