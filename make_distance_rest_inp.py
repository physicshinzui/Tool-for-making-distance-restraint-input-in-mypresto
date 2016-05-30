from Bio.PDB import PDBParser
import numpy as np
import math
import sys

def GetInfoSpecifiedAtoms(structure):
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
           if atom.get_name() == "CA":
             info = [[ChainId,ResNo,Resname,atom.get_name()], atom.get_coord()]
             info_list.append(info)
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
        distances.append([info_list[i][0],info_list[j][0], dist])
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
  

def MakingInpDistanceRest(distances):
   fout = open("restreint.inp", "w")
   fout.write("RDDSTC> LIST" + "\n")
   for line in distances:
     line[0][2] = ReplacePositiveChargedResNameForTplgeneFormat(line[0][2])
     line[1][2] = ReplacePositiveChargedResNameForTplgeneFormat(line[1][2])

     group1 = "  ".join(map(str,line[0]))
     group2 = "  ".join(map(str,line[1]))
     fout.write(" " +
                group1   + " " +
                group2   + " " +
                str(1.0) + " " +
                str(1.0) + " " +
                str(5.0) + " " +
                str(12.0)+ " " +
                "YES\n")
   fout.write("RDDSTC> STOP\n")
   fout.close()

if __name__ == "__main__":
   import sys
   import argparse
   
   parser = argparse.ArgumentParser()
   parser.add_argument("-i", "--reference", required = True) 
   parser.add_argument("-c", "--cutoff") 
   args = parser.parse_args()

   ref    = args.reference
   cutoff = float(args.cutoff) 
   print "Reference: ", args.reference
   print "Cutoff   : ", cutoff

   pParser = PDBParser()
   st      = pParser.get_structure("ref_system", ref)

   InfoSpecifiedAtoms  = GetInfoSpecifiedAtoms(st)
#   coords = GetCoordinates(InfoSpecifiedAtoms)

   distancesInfo = GetInfoDistances(InfoSpecifiedAtoms, cutoff)
   #for i in xrange(len(distancesInfo)):
   #  print distancesInfo[i]
   MakingInpDistanceRest(distancesInfo)

#    
