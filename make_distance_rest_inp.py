from Bio.PDB import PDBParser
from Bio.Data import IUPACData #.atom_weights
import numpy as np
import math
import sys

"""In here, should I include some code for simplicity? """
def GetInfoAtoms(structure): #SpecifiedAtomTypes):
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
           #if atom.get_name() in SpecifiedAtomTypes:
           info = [[ChainId,ResNo,Resname,atom.get_name()], atom.get_coord()]
           info_list.append(info)
   return info_list


def CalcDistance(coordA, coordB):
  deltaVector = coordA-coordB
  dist = math.sqrt(np.dot(deltaVector,deltaVector))
  return dist

def GetInfoDistances(info_list):
  distances_info = []
  natoms = len(info_list)
  for i in xrange(natoms):
    for j in xrange(i+1,natoms):
      dist = CalcDistance(info_list[i][1],info_list[j][1] )
      distances_info.append([info_list[i][0], info_list[j][0], dist])
  return distances_info 


def PrintError(selection):
  pass

def ExtractDistances(distances_info,cutoff, SelectionList):
  Selection1 = SelectionList[0]
  Selection2 = SelectionList[1]

  ChainId1  =Selection1[0] 
  Resi1     =Selection1[1] 
  AtomName1 =Selection1[2] 

  ChainId2  =Selection2[0] 
  Resi2     =Selection2[1] 
  AtomName2 =Selection2[2] 

  print "Selection 1: ",Selection1 
  print "Selection 2: ",Selection2 

  ListExtractedDistances = []
  for line in distances_info:
    if line[0][0] == ChainId1 and line[0][1] >= Resi1[0] and line[0][1] <= Resi1[1] and line[0][3] == AtomName1 and \
       line[1][0] == ChainId2 and line[1][1] >= Resi2[0] and line[1][1] <= Resi2[1] and line[1][3] == AtomName2 and \
       line[2] < cutoff:
      #print line
      ListExtractedDistances.append(line)
    else:
      pass
  return ListExtractedDistances

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

def MakingInpDistanceRest(DistancesInfo, FluctuationRange, filename):
   fout = open(filename, "w")
   fout.write("RDDSTC> LIST\n")
   for line in DistancesInfo:
     line[0][2] = ReplacePositiveChargedResNameForTplgeneFormat(line[0][2])
     line[1][2] = ReplacePositiveChargedResNameForTplgeneFormat(line[1][2])
     distance   = line[2]
     fluctuation_limit_list = DetermineRangeOfDistanceFluctuation(distance,FluctuationRange[0], FluctuationRange[1]) 
     LowerLim = fluctuation_limit_list[0]
     UpperLim = fluctuation_limit_list[1]

     group1 = line[0]
     group2 = line[1]

     #group1 = map(str,line[0])
     #group2 = map(str,line[1])

     #group1 = "  ".join(map(str,line[0]))
     #group2 = "  ".join(map(str,line[1]))

     if LowerLim < 0:
       print group1, group2, distance
       sys.exit("???Lower limit is negative. STOP???")

     fout.write("{0:>2}{1:>5}  {2:<5}{3:>4} ".format(
                group1[0],\
                group1[1],\
                group1[2],\
                group1[3]))
     fout.write("{0:>2}{1:>5}  {2:<5}{3:>4}  ".format(
                group2[0],\
                group2[1],\
                group2[2],\
                group2[3]))
     fout.write("{0:.1f}  {1:.1f}  ".format(
                1.0,\
                1.0))
     fout.write("%10.3f %10.3f" %(LowerLim,UpperLim))
     fout.write("  {0:>3}".format("YES\n"))

   fout.write("RDDSTC> STOP\n")
   fout.close()

def ReadSelectionList(filename):
  fin = open(filename, "r")
  SelectionList = [line.rstrip().split(",") for line in fin]
  for i, line in enumerate(SelectionList):
      SelectionList[i][1] = map(int,line[1].split("-"))
      SelectionList[i][2] = line[2].strip()
  return SelectionList


if __name__ == "__main__":
   import sys
   import argparse
   
   parser = argparse.ArgumentParser()
   parser.add_argument("-i", "--reference", required = True) 
   parser.add_argument("-o", "--output", required = True) 
   parser.add_argument("-c", "--cutoff",type=float, required = True) 
   parser.add_argument("-fr","--fluctuation_range",type=float, nargs=2) 
   parser.add_argument("-s", "--selection_parameter_file", required = True) 
   args = parser.parse_args()

   #***Parse here
   ref      = args.reference
   outfile  = args.output
   cutoff   = args.cutoff
   frange   = args.fluctuation_range
   SelectionParamFile = args.selection_parameter_file
   print "Reference: ", args.reference
   print "Cutoff   : ", cutoff
   print "fluctuation range:", frange
   print "Selection list: ",SelectionParamFile
   SelectionList = ReadSelectionList(SelectionParamFile)

   pParser = PDBParser()
   st      = pParser.get_structure("ref_system", ref)

   #***Store informatino of all atoms of a refrerence PDB
   InfoAtoms  = GetInfoAtoms(st)
   print "#DBG No of atoms   :", len(InfoAtoms)

   #***
   distancesInfo = GetInfoDistances(InfoAtoms)
   print "#DBG NO of pairs   :", len(distancesInfo)

   #***
   extracted_distances = ExtractDistances(distancesInfo, cutoff, SelectionList)
   print "#DBG No of pairs in cutoff distance:",len(extracted_distances)

   MakingInpDistanceRest(extracted_distances, frange, outfile)
