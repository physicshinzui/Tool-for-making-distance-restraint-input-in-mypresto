from Bio.PDB import PDBParser
from Bio.Data import IUPACData #.atom_weights
import numpy as np
import math
import sys

""" NOTE:
In *.list, selection have to be arranged to be molecule-order based on tpl file.
Residue number of each chain of the inupt PDB must be arranged from 1 to N.

e.g.-------------------
chain A
ATOM      1  CH3 ACE A   1       5.561  -4.399 -12.716 12.01 -0.37
...
ATOM      6  O   ACE A   1       5.373  -5.278 -14.918 16.00 -0.57
ATOM      7  N   MET A   2       3.501  -4.426 -13.987 14.01 -0.42
...
ATOM     12  HB2 MET A   2       2.628  -6.782 -14.514  1.01  0.02
.
.
.
ATOM   1456  CA  GLU-A  93      17.755  -4.468  -8.382 12.01  0.04
...
ATOM   1469  N   NHE A  94      17.511  -5.673 -10.516 14.01 -0.55
ATOM   1470  H2  NHE A  94      17.107  -6.420 -11.061  1.01  0.28
ATOM   1471  H3  NHE A  94      18.212  -5.062 -10.899  1.01  0.28

chain B
ATOM   1472  CH3 ACE B   1      24.160  33.215   6.267 12.01 -0.37
.
.
.
-------------------
"""

def GetInfoAtoms(structure):
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
           info = [[ChainId,ResNo,Resname,atom.get_name()], atom.get_coord()]
           info_list.append(info)
   return info_list


def CalcDistance(coordA, coordB):
  deltaVector = coordA - coordB
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

   #***Store information of all atoms of a refrerence PDB
   InfoAtoms  = GetInfoAtoms(st)
   print "#DBG No of atoms   :", len(InfoAtoms)

   #***Store distances of all atoms (This can make calculation speed slow.)
   distancesInfo = GetInfoDistances(InfoAtoms)
   print "#DBG NO of pairs   :", len(distancesInfo)

   #***
   extracted_distances = ExtractDistances(distancesInfo, cutoff, SelectionList)
   print "#DBG No of pairs in cutoff distance:",len(extracted_distances)
   print "#---Extracted distances---"
   for line in extracted_distances:
     print line

   MakingInpDistanceRest(extracted_distances, frange, outfile)
