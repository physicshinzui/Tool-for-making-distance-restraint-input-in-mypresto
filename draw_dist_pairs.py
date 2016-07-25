"""
23 June 2016
This was made by sinzi
"""
import sys
from pymol import cmd

def error_deal():
  pass

def draw_dist_pairs(drs_inp, molecule, draw_type="line"):
  f_drsinp = open(drs_inp, "r")
  obj = "network"
  cmd.create(obj, molecule)
  cmd.hide("everything", obj)
  for i, line in enumerate(f_drsinp):
    print "Pair No.: ", i
    if line.split()[1] != "LIST" and line.split()[1] != "STOP":
      chainid1, resid1, resname1, atomname1 = line.split()[0], line.split()[1], line.split()[2],line.split()[3]
      chainid2, resid2, resname2, atomname2 = line.split()[4], line.split()[5], line.split()[6],line.split()[7]
      print chainid1, resid1, resname1, atomname1, chainid2, resid2, resname2, atomname2
      if draw_type == "line":
        cmd.bond("%s and chain %s and resi %s and name %s"%(obj,chainid1,resid1,atomname1),\
                 "%s and chain %s and resi %s and name %s"%(obj,chainid2,resid2,atomname2) )
        cmd.show("lines", "%s and chain %s and resi %s and name %s"%(obj, chainid1,resid1,atomname1))
        cmd.show("lines", "%s and chain %s and resi %s and name %s"%(obj, chainid2,resid2,atomname2))
      elif draw_type == "dist":
        cmd.distance("dist%s"%i, "chain %s and resi %s and name %s"%(chainid1,str(resid1),atomname1),\
                                 "chain %s and resi %s and name %s"%(chainid2,str(resid2),atomname2) )
      print "chain %s and resi %s and name %s"%(chainid1,resid1,atomname1)
    elif line.split()[1] == "STOP":
      print "Reached last line"
      break

cmd.extend("draw_dist_pairs", draw_dist_pairs)
