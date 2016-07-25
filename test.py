from pymol import cmd

def test():
  cmd.distance("test","resi %s and name %s"%(2, "CA") , "resi %s and name %s"%(3, "CA") )

cmd.extend("test", test)
