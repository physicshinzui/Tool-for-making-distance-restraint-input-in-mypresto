# Tool-for-making-distance-restraint-input-in-mypresto

This program depends on BioPython.

Note:
Do not use chain id A,B,C,..., but 1,2,3,...whose order must be arranged based on a tpl file.

Example:
python make_distance_rest_inp.py --reference merged.pdb --cutoff 10.0 --fluctuation_range 1.5 4.5  -s selection.list -o restraint.inp

