

"""
Consider the conversion of cartesian coordinates
to a z-matrix.

Likely to achieve using openbabel.

Need to clarify: does a z-matrix contain all
potential dihedral angles in a molecule?

Need to parse a coordinate system and generate
a list of dihedral angles.

Vetee already has the ability to read a z-matrix
file and write a z-matrix file using openbabel.
All that remains to be done is read a z-matrix
file and pull out dihedrals.


"""




def get_dihedrals(input_file):
    with open(input_file, 'r') as f:
        
        
    
