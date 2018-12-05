
# NOTE: you need to have rmsd installed
# for the rmsd module to work

# NOTE: this function is going to change since I now
# know how to call the pieces of the relevant rmsd
# program (without piping input from the terminal).

from subprocess import run, PIPE

def calc_rmsd(f1, f2):
    """Calculate the root-mean-square deviation.

    Parameters
    ----------
    f1 : str
        The filename for the first geometry. Should
        be xyz or pdb.
    f2 : str
        The filename for the second geometry. Should
        be xyz or pdb.

    Returns
    -------
    rmsd : float
        The rmsd for the two molecular geometries.
        This value is after the rotation matrix
        is calculated and applied.

    """
    r = run(['calculate_rmsd', f1, f2], stdout=PIPE)
    output = r.stdout
    output = str(output)[2:]
    output = output.replace('\\n', ', ')
    output = output[:-3]
    return float(output)
