# general purpose tools or ways to use kaplan

from kaplan.inputs import Inputs
from kaplan.pmem import Pmem
import os
import csv
# profiling
import cProfile, pstats, io


# this is where the test files are located
TEST_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test/testfiles")



def generate_data(name, num_iter, **kwargs):
    """Generate a dataset from a molecule name.
    
    Parameters
    ----------
    name : str
        The name of the molecule to run. Should
        be in pubchem. Can be overwritten in
        kwargs.
    num_iter : int
        Number of pmems to write.
    kwargs : keyword arguments.
        Gets passed to the inputs object.

    Notes
    -----
    Writes an output file called dihedrals_energies.csv
    to the output directory. The num_geoms is set to
    1 by default.

    Returns
    -------
    None
    
    """
    inputs = Inputs()
    input_dict = {
        "struct_input": name,
        "num_geoms": 1,
    }
    for key, value in kwargs.items():
        input_dict[key] = value
    inputs.update_inputs(input_dict)
    outfile_name = "dihedrals_energies.csv"
    with open(os.path.join(inputs.output_dir, outfile_name), "w") as f:
        fcsv = csv.writer(f)
        header = ["dihed"+i for i in range(inputs.num_dihed)] + ["energy"]
        fcsv.writerow(header)
        for i in num_iter:
            p = Pmem(None, 0, inputs.num_geoms, inputs.num_dihed)
            for i, dihedrals in enumerate(p):
                p.set_energy_get_coords(i)
                row = [d for d in dihedrals] + p.energies[i]
                fcsv.writerow(row)



# profiling tools

def profile_function(func_name, output_file, *args, **kwargs):
    """Run a profile on a function, writing output to output_file.
    
    Notes
    -----
    The output_file is not readable (I believe it is a pstats object).
    To read the output in plain text, call the analyse_profile function
    in this module, where the first argument (profile) is the
    name of the output_file.
    
    """
    pr = cProfile.Profile()
    pr.enable()
    func_name(*args, **kwargs)
    pr.disable()
    s = io.StringIO()
    sortby = pstats.SortKey.CUMULATIVE
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    if output_file is None:
        ps.print_stats()
    else:
        ps.dump_stats(output_file)


def analyse_profile(profile, output_file, sortby="cumulative"):
    """Generate a plain text output of a profile.

    Parameters
    ----------
    profile : str
        The path/filename of the original pstats profile.
    output_file : str
        The path/filename to write of the profile.
    sortby : str
        How the profile should be ordered. Defaults to
        cumulative.
    
    Returns
    -------
    None

    """
    with open(output_file, "w") as fout:
        ps = pstats.Stats(profile, stream=fout)
        ps.strip_dirs().sort_stats(sortby).print_stats()
