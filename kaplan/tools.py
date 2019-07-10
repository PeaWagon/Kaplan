# general purpose tools or ways to use kaplan
# includes analysis tools for kaplan data

from vetee.coordinates import read_xyz
from kaplan.inputs import Inputs, InputError
from kaplan.pmem import Pmem
from kaplan.geometry import create_obmol
from kaplan.inputs import Inputs

import os
import csv
import pybel
import openbabel
import numpy as np

# plotting/heatmaps
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.lines import Line2D
from scipy.stats import linregress

# profiling
import cProfile
import pstats
import io


__all__ = [
    "amino_acids",
    "TEST_DIR",
    "units_by_prog",
    "constants",
    "energy_units",
    "get_bonds_list",
    "make_2d",
    "plot_2d",
    "generate_data",
    "profile_function",
    "analyse_profile",
    "energy_barplot",
    "energy_rmsd_scatter",
    "dihedrals_heatmap",
    "make_heatmap",
]


amino_acids = [
    "asparagine", "glutamine", "aspartate", "glycine",
    "tryptophan", "cysteine", "threonine", "alanine",
    "isoleucine", "leucine", "tyrosine", "glutamate",
    "proline", "histidine", "lysine", "serine",
    "arginine", "valine", "methionine", "phenylalanine",
]

# this is where the test files are located
TEST_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test/testfiles")

# the default units used by each program (for plotting purposes)
units_by_prog = {
    "psi4": "Ha",
    "openbabel": "kcalmol-1",
}

# http://www.psicode.org/psi4manual/1.1/autodoc_physconst.html?highlight=bohr
# https://www.convertunits.com/from/hartree/to/joule
constants = {
    "avogadro's number": 6.02214179e23,  # units mol-1
}

# usage example for energy_units:
# result = energy_units["Ha"]["kJmol-1"](1)
energy_units = {
    "Ha": {                                         # hartrees
        "Ha": lambda x: x,
        "kJmol-1": lambda x: x * 2625.500,          # kilojoules/mole
        "J": lambda x: x * 4.359744E-18,            # joules
        "eV": lambda x: x * 27.21138,               # electron-volts
        "kcal": lambda x: x * 1.0415567394524e-21,  # kilocalories
        "kcalmol-1": lambda x: x * 627.5095,        # kilocalories/mole
        "cm-1": lambda x: x * 219474.6,             # wavenumbers
        "MHz": lambda x: x * 6.579684e9,            # megahertz
    },
    "kJmol-1": {
        "kJmol-1": lambda x: x,
        "Ha": lambda x: x / 2625.500,
        "kcalmol-1": lambda x: x / 4.1858,
    },
    "J": {
        "J": lambda x: x,
        "Ha": lambda x: x / 4.359744E-18,
    },
    "eV": {
        "eV": lambda x: x,
        "Ha": lambda x: x / 27.21138,
    },
    "kcal": {
        "kcal": lambda x: x,
        "Ha": lambda x: x / 1.0415567394524e-21,
    },
    "kcalmol-1": {
        "kcalmol-1": lambda x: x,
        "Ha": lambda x: x / 627.5095,
        "kJmol-1": lambda x: x * 4.1858,
    },
    "cm-1": {
        "cm-1": lambda x: x,
        "Ha": lambda x: x / 219474.6,
    },
    "MHz": {
        "MHz": lambda x: x,
        "Ha": lambda x: x / 6.579684e9,
    },
}


def get_bonds_list(obmol):
    """Get a list of bonds from an openbabel molecule.

    Parameters
    ----------
    obmol : openbabel.OBMol object

    Returns
    -------
    bonds_list : list()
        Each element is a tuple of 3 integers:
        atom1 index, atom2 index, and bond order.

    """
    # this is how we iterate over the bonds in the
    # molecule and find out which atoms are connected
    # to one another
    bonds_list = []
    for obbond in openbabel.OBMolBondIter(obmol):
        # bond_len = obbond.GetLength()
        atom1 = obbond.GetBeginAtom().GetIdx() - 1
        atom2 = obbond.GetEndAtom().GetIdx() - 1
        bond_order = obbond.GetBondOrder()
        bonds_list.append((atom1, atom2, bond_order))
    return bonds_list


def make_2d():
    """Use openbabel to make a 2D image of the obmol.

    Notes
    -----
    The input parameter obmol is required to be set.
    Coordinates are made into 2D coordinates by openbabel.
    Note: setting update (to update the obmol's coordinates)
    to True does not work unless obmol does not have any
    hydrogens. The drawing function by default removes
    hydrogens.

    Returns
    -------
    None

    """
    inputs = Inputs()
    mol = pybel.Molecule(OBMol=inputs.obmol)
    mol.draw(
        filename=os.path.join(inputs.output_dir, "obmol2d.png"),
        show=False,
        update=False,
        usecoords=False
    )


def plot_2d(xyzfile=None):
    """Make a plot of coordinates in 2-dimensions given an xyz file.

    Parameters
    ----------
    xyzfile : str or None
        Defaults to None, in which case the output_dir is searched
        for input_coords.xyz. If a string is provided, then this
        function attemtps to read in the string as an input file
        and path.

    Notes
    -----
    Writes plot2d.png to the output_dir. Flattens the xyz coordinates
    such that only the x and y coordinates are shown. Each coordinate
    is labelled with the atom number and the type of atom (example: H).
    The purpose of this plot is to allow the user to crudely identify
    the dihedral angles of their input molecule. This function does
    not account for atom overlap, and so labels may become obscured.

    Inputs that must be set:
    * charge
    * multip
    * output_dir
    These values must correspond to the input file (if given).

    Returns
    -------
    None

    """
    # read in input data
    inputs = Inputs()
    if xyzfile is None:
        xyzfile = os.path.join(inputs.output_dir, "input_coords.xyz")
    data = read_xyz(xyzfile)
    xs = [data["_coords"][i][1] for i in range(data["_num_atoms"])]
    ys = [data["_coords"][i][2] for i in range(data["_num_atoms"])]
    labels = [f"{i},{data['_coords'][i][0]}" for i in range(data["_num_atoms"])]

    # generate basic scatterplot with atom labels
    _, ax = plt.subplots()
    plt.scatter(xs, ys, c="r")
    for i, x in enumerate(xs):
        label = labels[i]
        plt.annotate(
            label,  # this is the text
            (x, ys[i]),  # this is the point to label
            textcoords="offset points",  # how to position the text
            xytext=(0, 10),  # distance from text to points (x,y)
            ha="center"  # horizontal alignment can be left, right or center
        )

    # add lines where openbabel has determined there should be bonds
    obmol = create_obmol(xyzfile, inputs.charge, inputs.multip)
    bonds = get_bonds_list(obmol)
    for bond in bonds:
        x_coord = [xs[bond[0]], xs[bond[1]]]
        y_coord = [ys[bond[0]], ys[bond[1]]]
        bond_line = Line2D(x_coord, y_coord)
        ax.add_line(bond_line)

    plt.savefig(os.path.join(inputs.output_dir, "plot2d.png"), dpi=200)
    plt.close()


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
        header = ["dihed" + i for i in range(inputs.num_dihed)] + ["energy"]
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


def energy_barplot(pmem, inunits="Ha", outunits="kJmol-1", scale=True):
    """Create a barplot of energies using matplotlib.

    Parameters
    ----------
    pmem : kaplan.pmem.Pmem object
        Containts a list of energy values to plot.
    inunits : str
        The input units.
    outunits : str
        The output units to display on the plot.
    scale : bool
        Scale the y-axis to the energy values in the range
        95% of minimum energy to 105% of maximum energy.

    Returns
    -------
    None

    """
    energies = []
    labels = []
    # ignore None energies
    # create bar labels
    for i, energy in enumerate(pmem.energies):
        if energy is not None:
            labels.append(f"{i}")
            # convert energy units, otherwise raise error
            try:
                energies.append(energy_units[inunits][outunits](energy))
            except KeyError:
                raise ValueError(f"Unavailable unit conversion:\n{inunits} to {outunits}")
    # if no energies are present, raise error
    if energies == []:
        raise InputError("Pmem does not have any valid energies.")

    # make plot
    inputs = Inputs()
    x = range((inputs.num_geoms - (inputs.num_geoms - len(energies))))
    plt.bar(height=energies, x=x, tick_label=labels)
    plt.xlabel("Conformer Number")
    plt.ylabel(f"Energy ({outunits})")
    plt.title(f"Conformer Energies\n{inputs.struct_input}, Pmem{pmem.ring_loc}")

    if scale:
        min_scale = min(energies)
        max_scale = max(energies)
        if min_scale < 0:
            min_scale = 1.00001 * min_scale
        elif min_scale == 0:
            min_scale = -(max_scale - min_scale) * 0.05
        else:
            min_scale = 0.99999 * min_scale
        if max_scale < 0:
            max_scale = 0.99999 * max_scale
        else:
            max_scale = 1.00001 * max_scale
        axes = plt.gca()
        axes.set_ylim([min_scale, max_scale])
        axes.ticklabel_format(axis="y", style="sci", useOffset=False, useMathText=True)

    plt.tight_layout()
    outfile = os.path.join(inputs.output_dir, f"pmem{pmem.ring_loc}_energies.png")
    plt.savefig(outfile, dpi=300)
    plt.close()


def energy_rmsd_scatter(pmem, inunits="Ha", outunits="kJmol-1"):
    """Create a scatterplot of energies and rmsd values for each conformer pair using matplotlib.

    Parameters
    ----------
    pmem : kaplan.pmem.Pmem object
        Containts energy and rmsd values to plot.
    inunits : str
        The input units for the energy.
    outunits : str
        The output units for energy to display on the plot.

    Returns
    -------
    None

    """
    if pmem.num_geoms == 1:
        raise InputError("Need at least two geoms per pmem to make this plot.")
    # convert energy units, otherwise raise error
    try:
        assert inunits in energy_units
        assert outunits in energy_units[inunits]
    except AssertionError:
        raise ValueError(f"Unavailable unit conversion:\n{inunits} to {outunits}")
    try:
        assert all(x is not None for x in pmem.energies)
    except AssertionError:
        raise InputError("Pmem energies not set.")
    rmsds = []
    delta_es = []
    labels = []

    for conf1, conf2, rmsd in pmem.rmsds:
        e1 = energy_units[inunits][outunits](pmem.energies[conf1])
        e2 = energy_units[inunits][outunits](pmem.energies[conf2])
        delta_es.append(abs(e1 - e2))
        rmsds.append(rmsd)
        labels.append(f"{conf1},{conf2}")

    # regular size is 8x6
    # need at least 0.2 per pair
    # want a minimum of around 8 width
    if pmem.num_geoms > 8:
        fig, ax1 = plt.subplots(figsize=(pmem.num_pairs * 0.2, 6))
    else:
        fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    inputs = Inputs()
    ax1.set_xlabel("Conformer Pair", fontsize=20)
    ax1.set_ylabel("RMSD", fontsize=20)
    ax2.set_ylabel(r"$\Delta E_C$ " + f"{inunits}", fontsize=20)
    ax1.plot(labels, rmsds, 'g^', label=labels, markersize=10)
    for tick in ax1.get_xticklabels():
        tick.set_rotation(90)
    ax2.plot(labels, delta_es, 'ro', markersize=10)
    ax1.set_title(f"{inputs.struct_input}, Pmem{pmem.ring_loc}", fontsize=24)
    fig.legend(
        [Line2D([0], [0], marker="^", color="g"),
         Line2D([0], [0], marker="o", color="r")],
        ["RMSD", r"$\Delta E_C$"], loc="center"
    )
    outfile = os.path.join(inputs.output_dir, f"pmem{pmem.ring_loc}-delta-es-rmsds.png")
    plt.tight_layout()
    fig.savefig(outfile, dpi=300)
    plt.close()


def dihedrals_heatmap(ring):
    """Construct heatmap data containing R^2 values for each pair of dihedrals.

    Parameters
    ----------
    ring : kaplan.ring.Ring object
        Contains at least one pmem
        with dihedrals and fitness set.

    Notes
    -----
    Uses the linregress function from scipy.stats module.
    Makes a heatmap, where red represents the dihedrals
    that change together (high R^2) and blue represents
    the dihedrals that are not linearly correlated (low
    R^2). Sets of dihedrals are only considered if they
    can produce an energy (i.e. don't produce a convergence
    error). Use masked array from numpy to ignore None
    values.

    Returns
    -------
    r_sqr_matrix : np.array((inputs.num_dihed, inputs.num_dihed), float)
        The R^2 value for each pair of dihedral angles, as taken
        from doing a linear regression between all valid values
        for dihedral angle A and dihedral angle B, where
        0 <= A,B < inputs.num_dihed.

    """
    inputs = Inputs()
    dihedrals_matrix = np.zeros((inputs.num_dihed, ring.num_filled * inputs.num_geoms))
    dihedrals_mask = np.zeros(dihedrals_matrix.shape, int)
    for i in range(inputs.num_dihed):
        for j, slot in enumerate(ring.occupied):
            for k in range(inputs.num_geoms):
                # ignore dihedrals where energy was not calculable
                if ring[slot].energies[k] is not None:
                    dihedrals_matrix[i][j + k] = ring[slot].dihedrals[k][i]
                    dihedrals_mask[i][j + k] = 1

    masked_dihedrals = np.ma.MaskedArray(dihedrals_matrix, mask=dihedrals_mask)
    r_sqr_matrix = np.identity(inputs.num_dihed, float)

    for i in range(inputs.num_dihed):
        for j in range(i + 1, inputs.num_dihed):
            # results is a tuple of:
            # slope, intercept, r_value, p_value, std_err
            results = linregress(masked_dihedrals[i], masked_dihedrals[j])
            r_sqr = results[2] ** 2
            # upper right matrix
            r_sqr_matrix[i][j] = r_sqr
            # lower left matrix
            r_sqr_matrix[j][i] = r_sqr

    make_heatmap(r_sqr_matrix)

    return r_sqr_matrix


def make_heatmap(data):
    """Make a heatmap image from an input matrix.

    Parameters
    ----------
    data : np.array((N,N),float)
        The data to plot on the heatmap.

    Notes
    -----
    Saves a figure to the output directory called
    heatmap.png. This function is supposed to be called
    from dihedrals_heatmap, but can be used for
    other data (just change the axis labels).

    Returns
    -------
    None

    """
    # inputs is used for the outpur_dir location, nothing else
    inputs = Inputs()
    size = data.shape[0]
    # other colour maps can be chosen from here:
    # https://matplotlib.org/gallery/color/colormap_reference.html
    cmap_name = "coolwarm"
    # number of colour samples will change how many colours
    # there are on the colour bar (i.e. 128 would give a continuum,
    # whereas a small number gives colour blocks)
    num_colour_samples = 9

    mat = plt.matshow(data, cmap=cm.get_cmap(cmap_name, num_colour_samples),
                      interpolation="nearest", vmin=0, vmax=1)

    plt.colorbar(mat, fraction=0.046, pad=0.04)

    # set minor axes in between the labels
    ax = plt.gca()
    ax.set_xticks([x - 0.5 for x in range(1, size)], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, size)], minor=True)
    # plot grid on minor axes
    plt.grid(which="minor", ls="-", lw=1)
    # add axis labels
    ax.set_xlabel("Dihedral B")
    ax.xaxis.set_label_position("top")
    ax.set_ylabel("Dihedral A")
    # remove unnecessary ticks
    plt.tick_params(which="minor", top=False, left=False)
    plt.tick_params(which="both", bottom=False)

    # set x and y ticks and labels
    # row_labels = range(size)
    # col_labels = range(size)
    # plt.xticks(range(size), col_labels)
    # plt.yticks(range(size), row_labels)

    plt.savefig(os.path.join(inputs.output_dir, "heatmap.png"),
                bbox_inches="tight", dpi=300)
    plt.close()
