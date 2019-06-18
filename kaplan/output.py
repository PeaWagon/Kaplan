"""The purpose of this module is to handle writing output
after the conformer search has concluded.

This module decides where to put the output files, what
the output files will be called, and what the output
files contain.

The default output directory is:
os.cwd/kaplan_output/

Each time Kaplan is run, a new job number is created:
os.cwd/kaplan_output/job_0
os.cwd/kaplan_output/job_1
        ....
os.cwd/kaplan_output/job_n-1
where n is the number of times Kaplan is run.

Note: if run directories are deleted, Kaplan
will generate new job numbers depending on the
highest value it finds (so if job 0, 1, 5 and 10
are present, then the next job will be job_11).

Some new features that I'd like to add are
written as comments below."""

import os
import pickle

import numpy as np

from vetee.tools import periodic_table
from vetee.coordinates import write_xyz

from kaplan.inputs import Inputs
from kaplan.geometry import update_obmol

# OUTPUT_FORMAT = 'xyz'

# FEATURES TODO:
# add option to change output format
# make some images representing the population using matplotlib

units = {
    "psi4": "Ha",
    "openbabel": "kcal/mol",
}


def run_output(ring, save):
    """Run the output module.

    Parameters
    ----------
    ring : object
       The final ring data structure after evolution.
    save : bool
        If True, writes a pickle file of the ring and
        input objects.
        If False, does not write any pickle files.

    """
    inputs = Inputs()

    # find best pmem and its ring index
    best_pmem, best_fit = ring.best_pmem

    # write a stats file
    with open(os.path.join(inputs.output_dir, "stats-file.txt"), "w") as fout:        
        fout.write(f"average fitness:         {ring.mean_fitness} +/- {ring.std_fitness}\n")
        fout.write(f"median fitness:          {ring.median_fitness}\n")
        fout.write(f"best fitness:            {best_fit} (pmem{best_pmem})\n")
        fout.write(f"num conformers per pmem: {inputs.num_geoms}\n")
        fout.write(f"final slots filled:      {ring.num_filled}/{ring.num_slots}\n")
        fout.write(f"average energy ({units[inputs.prog]}): {ring.mean_energy} +/- {ring.std_energy}\n")
        fout.write(f"median energy ({units[inputs.prog]}):  {ring.median_energy}\n")
        fout.write(f"average rmsd:            {ring.mean_rmsd} +/- {ring.std_rmsd}\n")
        fout.write(f"median rmsd:             {ring.median_rmsd}\n")
        fout.write("--------------------------------------------------------\n")
        for conf in range(inputs.num_geoms):
            fout.write(f"Conformer {conf}:\n")
            fout.write(f"Dihedrals: {ring[best_pmem].dihedrals[conf]}\n")
            fout.write(f"Energy: {ring[best_pmem].energies[conf]}\n")
            fout.write("--------------------------------------------------------\n")
        fout.write(f"RMSDs:\n")
        for rmsd in ring[best_pmem].rmsds:
            fout.write(f"Conf{rmsd[0]} & Conf{rmsd[1]} => {rmsd[2]}\n")

    # generate the output file for the best pmem
    for i, geom in enumerate(ring[best_pmem].dihedrals):
        coords = update_obmol(inputs.obmol, inputs.min_diheds, geom)
        coords_with_atoms = []
        for j, atom in enumerate(coords):
            # need to convert from numpy's int64 to regular python int
            label = periodic_table(int(inputs.atomic_nums[j]))
            coords_with_atoms.append([label, atom[0], atom[1], atom[2]])
        comments = f"conformer {i}; energy {ring[best_pmem].energies[i]}"
        write_xyz({"_coords": coords_with_atoms, "_comments": comments},
                  os.path.join(inputs.output_dir, f"conf{i}.xyz"))
    
    # pickel if requested
    # note: cannot pickle obmol => TypeError: can't pickle SwigPyObject objects
    if save:
        with open(os.path.join(inputs.output_dir,"ring.pickle"), "wb") as f:
            pickle.dump(ring, f)
        with open(os.path.join(inputs.output_dir, "inputs.pickle"), "wb") as f:
            inputs.obmol = None
            pickle.dump(inputs, f)


"""Some basic functions for kaplan for post-processing."""

from vetee.coordinates import write_xyz
from vetee.tools import periodic_table
from kaplan.geometry import create_obmol

from kaplan.inputs import Inputs, InputError
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os
import pickle


# http://www.psicode.org/psi4manual/1.1/autodoc_physconst.html?highlight=bohr
# https://www.convertunits.com/from/hartree/to/joule
constants = {
    "avogadro's number": 6.02214179e23, # units mol-1
}
units = {
    
    "Ha": {                                      # hartrees
        "kJmol-1": lambda x: x*2625.500,         # kilojoules/mole
        "J": lambda x: x*4.359744E-18,           # joules
        "eV": lambda x: x*27.21138,              # electron-volts
        "kcal": lambda x: x*1.0415567394524e-21, # kilocalories
        "kcalmol-1": lambda x: x*627.5095,       # kilocalories/mole
        "cm-1": lambda x: x*219474.6,            # wavenumbers
        "MHz": lambda x: x*6.579684e9,           # megahertz
    },
    "kJmol-1": {
        "Ha": lambda x: x/2625.500,
        "kcalmol-1": lambda x: x/4.1858,
    },
    "J": {
        "Ha": lambda x: x/4.359744E-18,
    },
    "eV": {
        "Ha": lambda x: x/27.21138,
    },
    "kcal": {
        "Ha": lambda x: x/1.0415567394524e-21,
    },
    "kcalmol-1": {
        "Ha": lambda x: x/627.5095,
        "kJmol-1": lambda x: x*4.1858,
    },
    "cm-1": {
        "Ha": lambda x: x/219474.6,
    },
    "MHz": {
        "Ha": lambda x: x/6.579684e9,
    },
}

# usage example for units:
#result = units["Ha"]["kJmol-1"](1)    


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
    inputs = Inputs()
    # create bar labels
    labels = [f"conf{i}" for i in range(pmem.num_geoms)]
    # convert energy units, otherwise raise error
    try:
        energies = [units[inunits][outunits](e) for e in pmem.energies]
    except KeyError:
        raise ValueError(f"Unavailable unit conversion:\n{inunits} to {outunits}")
    # TypeError: unsupported operand type(s) for *: 'NoneType' and 'float'
    except TypeError:
        raise InputError("Pmem energies not set")
    x = range(pmem.num_geoms)
    #plt.bar(height=energies, x=x, label=labels)
    plt.bar(height=energies, x=x, tick_label=labels)
    plt.xlabel("Conformer Name")
    plt.ylabel("abs(energy) "+outunits)
    plt.title(f"Conformer Energies\n{inputs.struct_input}, Pmem{pmem.ring_loc}")
    if scale:
        plt.ylim(min(energies)*0.95, max(energies)*1.05)
    plt.tight_layout()
    outfile = os.path.join(inputs.output_dir, f"pmem{pmem.ring_loc}-energies.png")
    plt.savefig(outfile, dpi=300)

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
    inputs = Inputs()
    # convert energy units, otherwise raise error
    try:
        assert inunits in units
        assert outunits in units[inunits]
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
        e1 = units[inunits][outunits](pmem.energies[conf1])
        e2 = units[inunits][outunits](pmem.energies[conf2])
        delta_es.append(abs(e1-e2))
        rmsds.append(rmsd)
        labels.append(f"{conf1},{conf2}")

    #xt = range(len(rmsds))

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    ax1.set_xlabel("Conformer Pair", fontsize=20)
    ax1.set_ylabel("RMSD", fontsize=20)
    ax2.set_ylabel(r"$\Delta E_C$ "+f"{inunits}", fontsize=20)
    ax1.plot(labels, rmsds, 'g^', label=labels, markersize=10)
    ax2.plot(labels, delta_es, 'ro', markersize=10)
    ax1.set_title(f"{inputs.struct_input}, Pmem{pmem.ring_loc}", fontsize=24)
    fig.legend([Line2D([0], [0], marker="^", color="g"),
                Line2D([0], [0], marker="o", color="r")],
            ["RMSD", r"$\Delta E_C$"], loc="center")
    outfile = os.path.join(inputs.output_dir, f"pmem{pmem.ring_loc}-delta-es-rmsds.png")
    fig.savefig(outfile, dpi=300)



def apply_centroid_rotation(pmem):
    """Apply centering and rotating for input pmem.
    
    Notes
    -----
    This method should make the geometries for a pmem more
    comparable when viewing them in a visualiser (ex: VMD).
    They will be centred/rotated via rmsd module. Over-writes
    existing
    
    """
    pass