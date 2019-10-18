"""
A collection of helper functions that come in handy when working with NGLViewer

Note: may have to run the following commands in shell
$jupyter-nbextension enable --py --user widgetsnbextension
$jupyter-nbextension enable --py --user nglview
"""

# System stuff
import os
import sys
import tempfile
# Python stuff
import matplotlib.pyplot as plt
import numpy as np
# specific things
import mdtraj
import nglview

from . import rep_add

# ----------- Utilities -----------------

def intlist_to_alpha(intlist):
    """
    Convert a list of integers [0, 0, 1, 1, 5, ...] to a list of letters [A, A,
    B, B, F, ...]. Useful for converting monomer types to compartment names
    (which can then be used as atom names).

    Parameters
    ----------
    intlist : list of int
        list to be converted

    Notes
    -----
    If intlist has values >= 27, you will get the corresponding Unicode
    literals.
    """
    return [chr(65+i) for i in intlist]

# ------------ Displaying polymers with nglview -------------

def mdtop_for_polymer(N, exclude_bonds=np.array([]), add_bonds=np.array([]), chains=[(0, None, False)], atom_names="XXX"):
    """
    Generate an mdtraj.Topology object for a single polymer of length N.

    Parameters
    ----------
    N : int
        total number of monomers in the polymer, only mandatory argument.
    chains : list of 3-tuples (polychrom chain specification)
        chain specification, in the format (first_mon, last_mon, ring?), just
        as polychrom does. Note that every atom has to belong to exactly one
        chain, thus the chains have to be contiguous.
    exclude_bonds : array-like, (any)x1
        list of chain bonds to exclude from the topology. The bonds are
        numbered with the second atom in it, i.e. the bond between atom 59 and
        atom 60 has the number 60.
    add_bonds : array-like, (any)x2
        list of additional bonds to include. Every entry is read as two atoms
        that are to be bound together.
    atom_names : string (up to 3 characters) or list of such
        names for the atoms. This can be useful for labelling different types
        of monomers, such as compartments. If given as list, should have length
        N.

    Notes
    -----
        - use iron as element to prevent NGLView from calculating bonds
        - also don't try to make NGLView detect this as backbone, because then
          it would try to calculate bonds again (even though this would be nice
          for some representations)
        - we can exploit the way NGLView displays the residue name and number
          to give a seven digit atom number
    """
    # Check that chains are good
    npchains = np.sort(np.array(chains), axis=0)
    if not all(npchains[:, 0] == [0, *npchains[:-1, 1]]):
        print("WARNING: chains should be contiguous! Using one long chain.")
        npchains = np.array([(0, N, 0)])
    if npchains[-1, 1] is None:
        npchains[-1, 1] = N

    # Check that atom names are good
    if isinstance(atom_names, str):
        atom_names = [atom_names[:3] for _ in range(N)]
    elif not isinstance(atom_names, list) or not len(atom_names) == N:
        print("WARNING: atom_names should be a list of length N = {}. Using default." % N)
        atom_names = ["XXX" for _ in range(N)]

    # Generate topology
    top = mdtraj.Topology()
    resi = 0
    for chspec in npchains:
        ch = top.add_chain()
        res = top.add_residue("{0:03d}".format(int(np.floor(resi/10000))), ch)
        first_atom = top.add_atom(atom_names[resi], mdtraj.element.iron, res)
        resi += 1
        prev_atom = first_atom
        for _ in range(chspec[0]+1, chspec[1]):
            res = top.add_residue("{0:03d}".format(int(np.floor(resi/10000))), ch)
            cur_atom = top.add_atom(atom_names[resi], mdtraj.element.iron, res)
            resi += 1
            if not resi in exclude_bonds:
                top.add_bond(cur_atom, prev_atom)
            prev_atom = cur_atom
        if chspec[2]:
            top.add_bond(cur_atom, first_atom)

    for bond in add_bonds:
        top.add_bond(top.atom(bond[0]), top.atom(bond[1]))

    return top

def xyz2mdtraj(xyz_array, top=None):
    """
    Generate an mdtraj.Trajectory from a numpy array with coordinates.  This is
    done via transient writing of an .xyz file.
    
    Parameters
    ----------
    xyz_array : np.array
        either (n_atoms, 3) or (n_frames, n_atoms, 3). The coordinates of the
        atoms.
    top : mdtraj.Topology
        topology of the chain. If omitted, a single long chain is assumed.

    Notes
    -----
    This used to take an excluded_bonds argument as well. This use is
    discouraged now. Generate the topology yourself, using mdtop_for_polymer.
    """
    if len(xyz_array.shape) == 2:
        xyz_array = np.expand_dims(xyz_array, 0)
    if len(xyz_array.shape) != 3:
        error('xyz_array has to have shape (n_atoms, 3) or (n_frames, n_atoms, 3)!')
    
    if top is None:
        top = mdtop_for_polymer(xyz_array.shape[1])
     
    # use tempdir instead of file, because we want mdtraj to create the file
    with tempfile.TemporaryDirectory() as tmpdir:
        filename = os.path.join(tmpdir, 'poly.xyz')
        with mdtraj.formats.XYZTrajectoryFile(filename, mode='w') as file:
            file.write(xyz_array)
        traj = mdtraj.load(filename, top=top)
    
    return traj

def mdtraj2nglview(mdtrajs, cmaps=[plt.cm.viridis]):
    """
    Displays one or several mdtraj trajectories together, using the specified
    color maps
    
    Parameters
    ----------
    mdtrajs : list of mdtraj.Trajectory
    cmaps : list of colormaps
    """
    if not isinstance(mdtrajs, list):
        mdtrajs = [mdtrajs]
    if not isinstance(cmaps, list):
        cmaps = [cmaps]

    #traj = mdtrajs[0]
    #for i in range(1, len(mdtrajs)):
    #    traj = traj.stack(mdtrajs[i])
    
    #view = nglview.show_mdtraj(traj)
    view = nglview.NGLWidget()
    for i, traj in enumerate(mdtrajs):
        view.add_trajectory(traj)
        view.clear_representations(component=i)
        rep_add.colormap(view, cmaps[i%len(cmaps)], component=i)
    
    return view

def xyz2nglview(xyz_array, cmaps=plt.cm.viridis, top=None):
    """ 
    Display a trajectory / trajectories with nglview.

    Parameters
    ----------
    xyz_array : np.array of shape (any, any, 3) or (any, 3), or list of such
        the trajectory data. If given as a list, the single arrays will be
        added to the view as different components. If the array is
        three-dimensional, the first dimension is assumed to be time.
    cmaps : function, or list of such
        colormap(s) to use. Will be cycled through if shorter than number of
        trajectories. Should be a function taking a scalar value between 0 and
        1 and returning an rgb triplet, such as for example plt.cm.viridis
    top : mdtraj.Topology
        topology of the trajectories. mdtop_for_polymer comes in handy here.
    """
    if not isinstance(xyz_array, list):
        xyz_array = [xyz_array]
    if not isinstance(cmaps, list):
        cmaps = [cmaps]
        
    if top is None:
        traj = [xyz2mdtraj(xyz_array[i]) for i in range(len(xyz_array))]
    else:
        if not isinstance(top, list):
            top = [top]
        if not len(top) == len(xyz_array):
            raise ValueError("Need a topology for every trajectory!")

        traj = [xyz2mdtraj(xyz_array[i], top[i]) for i in range(len(xyz_array))]
        
    view = mdtraj2nglview(traj, cmaps)
    
    return view
