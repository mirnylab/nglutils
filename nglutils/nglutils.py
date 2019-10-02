"""
Small custom utility functions, including plotting with nglview.
Use this with

    import sys
    sys.path.append('/path/to/this/file/')
    import nglutils as ngu
    
    # Usage: given a list of 3d coordinates named xyz:
    view = ngu.xyz2nglview(xyz)
    view

Note: use nglview-2.1.0 for this to work. Also, maybe have to
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

# ----------- Handling of hexcolors -----------------

def hexscale_from_cmap(cmap, N):
    """
    Evaluate a colormap at N points.

    Parameters
    ----------
    cmap : function
        a function taking a scalar value between 0 and 1 and giving a color as
        rgb(a) with values between 0 and 1. These are for example the pyplot
        colormaps, like plt.cm.viridis
    N : int
        number of steps on which to evaluate the colormap

    Returns
    -------
    scale : a list of numbers representing the colors from the map, written in
        the format 0xrrggbb
    """
    rgb = [(round(255*col[0]), round(255*col[1]), round(255*col[2])) for col in map(cmap, np.arange(N)/(N-1))]
    return [0x010000*col[0] + 0x000100*col[1] + 0x000001*col[2] for col in rgb]

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
    chains.sort()
    chains = np.array(chains)
    if not all(chains[:, 0] == [0, *chains[:-1, 1]]):
        print("WARNING: chains should be contiguous! Using one long chain.")
        chains = np.array([(0, N, 0)])
    if chains[-1, 1] is None:
        chains[-1, 1] = N

    # Check that atom names are good
    if isinstance(atom_names, str):
        atom_names = [atom_names[:3] for _ in range(N)]
    else if not isinstance(atom_names, list) or not len(atom_names) == N:
        print("WARNING: atom_names should be a list of length N = {}. Using default." % N)
        atom_names = ["XXX" for _ in range(N)]

    # Generate topology
    top = mdtraj.Topology()
    resi = 0
    for chspec in chains:
        ch = top.add_chain()
        res = top.add_residue("{0:03d}".format(int(np.floor(resi/1000))), ch)
        first_atom = top.add_atom(atom_names[resi], mdtraj.element.iron, res)
        resi += 1
        prev_atom = first_atom
        for _ in range(chspec[0]+1, chspec[1]):
            res = top.add_residue("{0:03d}".format(int(np.floor(resi/1000))), ch)
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
    if not isinstance(cmaps, list):
        cmaps = [cmaps]

    traj = mdtrajs[0]
    for i in range(1, len(mdtrajs)):
        traj = traj.stack(mdtrajs[i])
    
    view = nglview.show_mdtraj(traj)
    view.clear_representations()
    for i in range(len(mdtrajs)):
        # add representations for all trajs
        view.add_representation('licorice', selection=':{0}'.format(chr(65+i)),
                                            colorScheme='atomindex',
                                            colorScale=hexscale_from_cmap(cmaps[i], mdtrajs[i].n_atoms))
    
    return view

def xyz2nglview(xyz_array, cmaps=plt.cm.viridis, top=None, scale=1, exclude_bonds=[np.array([])]):
    """ 
    Show an array containing xyz coordinates with nglview.  This is just a
    quick shortcut for the two commands shown.  If xyz_array is a list of such
    arrays, it is assumed that top is a list of the same length (or None).
    Optional arguments:
        - scale = 5: scale up all coordinates. Used to prevent nglview from
                     introducing unwanted bonds between nearby atoms. Also
                     usable to adjust relative thickness of the bonds
        - exclude_bonds = np.array([5, 65, 185]): bonds to exclude from
                     topology. See notes at 'mdtop_for_polymer'.  Should be a
                     list, if more than one chain is provided.  Ignored if
                     topology is given explicitly.
    """
    if not isinstance(xyz_array, list):
        xyz_array = [xyz_array]
        if not isinstance(cmaps, list):
            cmaps = [cmaps]
        if not isinstance(exclude_bonds, list):
            exclude_bonds = [exclude_bonds]
    else: # Check that exclude_bonds are compatible with list.
          # cmaps are automatically tiled, so that's okay
        if not isinstance(exclude_bonds, list):
            raise ValueError
        if len(exclude_bonds) == 1 and exclude_bonds[0].size == 0:
            exclude_bonds = [np.array([]) for jj in range(len(xyz_array))]
        
    if top is None:
        traj = [xyz2mdtraj(scale*xyz_array[i], exclude_bonds=exclude_bonds[i]) for i in range(len(xyz_array))]
    else:
        if not isinstance(top, list):
            top = [top]
        traj = [xyz2mdtraj(scale*xyz_array[i], top[i]) for i in range(len(xyz_array))]
        
    view = mdtraj2nglview(traj, cmaps)
    
    return view
