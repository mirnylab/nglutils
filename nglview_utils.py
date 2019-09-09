"""
Small custom utility functions, including plotting with nglview.
Use this with

    import sys
    sys.path.append('/path/to/this/file/')
    import nglview_utils as ngu
    
    # Usage: given a list of 3d coordinates named xyz:
    view = xyz2nglview(xyz)
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
# custom things
import openmmlib

# ----------- Handling of hexcolors -----------------

def color_as_hexstring(val, cmap=plt.cm.viridis, include_alpha=False):
    """ Maps a scalar in [0, 1] to a color and represents it as a string, e.g. 0xff440154 (for viridis(0)) """
    col = cmap(float(max(0, min(val, 1)))) # Clamping val to [0, 1]
    if include_alpha:
        return "0x{3:02x}{0:02x}{1:02x}{2:02x}".format(int(255*col[0]),
                                                       int(255*col[1]),
                                                       int(255*col[2]),
                                                       int(255*col[3]))
    else:
        return "0x{0:02x}{1:02x}{2:02x}".format(int(255*col[0]), int(255*col[1]), int(255*col[2]))
    
def hexcolors_polymers(lengths, maps=[plt.cm.viridis]):
    """
    Generate a list of colors for use in nglview.NGLWidget._set_color_by_residue
    The idea is to patch several colormaps (or the same one multiple times), one for each chain.
    
    Parameters
    ----------
    lengths : np.array
        lengths of the polymer chains (no. of residues in mdtraj language)
    maps : list of executable objects
        list of plt.cm colormaps to be cycled through
        Default: [plt.cm.viridis]
    """
    lengths = np.cumsum(np.insert(np.array(lengths), 0, 0))
    
    def chainNo(res):
        return np.min(np.argwhere(lengths > res))
    def pos_in_chain(res):
        chain = chainNo(res)
        return (res - lengths[chain-1])/(lengths[chain]-1 - lengths[chain-1])
    def cmap(res):
        return maps[(chainNo(res)-1) % len(maps)]
    
    return [color_as_hexstring(pos_in_chain(res), cmap=cmap(res)) for res in range(lengths[-1])]

# ------------ Loading stuff from disk ----------------

def load_traj_from_path(path, mon):
    """
    Load trajectory of monomer mon into an array.
    The folder in path is assumed to have the openmmlib generated structure of block*.dat
    """
    filenames = openmmlib.polymerutils.scanBlocks(path)["files"]

    # Extract trajectory of the pulled monomer
    traj = np.zeros((len(filenames), 3))
    for i in range(len(filenames)):
        traj[i, :] = openmmlib.polymerutils.load(filenames[i])[mon, :]
    
    return traj

# ------------ Displaying polymers with nglview -------------

def mdtop_for_polymer(N):
    """
    Generate an mdtraj.Topology object for a single polymer of length N.
    """
    top = mdtraj.Topology()
    ch = top.add_chain()
    res = top.add_residue("res", ch)
    top.add_atom('CA', mdtraj.element.carbon, res)
    for i in range(1, N):
        res = top.add_residue("res", ch)
        atom = top.add_atom('CA', mdtraj.element.carbon, res)
        top.add_bond(atom, top.atom(i-1))
    return top

def xyz2mdtraj(xyz_array, top=None):
    """
    Generate an mdtraj.Trajectory from a numpy array with coordinates.
    This is done via transient writing of an .xyz file.
    
    Parameters
    ----------
    xyz_array : np.array
        either (n_atoms, 3) or (n_frames, n_atoms, 3). The coordinates of the atoms.
    top : mdtraj.Topology
        topology of the chain. If omitted, a single long chain is assumed.
    """
    if len(xyz_array.shape) == 2:
        xyz_array = np.expand_dims(xyz_array, 0)
    if len(xyz_array.shape) != 3:
        error('xyz_array has to have shape (n_atoms, 3) or (n_frames, n_atoms, 3)!')
    
    if top is None:
        top = mdtop_for_polymer(xyz_array.shape[1])
     
    with tempfile.TemporaryDirectory() as tmpdir:
        filename = os.path.join(tmpdir, 'poly.xyz')
        with mdtraj.formats.XYZTrajectoryFile(filename, mode='w') as file:
            file.write(xyz_array)
        traj = mdtraj.load(filename, top=top)
    
    return traj

def mdtraj2nglview(mdtrajs, cmaps=[plt.cm.viridis]):
    """
    Displays one or several mdtraj trajectories together, using the specified color maps
    
    Parameters
    ----------
    mdtrajs : list of mdtraj.Trajectory
    cmaps : list of colormaps
    """
    traj = mdtrajs[0]
    for i in range(1, len(mdtrajs)):
        traj = traj.stack(mdtrajs[i])
    
    view = nglview.show_mdtraj(traj)
    view.clear_representations()
    view.add_representation('licorice')
    view._set_color_by_residue(hexcolors_polymers([t.n_residues for t in mdtrajs], cmaps))
    
    return view

def xyz2nglview(xyz_array, cmaps=[plt.cm.viridis], top=None):
    """ 
    Show an array containing xyz coordinates with nglview.
    This is just a quick shortcut for the two commands shown.
    If xyz_array is a list of such arrays, it is assumed that top is a list of the same length (or None).
    """
    if not isinstance(xyz_array, list):
        xyz_array = [xyz_array]
        
    if top is None:
        traj = [xyz2mdtraj(xyz_array[i]) for i in range(len(xyz_array))]
    else:
        if not isinstance(top, list):
            top = [top]
        traj = [xyz2mdtraj(xyz_array[i], top[i]) for i in range(len(xyz_array))]
        
    view = mdtraj2nglview(traj, cmaps)
    
    return view
