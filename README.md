# nglutils
Utility module for using nglview to visualize stuff

# Examples
See the included jupyter notebook for example usage.

# Viewer commands
 * Left click to rotate
 * Right click to move
 * click an atom or bond to center view on it
 * Mouse wheel to zoom

# Compartments
To visualize compartments, one could split the polymer into many small chains
(one for each compartment). This has *heavy* memory overhead, so should be
avoided. Instead, one can join all the monomers of one type into a long chain
and explicitly exclude the bonds between the different segments. This works way
better.
