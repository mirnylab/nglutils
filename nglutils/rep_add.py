import numpy as np

def hexscale_from_cmap_subchain(cmap, ind):
    """
    NOTE: this is pretty bad, conceptually...

    Same as hexscale_from_cmap, but patch left and right ends. This is used to
    color only a substretch of polymer.

    Parameters:
    -----------
    cmap : see hexscale_from_cmap. Example: plt.cm.viridis
    ind : tuple of int: (start, end, N)
        the indices of the stretch to evaluate the colormap for, together with
        the total number of monomers, such that we know how much to patch
    """

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

def uniform(view, color, selection='all', component=0, **kwargs):
    """
    A shortcut to add a uniformly colored representation to the view.

    Parameters
    ----------
    view : NGLWidget
        the view in question
    selection : string
        a selection string. For the format, see here:
        http://nglviewer.org/ngl/api/manual/selection-language.html
    color : int
        a color value, in the format 0xrrggbb
    """
    view.add_representation('licorice', component=component,
                                        selection=selection,
                                        colorScheme='uniform',
                                        colorValue=color,
                                        **kwargs)

def colormap(view, cmap, selection='all', component=0, Nscale=100, **kwargs):
    """
    A shortcut to color a polymer with a given colormap.

    Parameters
    ----------
    view : NGLWidget
        the view in question
    cmap : function or list of colors (either hex values or strings)
        a function taking a scalar value and giving an rgb color, such as for
        example plt.cm.viridis
        OR
        a list of colors, which will be handed to the colorScale argument of
        the representation. Examples: ["Red", "Yellow", "Green"]
                                      [0xff0000, 0xffff00, 0x0000ff]
    Nscale : int
        number of discretization steps for the colormap, if given as function
    selection : string
        selection string. For the format, see
        http://nglviewer.org/ngl/api/manual/selection-language.html
    component : int
        the component of the view we are talking about (usually there's only
        one)
    """
    if isinstance(cmap, list) or isinstance(cmap, str):
        cScale = cmap
    else:
        cScale = hexscale_from_cmap(cmap, Nscale)

    view.add_representation('licorice', component=component,
                                        selection=selection,
                                        colorScheme='residueindex',
                                        colorScale=cScale,
                                        **kwargs)

def colormap_subchain(view, cmap, ind, N, component=0, Nscale=None, **kwargs):
    """
    Colormap for subchain, which is annoying, because we can only color by
    residue number in the whole chain, so we have to patch the colormap.

    Parameters
    ----------
    ind : 2-tuple of int
        (start, end) of the subchain to color (open interval, as always in python)
    N : int
        total number of monomers in the chain (i.e. how much to patch
    for everything else see docstring of rep_add.colormap()
    """
    l_subchain = ind[1] - ind[0]
    if Nscale is None:
        Nscale = min(100, l_subchain // 10 + 1)

    if isinstance(cmap, list):
        Nscale = len(cmap)
    else:
        cmap = hexscale_from_cmap(cmap, Nscale)

    # Figuring out this paragraph was a bit of a pain, don't touch it if you
    # don't know exactly what you're doing
    kcur = N*(Nscale-1) // l_subchain
    nmax = lambda k: ind[0]*k // N
    nmin = lambda k: np.ceil((ind[1]-1)*k / N - Nscale + 1).astype(int)
    while nmax(kcur) < nmin(kcur):
        kcur -= 1

    n = nmax(kcur)
    m = max(0, kcur - n - Nscale + 1)

    cScale = n*[0x000000] + cmap + m*[0x000000]
    view.add_licorice('@'+','.join(np.arange(*ind).astype(str)),
                      component=component,
                      colorScheme='residueindex',
                      colorScale=cScale,
                      **kwargs)
