# nglutils
Utility module for using nglview to visualize stuff

# Examples
See the included jupyter notebook for example usage.

# Viewer commands
 * Left click to rotate
 * Right click to move
 * click an atom or bond to center view on it
 * Mouse wheel to zoom

# Installation
Clone the repo, then use pip or setup.py to install.

Note: it might be necessary to run the following two shell commands to enable
the NGL widget:
```sh
$jupyter-nbextension enable --py --user widgetsnbextension
$jupyter-nbextension enable --py --user nglview
```
