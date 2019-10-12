# Provide some nice colors to use in plots
# By 'nice' I mean the solarized colorscheme:
#
# https://ethanschoonover.com/solarized/
#
# NOTE: the dark background is 'base03', the light one is 'base3'. To set the
# background in the NGLWidget, use
# ```py
# view.stage.set_parameters(background_color=nglutils.colors.base3.ashex())
# ```

class color():
    def __init__(self, *args):
        if len(args) == 1: # It's a hex
            RGBhex = int(args[0])
            self._R = RGBhex >> 16
            self._G = (RGBhex % 0x010000) >> 8
            self._B = RGBhex % 0x000100
        elif len(args) == 3: # RGB, either 0-1 or 0-255
            if max(args) <= 1:
                self._R = int(255*args[0])
                self._G = int(255*args[1])
                self._B = int(255*args[2])
            else:
                self._R = int(args[0])
                self._G = int(args[1])
                self._B = int(args[2])
        else:
            raise ValueError("Could not interpret color: {}".format(args))

    def ashex(self):
        return 0x010000*self._R + 0x0000100*self._G + 0x000001*self._B

    def asrgb(self):
        return [self._R/255., self._G/255., self._B/255]

    def asRGB(self):
        return [self._R, self._G, self._B]

    def ashexstring(self):
        return "#{0:02x}{1:02x}{2:02x}".format(self._R, self._G, self._B)

# And now, generate some good colors
base03 = color(0x002b36)
base02 = color(0x073642)
base01 = color(0x586e75)
base00 = color(0x657b83)
base0 = color(0x839496)
base1 = color(0x93a1a1)
base2 = color(0xeee8d5)
base3 = color(0xfdf6e3)
yellow = color(0xb58900)
orange = color(0xcb4b16)
red = color(0xdc322f)
magenta = color(0xd33682)
violet = color(0x6c71c4)
blue = color(0x268bd2)
cyan = color(0x2aa198)
green = color(0x859900)
