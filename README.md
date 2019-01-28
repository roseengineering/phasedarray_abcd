# ABCD 2-port functions for phasing vertical antenna arrays. 

The examples are taken from Gehrke K2BT's "Vertical Phased Arrays"
series in Ham Radio Magazine, May-July 1983, October 1983, December 1983,
and from Capt. Paul Lee N6PL's "The Amateur Radio Vertical Antenna Handbook".

What is a ABCD 2-port matrix?  See https://en.wikipedia.org/wiki/Two-port_network#ABCD-parameters 

The one group of functions calculate the ABCD matrix.
The impedance of the components are passed in order from the
source side to the load side.
Another group of functions solve various impedance matching problems.
The solution parameter changes the components to their compliments.
The other groups are helper functions, functions that work on ABCD vectors,
and functions for printing results:

```
# generates ABCD matrix
######################################
def auto(ratio, xt, k=1):    # an 1:n autotransformer
def mutual(n, x1, k=1):      # a 1:n transformer
def trans(n):                # an ideal 1:n transformer
def series(z):               # a component in series
def shunt(y):                # a component in parallel
def halfpi(y, z):            # a shunt input l-match
def halftee(z, y):           # a series input l-match
def tline(deg, zo=50, loss=0): # a transmission line of length deg, db loss
def fulltee(z1, z2, z3):     # a tee section
def fullpi(z1, z2, z3):      # a pi section

# solvers whose results are passed into above ABCD functions
#############################################################
def to_halfwave(zs, za):     # match with a 90 degree tee/pi section
def to_halfpi(rin, za):      # match with a shunt input l net, rin > za.real
def to_halftee(rin, za):     # match with a series input l net, rin < za.real
def to_fullpi(deg, zo):      # shift phase with a pi section
def to_fulltee(deg, zo):     # shift phase with a tee section
def to_shunt(za):            # cancel reactance with a shunt section
def to_series(za):           # cancel reactance with a series section

# experimental
########################################
def to_stub1(za, zo=50, shorted=True): # match with a stub-series input 
def to_qwt1(za, zo=50):
def to_qwt2(za, zo=50):
def to_qwt3(za, z2, zo=50):

# ABCD vector functions
#####################################
def vec(e, i):               # returns a ABCD vector for E, I
def emag(v):                 # returns the magnitude of E
def ephase(v):               # returns the phase of E
def power(*vs):              # power of one or more lines together
def impedance(*vs):          # impedance of one or more lines together
def emax(v, zo=50):          # maximum rms voltage on transmission line

# helper functions
######################################
def qmin(zs, za):            # q for L network, or minimum q for tee or pi
def qmin2(zs, za):           # minimum q for a LL network
def opened_stub(deg, zo=50):
def shorted_stub(deg, zo=50):
def component_value(impedance, fd):
def qlosses(*impedances, q=200):
def parallel(*impedances):
def z2g(z, zo=50):
def g2z(gm, zo=50):
def swr(gm):
def s2p(z):                  # serial to parallel
def lmin(za, zo=50):         # distance to voltage min/max

# print functions
####################################
def polar(x):
def status(v, note=0):
def component(values, fd, precision=4):
def notation(x, precision=4, units=None):
```

# Example

For example, with a load of 1 amp and 5 volts flowing through it, transform
the load with 135 degrees of 50 ohm transmission line to get what the source
sees.


```
from abcd import *
line = vec(5, 1)
status(line)
line = tline(135) * line
status(line)
```

Gives

```
p(0) = 5.0000
z(0) = 5.0000+0.0000j
i(0) = 1.0000 / 0.0000
e(0) = 5.0000 / 0.0000

p(0) = 5.0000
z(0) = 9.9010-49.0099j
i(0) = 0.7106 / 174.2894
e(0) = 35.5317 / 95.7106
```

