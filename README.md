# ABCD 2-port functions for phasing vertical antenna arrays. 

If you don't know what ABCD 2-port functions are,  see
https://en.wikipedia.org/wiki/Two-port_network#ABCD-parameters.

The Jupyter notebook examples in the repo are taken from Gehrke K2BT's "Vertical Phased Arrays"
series in Ham Radio Magazine, May-July 1983, October 1983, December 1983,
and from Capt. Paul Lee N6PL's "The Amateur Radio Vertical Antenna Handbook".
Most of the equations used for the library were taken from the K2BT article.
The stub and quarter line solution functions were taken from
Sophocles J. Orfanidis's
[Electromagnetic Waves and Antennas](http://eceweb1.rutgers.edu/~orfanidi/ewa/).

Below are the list of functions provided by the library.
The first group of functions calculate the ABCD matrix itself.
The impedance of the components that make up the
matrix are passed to these functions in signal flow order, that
is starting from the source side to the load side.  So 
for a shunt input L match, the impedance of the shunt component 
would be the first argument and the impedance of the series component
would be the second.

The next group of functions calculates solutions for various impedance 
matching problems.  Their results are returned as a list.  This list can
then be passed directly using the Python star notation to the appropriate
ABCD 2-port function.  So, for example, the result of the input shunt 
L-match solver function can be passed as is using the 'star' to the input L-match
ABCD matrix function.
 
The other groups of functions are helper functions, functions that work on ABCD vectors,
and functions for printing results:

```
# generates ABCD matrix
######################################
def auto(ratio, xt, k=1, q=0): # an 1:n autotransformer
def mutual(n, x1, k=1, q=0):   # a 1:n transformer
def trans(n):                  # an ideal 1:n transformer
def series(z):                 # a component in series
def shunt(y):                  # a component in parallel
def halfpi(y, z):              # a shunt input l-match
def halftee(z, y):             # a series input l-match
def tline(deg, zo=50, loss=0): # a transmission line of length deg, db loss
def fulltee(z1, z2, z3):       # a tee section
def fullpi(z1, z2, z3):        # a pi section

# solvers whose results are passed into the above ABCD functions
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
def to_lmin(za, zo=50):      # distance to voltage min/max
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
def opened_stub(deg, zo=50):
def shorted_stub(deg, zo=50):
def component_value(impedance, fd):
def qlosses(impedances, q=200):
def parallel(*impedances):
def z2g(z, zo=50):
def g2z(gm, zo=50):
def swr(gm):
def s2p(z):                  # series to parallel

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

