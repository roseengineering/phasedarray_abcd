# ABCD 2-port functions for phasing vertical antenna arrays. 

The examples are taken from Gehrke K2BT's "Vertical Phased Arrays"
series in Ham Radio Magazine, May-July 1983, October 1983, December 1983,
and from Capt. Paul Lee N6PL's "The Amateur Radio Vertical Antenna Handbook".

What is a ABCD 2-port matrix?  See https://en.wikipedia.org/wiki/Two-port_network#ABCD-parameters 

The first group of functions calculate the ABCD matrix.
The impedance of the components are passed in order from the
source side to the load side.

```
# generates ABCD matrix
######################################
def auto(ratio, xt, k=1):       an autotransformer; 1:n is the ratio
def mutual(n, x1, k=1):         a transformer; 1:n is the ratio
def trans(n):                   an ideal transformer; 1:n is the ratio
def series(z):                  a component in series
def shunt(y):                   a component in parallel
def halfpi(y, z):               a shunt input l-match
def halftee(z, y):              a series input l-match
def tline(deg, zo=50, loss=0):  a transmission line z of length deg, loss in db
def fulltee(z1, z2, z3):        a tee section
def fullpi(z1, z2, z3):         a pi section
```

The second group of functions solve various impedance matching problems.
The solution parameter changes the components to their compliments.

```
# solvers whose results are passed into above ABCD functions
#############################################################
def to_halfwave(zs, za):        match with a 90 degree tee/pi section
def to_halfpi(rin, za):         match with a shunt input l-match
def to_halftee(rin, za):        match with a series input l-match
def to_fullpi(deg, zo):         shift phase with a pi section
def to_fulltee(deg, zo):        shift phase with a tee section
def to_shunt(za):               cancel reactance with a shunt section
def to_series(za):              cancel reactance with a series section
def to_stub1(za, zo=50, shorted=True):  match with a stub-series input 
```

The third group of functions are basically helper functions

```
# helper functions
######################################
def qmin(zs, za):
def qmin2(zs, za):
def open_stub(deg, zo=50):
def shorted_stub(deg, zo=50):
def reactance_value(component, fd):
def component_value(impedance, fd):
def qlosses(*impedances, q=200):
def parallel(*impedances):
def z2g(z, zo=50):
def g2z(gm, zo=50):
def swr(gm):
def s2p(z):
```

The fourth group of functions work on ABCD vectors:

```
# ABCD vector functions
#####################################
def vec(e, i):                    returns a ABCD vector for E, I
def emag(v):                      returns the magnitude of E
def ephase(v):                    returns the phase of E
def power(*vs):                   power of one or more lines together
def impedance(*vs):               impedance of one or more lines together
def emax(v, zo=50):               maximum rms voltage on transmission line
```

The last group of functions are for printing results:

```
# print functions
####################################
def polar(x):
def status(v, note=0):
def reactance(values, fd, precision=4):
def component(values, fd, precision=4):
def notation(x, precision=4):
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

