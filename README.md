# ABCD 2-port functions for phasing vertical antenna arrays. 

The examples are taken from Gehrke K2BT's "Vertical Phased Arrays"
series in Ham Radio Magazine, May-July 1983, October 1983, December 1983.

The first group of functions calculate the ABCD matrix.
The order of the components are entered in order from the
source side to the load side.

```
auto(ratio, xt, k=1):       an autotransformer; 1:n is the ratio
mutual(n, x1, k=1):         a transformer; 1:n is the ratio
trans(n):                   an ideal transformer; 1:n is the ratio
series(z):                  a component in series
shunt(y):                   a component in parallel
halfpi(y, z):               a shunt input l-match
halftee(z, y):              a series input l-match
tline(deg, zo=50):          a transmission line z of length deg
fulltee(z1, z2, z3):        a tee section
fullpi(z1, z2, z3):         a pi section
halftee2(z1, y1, z2, y2):   a double series input l-match
halfpi2(z1, y1, z2, y2):    a double shunt input l-match
```

The second group of functions solve various impedance matching problems.
The solution parameter changes the components to their compliments.

```
to_halfpi(rin, za, solution=0):     match with a shunt input l-match
to_halftee(rin, za, solution=0):    match with a series input l-match
to_halfwave(zs, za, solution=0):    match with a 90 degree tee/pi section
to_shunt(za):                       cancel reactance with a shunt section
to_series(za):                      cancel reactance with a series section
to_fullpi(deg, zo):                 shift phase with a pi section
to_fulltee(deg, zo):                      shift phase with a tee section
to_halfpi2(zs, za, solution=(0,0)):       match with a double shunt l-match
to_halftee2(zs, za, solution=(0,0)):      match with a double series l-match
to_fullpi2(zs, za, q=0, solution=(0,0)):  match with a pi section
to_fulltee2(zs, za, q=0, solution=(0,0)): match with a tee section
```

The third group of functions are basically helper functions

```
qmin(zs, za):
qmin2(zs, za):
openstub(deg, zo):
shortstub(deg, zo):
reactance_value(component, fd):
component_value(impedance, fd):
parallel(*impedances):
vec(e, i):
emag(v):
ephase(v):
power(*vectors):
impedance(*vectors):
```

The fourth group of functions are for printing results:

```
polar(x):
status(v, note=0):
component(values, fd, precision=4):
notation(x, precision=4):
```

# Example

For example, with a load of 1 amp and 5 volts flowing through it, transform
it with 90 degrees of 50 ohm transmission line.


```
from abcd import *
line = vec(50, 1)
line = tline(135) * line
status(line)
```

Gives

```
p(0) = 50.0000
z(0) = 50.0000+0.0000j
i(0) = 1.0000 / 135.0000
e(0) = 50.0000 / 135.0000
```






