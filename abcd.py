
import numpy as np


# generates ABCD matrix
######################################

def auto(ratio, xt, k=1):    # an autotransformer; 1:n is the ratio
    n = ratio / (1 - ratio)
    x1 = xt / (1 + n**2 + 2 * k * n)
    x2 = x1 * n**2
    xm = k * n * x1
    return fulltee((x1 + xm) * 1j, (x2 + xm) * 1j, -xm * 1j)

def mutual(n, x1, k=1):      # a transformer; 1:n is the ratio
    x2 = x1 * n**2
    xm = k * n * x1
    return fulltee((x1 - xm) * 1j, xm * 1j, (x2 - xm) * 1j)

def trans(n):                # an ideal transformer; 1:n is the ratio
    return np.matrix([[ 1/n, 0], [0, n]])
    
def series(z):               # a component in series
    """
    <---Z---< ZA
    """
    return np.matrix([[1, z], [0, 1]])

def shunt(y):                # a component in parallel
    """
    <---+---< ZA
        Y
    """
    return np.matrix([[1, 0], [1/y, 1]])

def halfpi(y, z):            # a shunt input l-match
    """
    <---+---Z---< ZA
        Y
    """
    return np.matrix([
        [ 1, z ],
        [ 1/y , 1+z/y ]
    ])

def halftee(z, y):           # a series input l-match
    """
    <---Z---+---<  ZA
            Y
    """
    return np.matrix([
        [ 1+z/y, z ],
        [ 1/y , 1 ]
    ])

def tline(deg, zo=50, loss=0): # a transmission line of length deg, db loss
    """
    <----O=======O---< ZA
    """
    theta = -loss / 8.688 + 1j * np.deg2rad(deg)
    return np.matrix([
        [ np.cosh(theta), zo * np.sinh(theta)],
        [ np.sinh(theta) / zo, np.cosh(theta) ]
    ])

def fulltee(z1, z2, z3):     # a tee section
    """
    <---Z1---+---Z3---< ZA
             z2     
    """
    return np.matrix([
        [ 1 + z1/z2, z1 + z3 + z1*z3/z2 ],
        [ 1/z2, 1 + z3/z2 ],
    ])

def fullpi(z1, z2, z3):      # a pi section
    """
    <---+---Z2---+---< ZA
        Z1       Z3
    """
    return np.matrix([
        [ 1 + z2/z3, z2 ],
        [ (z1 + z2 + z3) / (z1*z3), 1 + z2/z1 ]
    ])


# solvers whose results are passed into above ABCD functions
#############################################################

def to_halfwave(zs, za):     # match with a 90 degree tee/pi section
    """
    note, a halfwave tee or pi keeps the phase the same!
    """
    r1, r2 = zs.real, za.real
    x = np.sqrt(r1 * r2) * 1j
    return [[x, -x, x], [-x, x, -x]]
    
def to_halfpi(rin, za):      # match with a shunt input l net, rin > za.real
    """
    """
    ra, xa = za.real, za.imag
    xd = np.sqrt(ra * (rin - ra))
    if np.iscomplex(xd): raise ValueError
    x2 = np.array([-xa - xd, -xa + xd])
    x1 = -(ra**2 + (x2 + xa)**2) / (x2 + xa)
    return np.transpose([x1 * 1j, x2 * 1j]).tolist()

def to_halftee(rin, za):     # match with a series input l net, rin < za.real
    ra, xa = za.real, za.imag
    xd = np.sqrt(rin * ra * (ra**2 + xa**2 - rin * ra))
    if np.iscomplex(xd): raise ValueError
    x2 = np.array([(-rin * xa + xd) / (rin - ra),
                   (-rin * xa - xd) / (rin - ra)])
    x1 = -x2 * (ra**2 + xa * (x2 + xa)) / (ra**2 + (x2 + xa)**2)
    return np.transpose([x1 * 1j, x2 * 1j]).tolist()

def to_fullpi(deg, zo):      # shift phase with a pi section
    zo = zo.real
    theta = np.deg2rad(deg)
    x2 = zo * np.sin(theta)
    x1 = -zo * np.sin(theta) / (1 - np.cos(theta))
    return [x1 * 1j, x2 * 1j, x1 * 1j]

def to_fulltee(deg, zo):     # shift phase with a tee section
    zo = zo.real
    theta = np.deg2rad(deg)
    x2 = -zo / np.sin(theta)
    x1 = zo * (1 - np.cos(theta)) / np.sin(theta)
    return [x1 * 1j, x2 * 1j, x1 * 1j]

def to_shunt(za):            # cancel reactance with a shunt section
    ra, xa = za.real, za.imag
    x1 = -(xa + ra**2 / xa)
    return [x1 * 1j]

def to_series(za):           # cancel reactance with a series section
    x1 = -za.imag
    return [x1 * 1j]


# experimental
########################################

def to_stub1(za, zo=50, shorted=True): # match with a stub-series input 
    """
    -----------------/-----------|
    main line zo    /            za
    ---------------/---/----l----|
                  /   d
                 /___/
    """
    GL = z2g(za, zo)
    thL = np.angle(GL)
    bl = thL / 2 + np.array([ 1, -1 ]) * np.arccos(-abs(GL)) / 2
    if shorted:
        bd = np.arctan(-np.tan(2 * bl - thL) / 2)
    else:
        bd = np.arctan(1 / (np.tan(2 * bl - thL) / 2))
    return np.transpose(np.rad2deg(unwrap([ bd, bl ]))).tolist()


def to_qwt1(za, zo=50):
    """
    ---------------==========---------|
    main line zo       z1       zo    za
    ---------------==========---------|
                   l1=1/4     lm
    """
    lm, zm = lmin(za, zo)
    return np.transpose([ np.sqrt(zo * zm), lm ]).tolist()

def to_qwt2(za, zo=50):
    """
    ---------------==========----|--|
    main line zo       z1        |  za
    ---------------==========-|--|--|
                     L1=1/4   |  |
                              |z2| L2=1/8 shorted or
                              |__|    3/8 opened
    """
    ya = 1 / za
    gl, bl = ya.real, ya.imag
    z1 = np.sqrt(zo / gl)
    z2 = 1 / bl
    return z1, z2

def to_qwt3(za, z2, zo=50):
    """
    ---------------==========----|--|
    main line zo       z1        |  za
    ---------------==========-|--|--|
                     L1=1/4   |  |
                              |z2| d                
                              |__| shorted or opened
    """
    ya = 1 / za
    gl, bl = ya.real, ya.imag
    z1 = sqrt(zo / gl);
    ds = arctan(1 / (bl * z2)) / (2 * np.pi)
    do = arctan(-bl * z2) / (2 * np.pi)
    d = np.mod([ds, do], 0.5)
    return z1, d


# ABCD vector functions
#####################################

def vec(e, i):               # returns a ABCD vector for E, I
    return np.matrix([complex(e), complex(i)]).T

def emag(v):                 # returns the magnitude of E
    return float(np.absolute(v[0]))

def ephase(v):               # returns the phase of E
    return float(np.angle(v[0], deg=True))

def power(*vs):              # power of one or more lines together
    return sum(float(np.absolute(v[1])**2 * impedance(v).real) for v in vs)

def impedance(*vs):          # impedance of one or more lines together
    return 1 / sum(complex(v[1] / v[0]) for v in vs)

def emax(v, zo=50):          # maximum rms voltage on transmission line
    gm = z2g(impedance(v), zo)
    return np.sqrt(power(v) * zo * np.abs(swr(gm)))


# helper functions
######################################

def qmin(zs, za):            # q for L network, or minimum q for tee or pi
    rp = max(zs.real, za.real)
    rs = min(zs.real, za.real)
    return np.sqrt(rp / rs - 1)

def qmin2(zs, za):           # minimum q for a LL network
    rp = max(zs.real, za.real)
    rv = np.sqrt(rp * rs)
    return np.sqrt(rp / rv - 1)

def opened_stub(deg, zo=50):
    theta = np.deg2rad(deg)
    return -1j * zo / np.tan(theta)

def shorted_stub(deg, zo=50):
    theta = np.deg2rad(deg)
    return 1j * zo * np.tan(theta)

def component_value(impedance, fd):
    w = 2 * np.pi * fd
    x = impedance.imag
    return 1 / (w * x) if x < 0 else x / w

def qlosses(*impedances, q=200):
    return [ z + (z.imag / q if z.imag > 0 else 0) for z in impedances ]

def parallel(*impedances):
    return 1 / sum(1 / x for x in impedances)

def z2g(z, zo=50):
    return (z - zo) / (z + zo)

def g2z(gm, zo=50):
    return zo * (1 + gm) / (1 - gm)

def swr(gm):
    return (1 + abs(gm)) / (1 - abs(gm))

def s2p(z):                  # serial to parallel
    zp = 1/z
    return 1/zp.real - 1j/zp.imag

def unwrap(theta):           # convert theta rads to between 0 and 2*pi
    return np.mod(theta, 2 * np.pi)

def lmin(za, zo=50):         # distance to voltage min/max
    gm = z2g(za, zo)
    th = np.angle(gm)
    lm = np.array([ (th + np.pi) / 2, unwrap(th) / 2 ])
    zm = np.array([ zo / swr(gm), zo * swr(gm) ])
    return np.rad2deg(lm), zm


# print functions
####################################

def polar(x):
    return "%.4f / %.4f" % (np.absolute(x), np.angle(x, deg=True))

def status(v, note=0):
    print("p({}) = {:.4f}".format(note, power(v)))
    print("z({}) = {:.4f}".format(note, complex(v[0]/v[1])))
    print("i({}) = {}".format(note, polar(v[1])))
    print("e({}) = {}".format(note, polar(v[0])))
    print()

def component(values, fd, precision=4):
    return [ notation(component_value(z, fd), precision, units='FH')
             for z in values ]

def notation(x, precision=4, units=None):
    SUFFIX = ["p", "n", "u", "m", "", "k", "M", "G"]
    exp = np.floor(np.log10(np.absolute(x)))
    mant = round(x / 10**exp, precision-1)
    p = int(exp // 3)
    value = (mant * 10**exp) / 10**(3 * p)
    unit = units[0 if x < 0 else 1] if units else ''
    return "%g%s%s" % (np.absolute(value), SUFFIX[p-4], unit)




# undocumented   
########################################

# fix
def to_halfpi2(zs, za, solution=(0,0)):  # match with a double shunt l-match
    rin = np.sqrt(zs.real * za.real)
    x = to_halftee(rin, zs, solution=solution[0])
    y = to_halfpi(rin, za, solution=solution[1])
    return x[1], x[0], y[0], y[1]

# fix
def to_halftee2(zs, za, solution=(0,0)): # match with a double series l-match
    rin = np.sqrt(zs.real * za.real)
    x = to_halfpi(rin, zs, solution=solution[0])
    y = to_halftee(rin, za, solution=solution[1])
    return x[1], x[0], y[0], y[1]
    
# fix
def to_fullpi2(zs, za, q=0, solution=(0,0)):   # match with a pi section
    q = q or qmin(zs, za) + 1e-9
    rin = max(zs.real, za.real) / (q**2 + 1)
    x = to_halftee(rin, zs, solution=solution[0])
    y = to_halftee(rin, za, solution=solution[1])
    return x[1], x[0] + y[0], y[1]

# fix
def to_fulltee2(zs, za, q=0, solution=(0,0)):  # match with a tee section
    q = q or qmin(zs, za) + 1e-9
    rin = min(zs.real, za.real) * (q**2 + 1)
    x = to_halfpi(rin, zs, solution=solution[0])
    y = to_halfpi(rin, za, solution=solution[1])
    return x[1], parallel(x[0], y[0]), y[1]
  
def to_resistive_halftee(rin, ra):    # rin > ra
    r2 = ra - np.sqrt(rin / (rin - ra))
    r1 = rin - (ra * r2) / (ra + r2)
    return [r1, r2]

def to_resistive_halfpi(rin, ra):     # rin < ra
    return to_resistive_halftee(ra, rin)[::-1]

def halftee2(z1, y1, z2, y2): # a double series input l-match
    return halftee(z1, y1) * halftee(z2, y2)

def halfpi2(z1, y1, z2, y2):  # a double shunt input l-match
    return halfpi(z1, y1) * halfpi(z2, y2)

