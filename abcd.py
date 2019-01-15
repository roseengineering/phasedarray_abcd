import numpy as np

# generates ABCD matrix
######################################

def auto(ratio, xt, k=1):
    n = ratio / (1 - ratio)
    x1 = xt / (1 + n**2 + 2 * k * n)
    x2 = x1 * n**2
    xm = k * n * x1
    return fulltee((x1 + xm) * 1j, (x2 + xm) * 1j, -xm * 1j)

def mutual(n, x1, k=1):
    x2 = x1 * n**2
    xm = k * n * x1
    return fulltee((x1 - xm) * 1j, xm * 1j, (x2 - xm) * 1j)

def trans(n):
    return np.matrix([[ 1/n, 0], [0, n]])
    
# <---Z---< ZL
def series(z):
    return np.matrix([[1, z], [0, 1]])

# <---+---< ZL
#     Y
def shunt(y):
    return np.matrix([[1, 0], [1/y, 1]])

# <---+---Z---< ZL
#     Y
def halfpi(y, z):
    return np.matrix([
        [ 1, z ],
        [ 1/y , 1+z/y ]
    ])

# <---Z---+---<  ZL
#         Y
def halftee(z, y):
    return np.matrix([
        [ 1+z/y, z ],
        [ 1/y , 1 ]
    ])

# <----O=======O---< ZL
def tline(deg, zo=50):
    theta = np.deg2rad(deg)
    return np.matrix([
        [ np.cos(theta), 1j*zo*np.sin(theta)],
        [ 1j*np.sin(theta)/zo, np.cos(theta) ]
    ])

# <---Z1---+---Z3---< ZL
#          z2     
def fulltee(z1, z2, z3):
    return np.matrix([
        [ 1 + z1/z2, z1 + z3 + z1*z3/z2 ],
        [ 1/z2, 1 + z3/z2 ],
    ])

# <---+---Z2---+---< ZL
#     Z1       Z3
def fullpi(z1, z2, z3):
    return np.matrix([
        [ 1 + z2/z3, z2 ],
        [ (z1 + z2 + z3) / (z1*z3), 1 + z2/z1 ]
    ])

def gain(db):
    g = 10 ** (db / 10)
    return np.matrix([
        [ np.sqrt(g), 0 ],
        [ 0, np.sqrt(g) ]
    ])

def halftee2(z1, y1, z2, y2):
    return halftee(z1, y1) * halftee(z2, y2)

def halfpi2(z1, y1, z2, y2):
    return halfpi(z1, y1) * halfpi(z2, y2)

# circuit solutions
######################################

# it keeps the phase the same!
def to_halfwave(zs, za, solution=0):
    r1, r2 = zs.real, za.real
    x = np.sqrt(r1 * r2) * 1j
    if solution: x = -x
    return x, -x, x
    
def to_halfpi(rin, za, solution=0):
    ra, xa = za.real, za.imag
    xd = np.sqrt(ra * (rin - ra))
    if np.iscomplex(xd): raise ValueError
    x2 = np.array([-xa - xd, -xa + xd])
    x2 = x2[solution]
    x1 = -(ra**2 + (x2 + xa)**2) / (x2 + xa)
    return x1 * 1j, x2 * 1j

def to_halftee(rin, za, solution=0):
    ra, xa = za.real, za.imag
    xd = np.sqrt(rin * ra * (ra**2 + xa**2 - rin * ra))
    if np.iscomplex(xd): raise ValueError
    x2 = np.array([(-rin * xa + xd) / (rin - ra),
                   (-rin * xa - xd) / (rin - ra)])
    x2 = x2[solution]
    x1 = -x2 * (ra**2 + xa * (x2 + xa)) / (ra**2 + (x2 + xa)**2)
    return x1 * 1j, x2 * 1j

def to_fullpi(deg, zo):
    zo = zo.real
    theta = np.deg2rad(deg)
    x2 = zo * np.sin(theta)
    x1 = -zo * np.sin(theta) / (1 - np.cos(theta))
    return x1 * 1j, x2 * 1j, x1 * 1j

def to_fulltee(deg, zo):
    zo = zo.real
    theta = np.deg2rad(deg)
    x2 = -zo / np.sin(theta)
    x1 = zo * (1 - np.cos(theta)) / np.sin(theta)
    return x1 * 1j, x2 * 1j, x1 * 1j

def to_halfpi2(zs, za, solution=(0,0)):
    rin = np.sqrt(zs.real * za.real)
    x = to_halftee(rin, zs, solution=solution[0])
    y = to_halfpi(rin, za, solution=solution[1])
    return x[1], x[0], y[0], y[1]

def to_halftee2(zs, za, solution=(0,0)):
    rin = np.sqrt(zs.real * za.real)
    x = to_halfpi(rin, zs, solution=solution[0])
    y = to_halftee(rin, za, solution=solution[1])
    return x[1], x[0], y[0], y[1]
    
def to_fullpi2(zs, za, q=0, solution=(0,0)):
    q = q or qmin(zs, za) + 1e-9
    rin = max(zs.real, za.real) / (q**2 + 1)
    x = to_halftee(rin, zs, solution=solution[0])
    y = to_halftee(rin, za, solution=solution[1])
    return x[1], x[0] + y[0], y[1]

def to_fulltee2(zs, za, q=0, solution=(0,0)):
    q = q or qmin(zs, za) + 1e-9
    rin = min(zs.real, za.real) * (q**2 + 1)
    x = to_halfpi(rin, zs, solution=solution[0])
    y = to_halfpi(rin, za, solution=solution[1])
    return x[1], parallel(x[0], y[0]), y[1]
     
def to_shunt(za):
    ra, xa = za.real, za.imag
    x1 = -(xa + ra**2 / xa)
    return (x1 * 1j,)

def to_series(za):
    x1 = -za.imag
    return (x1 * 1j,)

######################################

# q for L network, or minimum q for tee or pi
def qmin(zs, za):
    rp = max(zs.real, za.real)
    rs = min(zs.real, za.real)
    return np.sqrt(rp / rs - 1)

# minimum q for a LL network
def qmin2(zs, za):
    rp = max(zs.real, za.real)
    rv = np.sqrt(rp * rs)
    return np.sqrt(rp / rv - 1)

def openstub(deg, zo):
    theta = np.deg2rad(deg)
    return -1j * zo * np.cot(theta)

def shortstub(deg, zo):
    theta = np.deg2rad(deg)
    return 1j * zo * np.tan(theta)

def reactance_value(component, fd):
    w = 2 * np.pi * fd
    return 1j / (w * component) if component < 0 else 1j * w * component

def component_value(impedance, fd):
    w = 2 * np.pi * fd
    return 1 / (w * impedance.imag) if impedance.imag < 0 else impedance.imag / w  

def parallel(*impedances):
    return 1 / sum(1 / x for x in impedances)

# for line

def vec(e, i):
    return np.matrix([complex(e), complex(i)]).T

def emag(v):
    return float(np.absolute(v[0]))

def ephase(v):
    return float(np.angle(v[0], deg=True))

def power(*vectors):
    return sum(float(np.absolute(v[1])**2 * impedance(v).real) for v in vectors)

def impedance(*vectors):
    return 1 / sum(complex(v[1] / v[0]) for v in vectors)


# for printing

def polar(x):
    return "%.4f / %.4f" % (np.absolute(x), np.angle(x, deg=True))

def status(v, note=0):
    print("p({}) = {:.4f}".format(note, power(v)))
    print("z({}) = {:.4f}".format(note, complex(v[0]/v[1])))
    print("i({}) = {}".format(note, polar(v[1])))
    print("e({}) = {}".format(note, polar(v[0])))
    print()

def reactance(values, fd, precision=4):
    return [ notation(reactance_value(z, fd), precision) for z in values ]

def component(values, fd, precision=4):
    return [ notation(component_value(z, fd), precision) for z in values ]

def notation(x, precision=4):
    SUFFIX = ["p", "n", "u", "m", "", "k", "M", "G"]
    exp = np.floor(np.log10(np.absolute(x)))
    mant = round(x / 10**exp, precision-1)
    p = int(exp // 3)
    value = (mant * 10**exp) / 10**(3 * p)
    return "%g%s%s" % (np.absolute(value), SUFFIX[p-4], 'F' if x < 0 else 'H')

