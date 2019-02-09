
import numpy as np
from numpy import inf

# generates ABCD matrices
######################################

def auto(ratio, xt, k=1):       # a 1:n autotransformer
    n = ratio / (1 - ratio)
    x1 = xt / (1 + n**2 + 2 * k * n)
    x2 = x1 * n**2
    xm = k * n * x1
    return fulltee((x1 + xm) * 1j, (x2 + xm) * 1j, -xm * 1j)

def mutual(n, x1, k=1, q=None): # a 1:n transformer
    x2 = x1 * n**2
    xm = k * n * x1
    r1 = x1 / q if q else 0
    r2 = x2 / q if q else 0
    return fulltee(r1 + (x1 - xm) * 1j, xm * 1j, r2 + (x2 - xm) * 1j)

def trans(n):                   # an ideal 1:n transformer
    return np.matrix([[ 1/n, 0], [0, n]])
    
def series(z):                  # a component in series
    """
    <---Z---< ZA
    """
    return np.matrix([[1, z], [0, 1]])

def shunt(y):                   # a component in parallel
    """
    <---+---< ZA
        Y
    """
    return np.matrix([[1, 0], [1/y, 1]])

def halfpi(y, z):               # a shunt input l-match
    """
    <---+---Z---< ZA
        Y
    """
    return np.matrix([
        [ 1, z ],
        [ 1/y , 1+z/y ]
    ])

def halftee(z, y):              # a series input l-match
    """
    <---Z---+---<  ZA
            Y
    """
    return np.matrix([
        [ 1+z/y, z ],
        [ 1/y , 1 ]
    ])

def tline(deg, zo=50, loss=0):  # a transmission line of length deg, db loss
    """
    <----O=======O---< ZA
    """
    theta = -loss / 8.688 + 1j * np.deg2rad(deg)
    return np.matrix([
        [ np.cosh(theta), zo * np.sinh(theta)],
        [ np.sinh(theta) / zo, np.cosh(theta) ]
    ])

def fulltee(z1, z2, z3):        # a tee section
    """
    <---Z1---+---Z3---< ZA
             z2     
    """
    return np.matrix([
        [ 1 + z1/z2, z1 + z3 + z1*z3/z2 ],
        [ 1/z2, 1 + z3/z2 ],
    ])

def fullpi(z1, z2, z3):         # a pi section
    """
    <---+---Z2---+---< ZA
        Z1       Z3
    """
    return np.matrix([
        [ 1 + z2/z3, z2 ],
        [ (z1 + z2 + z3) / (z1*z3), 1 + z2/z1 ]
    ])


# solvers whose results are passed into the above ABCD functions
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


# experimental solvers
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
    bl = thL / 2 + np.array([1, -1]) * np.arccos(-abs(GL)) / 2
    if shorted:
        bd = np.arctan(-np.tan(2 * bl - thL) / 2)
    else:
        bd = np.arctan(1 / (np.tan(2 * bl - thL) / 2))
    d = np.mod([ bd, bl ], np.pi)
    d = np.rad2deg(d)
    return np.transpose(d).tolist()

def to_lmin(za, zo=50):         # distance to voltage min/max
    gm = z2g(za, zo)
    th = np.angle(gm)
    zm = np.array([ zo / swr(gm), zo * swr(gm) ])
    lm = np.array([ (th + np.pi) / 2, np.mod(th, np.pi) / 2 ])
    lm = np.rad2deg(lm)
    return lm, zm

def to_qwt1(za, zo=50):
    """
    ---------------==========---------|
    main line zo       z1       zo    za
    ---------------==========---------|
                     l1=1/4     lm
    """
    lm, zm = to_lmin(za, zo)
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
    return [[ z1, z2, 45 ],  # shorted
            [ z1, z2, 135 ]] # opened

def to_qwt3(za, z2, zo=50):
    """
    ---------------==========----|--|
    main line zo       z1        |  za
    ---------------==========-|--|--|
                     L1=1/4   |  |
                              |z2| d=[shorted,opened]
                              |__|
    """
    ya = 1 / za
    gl, bl = ya.real, ya.imag
    z1 = np.sqrt(zo / gl)
    d = np.arctan([ 1 / (bl * z2), -bl * z2 ])
    d = np.mod(d, np.pi)
    d = np.rad2deg(d)
    return [[ z1, d[0] ], # shorted
            [ z1, d[1] ]] # opened


# ABCD vector functions
#####################################

def vec(e, i):               # returns a ABCD vector for E, I
    return np.matrix([complex(e), complex(i)]).T

def emag(v):                 # returns the magnitude of E
    return float(np.absolute(v[0]))

def ephase(v):               # returns the phase of E
    return float(np.angle(v[0], deg=True))

def power(*vs):              # power of one or more lines together in parallel
    return sum(float(np.absolute(v[1])**2 * impedance(v).real) for v in vs)

def impedance(*vs):          # impedance of one or more lines together in parallel
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

def qlosses(impedances, q=200):
    return [ z + (z.imag / q if z.imag > 0 else 0) for z in impedances ]

def parallel(*impedances):
    return 1 / sum(1 / x for x in impedances)

def z2g(z, zo=50):
    return (z - zo) / (z + zo)

def g2z(gm, zo=50):
    return zo * (1 + gm) / (1 - gm)

def swr(gm):
    return (1 + abs(gm)) / (1 - abs(gm))

def s2p(z):                  # series to parallel
    zp = 1 / z
    return 1 / zp.real - 1j / zp.imag


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
    return v

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


# fix
########################################

def to_halfpi2(zs, za, solution=(0,0)):  # match with a shunt LL-match
    rin = np.sqrt(zs.real * za.real)
    x = to_halftee(rin, zs, solution=solution[0])
    y = to_halfpi(rin, za, solution=solution[1])
    return x[1], x[0], y[0], y[1]

def to_halftee2(zs, za, solution=(0,0)): # match with a series LL-match
    rin = np.sqrt(zs.real * za.real)
    x = to_halfpi(rin, zs, solution=solution[0])
    y = to_halftee(rin, za, solution=solution[1])
    return x[1], x[0], y[0], y[1]
    
def to_fullpi2(zs, za, q=0, solution=(0,0)):   # match with a pi section
    q = q or qmin(zs, za) + 1e-9
    rin = max(zs.real, za.real) / (q**2 + 1)
    x = to_halftee(rin, zs, solution=solution[0])
    y = to_halftee(rin, za, solution=solution[1])
    return x[1], x[0] + y[0], y[1]

def to_fulltee2(zs, za, q=0, solution=(0,0)):  # match with a tee section
    q = q or qmin(zs, za) + 1e-9
    rin = min(zs.real, za.real) * (q**2 + 1)
    x = to_halfpi(rin, zs, solution=solution[0])
    y = to_halfpi(rin, za, solution=solution[1])
    return x[1], parallel(x[0], y[0]), y[1]
  
def qmin2(zs, za):           # minimum q for a LL network
    rp = max(zs.real, za.real)
    rv = np.sqrt(rp * rs)
    return np.sqrt(rp / rv - 1)

# undocumented
########################################

def to_resistive_halftee(rin, ra):    # rin > ra
    r2 = ra - np.sqrt(rin / (rin - ra))
    r1 = rin - (ra * r2) / (ra + r2)
    return [r1, r2]

def to_resistive_halfpi(rin, ra):     # rin < ra
    return to_resistive_halftee(ra, rin)[::-1]

def halftee2(z1, y1, z2, y2): # a series input LL-match
    return halftee(z1, y1) * halftee(z2, y2)

def halfpi2(z1, y1, z2, y2):  # a shunt input LL-match
    return halfpi(z1, y1) * halfpi(z2, y2)

# amp models

def transconductance(mode, IC=None, ID=None, VP=None, IDSS=None):
    if mode == 'bjt':
        return IC / 26
    elif mode == 'fet':
        return -2 * np.sqrt(ID * IDSS) / VP
    raise ValueError

def hybrid(ai=None, av=None, rin=None, rout=None):
    return np.matrix([
        [ rin,  ai ],
        [ 1/av, 1/rout ]])

def common_source(RD, RS=0, ID=1, VP=-6, IDSS=8):
    gm = transconductance('fet', ID=ID, VP=VP, IDSS=IDSS)
    rs = RS + 1 / gm
    return hybrid(ai=inf, av=-RD/rs, rin=inf, rout=RD)

def common_drain(RS, ID=1, VP=-6, IDSS=8):
    gm = transconductance('fet', ID=ID, VP=VP, IDSS=IDSS)
    rs = RS + 1 / gm
    return hybrid(ai=-inf, av=RS/rs, rin=inf, rout=1/gm)

def common_gate(RD, ID=1, VP=-6, IDSS=8):
    gm = transconductance('fet', ID=ID, VP=VP, IDSS=IDSS)
    return hybrid(ai=-1, av=gm*RD, rin=1/gm, rout=RD)

def common_emitter(RC, RE=0, IC=1, beta=100, ft=300, f=0):
    gm = transconductance('bjt', IC=IC)
    beta /= (1 + 1j * beta * f / ft)
    re = RE + 1 / gm
    return hybrid(ai=beta, av=-RC/re, rin=beta*re, rout=RC)

def common_collector(RE, IC=1, beta=100, ft=300, f=0):
    gm = transconductance('bjt', IC=IC)
    beta /= (1 + 1j * beta * f / ft)
    re = RE + 1 / gm
    return hybrid(ai=-beta, av=RE/re, rin=beta*re, rout=1/gm)

def common_base(RC, IC=1):
    gm = transconductance('bjt', IC=IC)
    return hybrid(ai=-1, av=gm*RC, rin=1/gm, rout=RC)

def feedback_amplifier(RE=0, RF=0, RL=0, RS=0, mode='bjt',
                       IC=1, ID=1, VP=-6, IDSS=8):
    gm = transconductance(mode, IC=IC, ID=ID, IDSS=8, VP=VP)
    RD = RE + 1 / gm
    Gv = -RL * (RF - RD) / RD / (RL + RF)
    Zin = RD * (RL + RF) / (RL + RD)
    Zout = RD * (RF + RS) / (RD + RS)
    return Gv, Zin, Zout

# biasing

def fet_self_bias(ID=1, VP=-6, IDSS=8):
    RS = VP / ID * (np.sqrt(ID / IDSS) - 1)
    return RS

def npn_feedback_bias(RC, RE=0, RB=inf, IC=1, VCC=12, beta=100):
    IC = IC / 1000
    ib = IC / beta
    vb = (IC + ib) * RE + .7
    ibb = ib + vb / RB
    vc = VCC - RC * (IC + ibb) 
    RF = (vc - vb) / ibb
    return RF

def fet_divider_bias(VG, ID=1, VP=-6, IDSS=8):
    RS = VG / ID
    return RS



