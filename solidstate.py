
import numpy as np
from numpy import inf

def hybrid(ai=None, av=None, rin=None, rout=None):
    return np.matrix([
        [ rin,  ai ],
        [ 1/av, 1/rout ]])

# amp models

def transconductance(mode='bjt', **kw):
    if mode == 'bjt':
        IC = kw.get('IC', 1)     # ma
        gm = IC / 26
        return gm
    elif mode == 'fet':
        ID = kw.get('ID', 1)     # ma
        IDSS = kw.get('IDSS', 8) # ma
        VP = kw.get('VP', -6)    # v
        gm = -2 * np.sqrt(ID * IDSS) / VP
        return gm / 1000
    elif mode == 'mos':
        ID = kw.get('ID', 1)     # ma
        K = kw.get('K', 1.5)     # ma/v^2
        gm = 2 * np.sqrt(K * ID)
        return gm / 1000
    raise ValueError

def common_source(RD, RS=0, **kw):
    gm = transconductance('fet', **kw)
    rs = RS + 1 / gm
    return hybrid(ai=inf, av=-RD/rs, rin=inf, rout=RD)

def common_drain(RS, **kw):
    gm = transconductance('fet', **kw)
    rs = RS + 1 / gm
    return hybrid(ai=-inf, av=RS/rs, rin=inf, rout=1/gm)

def common_gate(RD, **kw):
    gm = transconductance('fet', **kw)
    return hybrid(ai=-1, av=gm*RD, rin=1/gm, rout=RD)

def common_emitter(RC, RE=0, beta=100, f=0, ft=300, **kw):
    gm = transconductance('bjt', **kw)
    beta /= (1 + 1j * beta * f / ft)
    re = RE + 1 / gm
    return hybrid(ai=beta, av=-RC/re, rin=beta*re, rout=RC)

def common_collector(RE, beta=100, f=0, ft=300, **kw):
    gm = transconductance('bjt', **kw)
    beta /= (1 + 1j * beta * f / ft)
    re = RE + 1 / gm
    return hybrid(ai=-beta, av=RE/re, rin=beta*re, rout=1/gm)

def common_base(RC, **kw):
    gm = transconductance('bjt', **kw)
    return hybrid(ai=-1, av=gm*RC, rin=1/gm, rout=RC)

def fba(RD=0, RF=0, RL=0, RS=0, N=1):
    RD /= N
    Gv = -RL * (RF - RD) / RD / (RL + RF)
    Zin = RD * (RL + RF) / (RL + RD)
    Zout = RD * (RF + RS) / (RD + RS)
    return Gv, Zin, Zout

# biasing
# XCB << RS so let XCB < RS/10 at fmin (common emitter)
# XCE << RS+re||RE so let XCE < (RS+re)/10 at fmin (common base)
# RF||RB << (BETA+1)*RE so let RB||RF < (BETA+1)*RE/10

def bias_bjt_feedback(RC, RE=0, RB=inf, IC=1, VCC=12, beta=100):
    VBE = .7
    IC = IC / 1000
    ib = IC / beta
    vb = (IC + ib) * RE + VBE
    ibb = ib + vb / RB
    vc = VCC - RC * (IC + ibb) 
    RF = (vc - vb) / ibb
    N = (beta + 1) * RE / parallel(RF, RB)
    return RF, N

def bias_fet_divider(N=10, ID=1, VP=-6, IDSS=8):
    # ID = IDSS * (1 - VGS / VP)^2, solving for VGS
    VGS = VP * (1 - np.sqrt(ID / IDSS))
    VG = -VGS * N + VGS
    RS = -VGS * N / ID
    return VG, RS

def bias_mos_divider(N=10, ID=1, VTH=1, K=1.5):
    # ID = K / 2 * (VGS - VTH)**2, solving for VGS
    VGS = np.sqrt(2 * ID / K) + VTH
    VG = VGS * N + VGS
    RS = VGS * N / ID
    return VG, RS


