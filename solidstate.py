
import numpy as np
from numpy import inf

DEFAULT_IC = 1     # ma
DEFAULT_IDSS = 10  # ma (mpf102)
DEFAULT_VP = -4    # v (mpf102)
DEFAULT_K = 100    # ma/v^2 (2n7000)
DEFAULT_VTH = 1
DEFAULT_BETA = 100
DEFAULT_FT = 300

def hybrid(ai=0, av=0, rin=0, rout=0):
    return np.matrix([
        [ rin,  ai ],
        [ 1/av, 1/rout ]])

# amp models

def transconductance(mode='bjt', **kw):
    # returns gm in siemens
    # note, ic and id are in ma
    IC = kw.get('IC') or DEFAULT_IC
    ID = kw.get('ID') or IC
    if mode == 'bjt':
        gm = IC / 26
        return gm
    elif mode == 'fet':
        IDSS = kw.get('IDSS') or DEFAULT_IDSS
        VP = kw.get('VP') or DEFAULT_VP
        gm = -2 * np.sqrt(ID * IDSS) / VP
        return gm / 1000
    elif mode == 'mos':
        # ID = K * (VG - VTH)^2
        # so ID^.5 = K^.5 * (VG - VTH)
        # therefore K^.5 = d(ID^.5) / d(VG)
        K = kw.get('K') or DEFAULT_K
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

def common_emitter(RC, RE=0, beta=DEFAULT_BETA, f=0, ft=DEFAULT_FT, **kw):
    gm = transconductance('bjt', **kw)
    beta /= (1 + 1j * beta * f / ft)
    re = RE + 1 / gm
    return hybrid(ai=beta, av=-RC/re, rin=beta*re, rout=RC)

def common_collector(RE, beta=DEFAULT_BETA, f=0, ft=DEFAULT_FT, **kw):
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

# bjt feedback bias
# XCB << RS so let XCB < RS/10 at fmin (common emitter)
# XCE << RS+re||RE so let XCE < (RS+re)/10 at fmin (common base)
# RF||RB << (BETA+1)*RE so let RB||RF < (BETA+1)*RE/10

def bias_bjt_feedback(RC, RE=0, RB=inf, IC=DEFAULT_IC, VCC=12, beta=DEFAULT_BETA):
    VBE = .7
    IC = IC / 1000
    ib = IC / beta
    vb = (IC + ib) * RE + VBE
    ibb = ib + vb / RB
    vc = VCC - RC * (IC + ibb) 
    RF = (vc - vb) / ibb
    N = (beta + 1) * RE / parallel(RF, RB)
    return RF, N

# fet divider bias
# RS * ID >> VGS for best feedback
# if RS * ID >> VGS and VG = RS * ID + VGS 
#       then VG ~= RS * ID
#       so ID ~= VG / RS 
#       so RS ~= VG / ID

def bias_fet_divider(N=10, ID=DEFAULT_IC, VP=DEFAULT_VP, IDSS=DEFAULT_IDSS):
    # ID = IDSS * (1 - VGS / VP)^2, solving for VGS
    VGS = VP * (1 - np.sqrt(ID / IDSS))
    VG = -VGS * N + VGS
    RS = -VGS * N / ID
    return VG, RS

# mosfet divider bias
# RS * ID >> VGS for best feedback
# if RS * ID >> VGS and VG = RS * ID + VGS 
#       then VG ~= RS * ID
#       so ID ~= VG / RS 
#       so RS ~= VG / ID

def bias_mos_divider(N=10, ID=DEFAULT_IC, VTH=DEFAULT_VTH, K=DEFAULT_K):
    # ID = K / 2 * (VGS - VTH)**2, solving for VGS
    VGS = np.sqrt(2 * ID / K) + VTH
    VG = VGS * N + VGS
    RS = VGS * N / ID
    return VG, RS

