# This file was automatically created by FeynRules 2.0.23
# Mathematica version: 9.0 for Mac OS X x86 (64-bit) (November 20, 2012)
# Date: Sat 20 Sep 2014 16:11:37



from object_library import all_parameters, Parameter


from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot

# This is a default parameter object representing 0.
ZERO = Parameter(name = 'ZERO',
                 nature = 'internal',
                 type = 'real',
                 value = '0.0',
                 texname = '0')

# User-defined parameters.
cabi = Parameter(name = 'cabi',
                 nature = 'external',
                 type = 'real',
                 value = 0.488,
                 texname = '\\theta _c',
                 lhablock = 'CKMBLOCK',
                 lhacode = [ 1 ])

mZinput = Parameter(name = 'mZinput',
                    nature = 'external',
                    type = 'real',
                    value = 91.188,
                    texname = '\\text{mZinput}',
                    lhablock = 'GAUGEMASS',
                    lhacode = [ 1 ])

mZDinput = Parameter(name = 'mZDinput',
                     nature = 'external',
                     type = 'real',
                     value = 20,
                     texname = '\\text{mZDinput}',
                     lhablock = 'HIDDEN',
                     lhacode = [ 1 ])

MHSinput = Parameter(name = 'MHSinput',
                     nature = 'external',
                     type = 'real',
                     value = 200,
                     texname = '\\text{MHSinput}',
                     lhablock = 'HIDDEN',
                     lhacode = [ 2 ])

epsilon = Parameter(name = 'epsilon',
                    nature = 'external',
                    type = 'real',
                    value = 0.01,
                    texname = '\\epsilon',
                    lhablock = 'HIDDEN',
                    lhacode = [ 3 ])

kap = Parameter(name = 'kap',
                nature = 'external',
                type = 'real',
                value = 1.e-9,
                texname = '\\text{kap}',
                lhablock = 'HIDDEN',
                lhacode = [ 4 ])

aXM1 = Parameter(name = 'aXM1',
                 nature = 'external',
                 type = 'real',
                 value = 127.9,
                 texname = '\\text{aXM1}',
                 lhablock = 'HIDDEN',
                 lhacode = [ 5 ])

MHinput = Parameter(name = 'MHinput',
                    nature = 'external',
                    type = 'real',
                    value = 125,
                    texname = '\\text{MHinput}',
                    lhablock = 'HIGGS',
                    lhacode = [ 1 ])

swsq = Parameter(name = 'swsq',
                 nature = 'external',
                 type = 'real',
                 value = 0.225,
                 texname = '\\text{swsq}',
                 lhablock = 'SMINPUTS',
                 lhacode = [ 1 ])

aEWM1 = Parameter(name = 'aEWM1',
                  nature = 'external',
                  type = 'real',
                  value = 127.9,
                  texname = '\\text{aEWM1}',
                  lhablock = 'SMINPUTS',
                  lhacode = [ 2 ])

Gf = Parameter(name = 'Gf',
               nature = 'external',
               type = 'real',
               value = 0.000011663900000000002,
               texname = '\\text{Gf}',
               lhablock = 'SMINPUTS',
               lhacode = [ 3 ])

aS = Parameter(name = 'aS',
               nature = 'external',
               type = 'real',
               value = 0.118,
               texname = '\\text{aS}',
               lhablock = 'SMINPUTS',
               lhacode = [ 4 ])

ymc = Parameter(name = 'ymc',
                nature = 'external',
                type = 'real',
                value = 1.42,
                texname = '\\text{ymc}',
                lhablock = 'YUKAWA',
                lhacode = [ 4 ])

ymb = Parameter(name = 'ymb',
                nature = 'external',
                type = 'real',
                value = 4.7,
                texname = '\\text{ymb}',
                lhablock = 'YUKAWA',
                lhacode = [ 5 ])

ymt = Parameter(name = 'ymt',
                nature = 'external',
                type = 'real',
                value = 174.3,
                texname = '\\text{ymt}',
                lhablock = 'YUKAWA',
                lhacode = [ 6 ])

ymel = Parameter(name = 'ymel',
                 nature = 'external',
                 type = 'real',
                 value = 0.000511,
                 texname = '\\text{ymel}',
                 lhablock = 'YUKAWA',
                 lhacode = [ 11 ])

ymmu = Parameter(name = 'ymmu',
                 nature = 'external',
                 type = 'real',
                 value = 0.1057,
                 texname = '\\text{ymmu}',
                 lhablock = 'YUKAWA',
                 lhacode = [ 13 ])

ymtau = Parameter(name = 'ymtau',
                  nature = 'external',
                  type = 'real',
                  value = 1.777,
                  texname = '\\text{ymtau}',
                  lhablock = 'YUKAWA',
                  lhacode = [ 15 ])

ME = Parameter(name = 'ME',
                nature = 'external',
                type = 'real',
                value = 0.000511,
                texname = '\\text{ME}',
                lhablock = 'MASS',
                lhacode = [ 11 ])
                
MM = Parameter(name = 'MM',
                nature = 'external',
                type = 'real',
                value = 0.1057,
                texname = '\\text{MM}',
                lhablock = 'MASS',
                lhacode = [ 13 ])
                
MTA = Parameter(name = 'MTA',
                nature = 'external',
                type = 'real',
                value = 1.777,
                texname = '\\text{MTA}',
                lhablock = 'MASS',
                lhacode = [ 15 ])

MC = Parameter(name = 'MC',
               nature = 'external',
               type = 'real',
               value = 1.42,
               texname = '\\text{MC}',
               lhablock = 'MASS',
               lhacode = [ 4 ])

MT = Parameter(name = 'MT',
               nature = 'external',
               type = 'real',
               value = 174.3,
               texname = '\\text{MT}',
               lhablock = 'MASS',
               lhacode = [ 6 ])

MB = Parameter(name = 'MB',
               nature = 'external',
               type = 'real',
               value = 4.7,
               texname = '\\text{MB}',
               lhablock = 'MASS',
               lhacode = [ 5 ])

WT = Parameter(name = 'WT',
               nature = 'external',
               type = 'real',
               value = 1.50833649,
               texname = '\\text{WT}',
               lhablock = 'DECAY',
               lhacode = [ 6 ])

WZ = Parameter(name = 'WZ',
               nature = 'external',
               type = 'real',
               value = 2.44140351,
               texname = '\\text{WZ}',
               lhablock = 'DECAY',
               lhacode = [ 23 ])

WZp = Parameter(name = 'WZp',
                nature = 'external',
                type = 'real',
                value = 0.0008252,
                texname = '\\text{WZp}',
                lhablock = 'DECAY',
                lhacode = [ 1023 ])

WW = Parameter(name = 'WW',
               nature = 'external',
               type = 'real',
               value = 2.04759951,
               texname = '\\text{WW}',
               lhablock = 'DECAY',
               lhacode = [ 24 ])

WH = Parameter(name = 'WH',
               nature = 'external',
               type = 'real',
               value = 0.00282299,
               texname = '\\text{WH}',
               lhablock = 'DECAY',
               lhacode = [ 25 ])

WHS = Parameter(name = 'WHS',
                nature = 'external',
                type = 'real',
                value = 5.23795,
                texname = '\\text{WHS}',
                lhablock = 'DECAY',
                lhacode = [ 35 ])

cw = Parameter(name = 'cw',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(1 - swsq)',
               texname = 'c_w')

sw = Parameter(name = 'sw',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(swsq)',
               texname = 's_w')

aEW = Parameter(name = 'aEW',
                nature = 'internal',
                type = 'real',
                value = '1/aEWM1',
                texname = '\\text{aEW}')

G = Parameter(name = 'G',
              nature = 'internal',
              type = 'real',
              value = '2*cmath.sqrt(aS)*cmath.sqrt(cmath.pi)',
              texname = 'G')

aX = Parameter(name = 'aX',
               nature = 'internal',
               type = 'real',
               value = '1/aXM1',
               texname = '\\text{aX}')

MZ = Parameter(name = 'MZ',
               nature = 'internal',
               type = 'real',
               value = 'mZinput',
               texname = '\\text{MZ}')

MZp = Parameter(name = 'MZp',
                nature = 'internal',
                type = 'real',
                value = 'mZDinput',
                texname = '\\text{MZp}')

MH = Parameter(name = 'MH',
               nature = 'internal',
               type = 'real',
               value = 'MHinput',
               texname = '\\text{MH}')

MHS = Parameter(name = 'MHS',
                nature = 'internal',
                type = 'real',
                value = 'MHSinput',
                texname = '\\text{MHS}')

v = Parameter(name = 'v',
              nature = 'internal',
              type = 'real',
              value = '1/(2**0.25*cmath.sqrt(Gf))',
              texname = 'v')

CKM1x1 = Parameter(name = 'CKM1x1',
                   nature = 'internal',
                   type = 'complex',
                   value = 'cmath.cos(cabi)',
                   texname = '\\text{CKM1x1}')

CKM1x2 = Parameter(name = 'CKM1x2',
                   nature = 'internal',
                   type = 'complex',
                   value = 'cmath.sin(cabi)',
                   texname = '\\text{CKM1x2}')

CKM2x1 = Parameter(name = 'CKM2x1',
                   nature = 'internal',
                   type = 'complex',
                   value = '-cmath.sin(cabi)',
                   texname = '\\text{CKM2x1}')

CKM2x2 = Parameter(name = 'CKM2x2',
                   nature = 'internal',
                   type = 'complex',
                   value = 'cmath.cos(cabi)',
                   texname = '\\text{CKM2x2}')

eta = Parameter(name = 'eta',
                nature = 'internal',
                type = 'real',
                value = 'epsilon/(cw*cmath.sqrt(1 - epsilon**2/cw**2))',
                texname = '\\eta')

ee = Parameter(name = 'ee',
               nature = 'internal',
               type = 'real',
               value = '2*cmath.sqrt(aEW)*cmath.sqrt(cmath.pi)',
               texname = 'e')

GH = Parameter(name = 'GH',
               nature = 'internal',
               type = 'real',
               value = '-(G**2*(1 + (13*MH**6)/(16800.*MT**6) + MH**4/(168.*MT**4) + (7*MH**2)/(120.*MT**2)))/(12.*cmath.pi**2*v)',
               texname = 'G_H')

Gphi = Parameter(name = 'Gphi',
                 nature = 'internal',
                 type = 'real',
                 value = '-(G**2*(1 + MH**6/(560.*MT**6) + MH**4/(90.*MT**4) + MH**2/(12.*MT**2)))/(8.*cmath.pi**2*v)',
                 texname = 'G_h')

gX = Parameter(name = 'gX',
               nature = 'internal',
               type = 'real',
               value = '2*cmath.sqrt(aX)*cmath.sqrt(cmath.pi)',
               texname = 'g_X')

yb = Parameter(name = 'yb',
               nature = 'internal',
               type = 'real',
               value = '(ymb*cmath.sqrt(2))/v',
               texname = '\\text{yb}')

yc = Parameter(name = 'yc',
               nature = 'internal',
               type = 'real',
               value = '(ymc*cmath.sqrt(2))/v',
               texname = '\\text{yc}')

ye = Parameter(name = 'ye',
               nature = 'internal',
               type = 'real',
               value = '(ymel*cmath.sqrt(2))/v',
               texname = '\\text{ye}')

ym = Parameter(name = 'ym',
               nature = 'internal',
               type = 'real',
               value = '(ymmu*cmath.sqrt(2))/v',
               texname = '\\text{ym}')

yt = Parameter(name = 'yt',
               nature = 'internal',
               type = 'real',
               value = '(ymt*cmath.sqrt(2))/v',
               texname = '\\text{yt}')

ytau = Parameter(name = 'ytau',
                 nature = 'internal',
                 type = 'real',
                 value = '(ymtau*cmath.sqrt(2))/v',
                 texname = '\\text{ytau}')

chi = Parameter(name = 'chi',
                nature = 'internal',
                type = 'real',
                value = 'eta/cmath.sqrt(1 + eta**2)',
                texname = '\\chi')

DZ = Parameter(name = 'DZ',
               nature = 'internal',
               type = 'real',
               value = '(mZDinput**4 + mZinput**4 + (2*mZinput**2*cmath.atan(10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000*(mZDinput - mZinput))*cmath.sqrt(mZDinput**4 + mZinput**4 - 2*mZDinput**2*mZinput**2*(1 + 2*eta**2*sw**2)))/cmath.pi + mZDinput**2*(-2*eta**2*mZinput**2*sw**2 + (2*cmath.atan(10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000*(mZDinput - mZinput))*cmath.sqrt(mZDinput**4 + mZinput**4 - 2*mZDinput**2*mZinput**2*(1 + 2*eta**2*sw**2)))/cmath.pi))/(2.*mZDinput**2*mZinput**2)',
               texname = '\\text{DZ}')

MZ0 = Parameter(name = 'MZ0',
                nature = 'internal',
                type = 'real',
                value = 'cmath.sqrt((mZDinput**2 + mZinput**2 - (2*cmath.atan(10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000*(mZDinput - mZinput))*cmath.sqrt(mZDinput**4 + mZinput**4 - 2*mZDinput**2*mZinput**2*(1 + 2*eta**2*sw**2)))/cmath.pi)/(2 + 2*eta**2*sw**2))',
                texname = '\\text{MZ0}')

g1 = Parameter(name = 'g1',
               nature = 'internal',
               type = 'real',
               value = 'ee/cw',
               texname = 'g_1')

gw = Parameter(name = 'gw',
               nature = 'internal',
               type = 'real',
               value = 'ee/sw',
               texname = 'g_w')

MW = Parameter(name = 'MW',
               nature = 'internal',
               type = 'real',
               value = 'cw*MZ0',
               texname = 'M_W')

MX = Parameter(name = 'MX',
               nature = 'internal',
               type = 'real',
               value = 'MZ0*cmath.sqrt(DZ)',
               texname = '\\text{MX}')

ta = Parameter(name = 'ta',
               nature = 'internal',
               type = 'real',
               value = '-(-1 + DZ + eta**2*sw**2 - (2*cmath.atan(10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000*(-1 + DZ))*cmath.sqrt(4*eta**2*sw**2 + (-1 + DZ + eta**2*sw**2)**2))/cmath.pi)/(2.*eta*sw)',
               texname = 't_{\\alpha }')

ca = Parameter(name = 'ca',
               nature = 'internal',
               type = 'real',
               value = '1/cmath.sqrt(1 + ta**2)',
               texname = 'c_{\\alpha }')

sa = Parameter(name = 'sa',
               nature = 'internal',
               type = 'real',
               value = 'ta/cmath.sqrt(1 + ta**2)',
               texname = 's_{\\alpha }')

AH = Parameter(name = 'AH',
               nature = 'internal',
               type = 'real',
               value = '(47*ee**2*(1 - (2*MH**4)/(987.*MT**4) - (14*MH**2)/(705.*MT**2) + (213*MH**12)/(2.634632e7*MW**12) + (5*MH**10)/(119756.*MW**10) + (41*MH**8)/(180950.*MW**8) + (87*MH**6)/(65800.*MW**6) + (57*MH**4)/(6580.*MW**4) + (33*MH**2)/(470.*MW**2)))/(72.*cmath.pi**2*v)',
               texname = 'A_H')

xi = Parameter(name = 'xi',
               nature = 'internal',
               type = 'real',
               value = '(MX*cmath.sqrt(1 - chi**2))/gX',
               texname = '\\xi')

th = Parameter(name = 'th',
               nature = 'internal',
               type = 'real',
               value = '(-MHinput**2 + MHSinput**2 + (2*cmath.atan(10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000*(MHinput - MHSinput))*cmath.sqrt((MHinput**2 - MHSinput**2)**2 - 4*kap**2*v**2*xi**2))/cmath.pi)/(2.*kap*v*xi)',
               texname = 't_h')

lam = Parameter(name = 'lam',
                nature = 'internal',
                type = 'real',
                value = '(MHinput**2 + MHSinput**2 + (2*cmath.atan(10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000*(MHinput - MHSinput))*cmath.sqrt((MHinput**2 - MHSinput**2)**2 - 4*kap**2*v**2*xi**2))/cmath.pi)/(4.*v**2)',
                texname = '\\text{lam}')

rho = Parameter(name = 'rho',
                nature = 'internal',
                type = 'real',
                value = '(MHinput**2 + MHSinput**2 - (2*cmath.atan(10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000*(MHinput - MHSinput))*cmath.sqrt((MHinput**2 - MHSinput**2)**2 - 4*kap**2*v**2*xi**2))/cmath.pi)/(4.*xi**2)',
                texname = '\\rho')

ch = Parameter(name = 'ch',
               nature = 'internal',
               type = 'real',
               value = '1/cmath.sqrt(1 + th**2)',
               texname = 'c_h')

muH2 = Parameter(name = 'muH2',
                 nature = 'internal',
                 type = 'real',
                 value = '(kap*v**2)/2. + rho*xi**2',
                 texname = '\\text{muH2}')

muSM2 = Parameter(name = 'muSM2',
                  nature = 'internal',
                  type = 'real',
                  value = 'lam*v**2 + (kap*xi**2)/2.',
                  texname = '\\text{muSM2}')

sh = Parameter(name = 'sh',
               nature = 'internal',
               type = 'real',
               value = 'th/cmath.sqrt(1 + th**2)',
               texname = 's_h')

