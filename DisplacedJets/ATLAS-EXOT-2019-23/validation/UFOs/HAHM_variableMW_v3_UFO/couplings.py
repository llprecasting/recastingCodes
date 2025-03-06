# This file was automatically created by FeynRules 2.0.23
# Mathematica version: 9.0 for Mac OS X x86 (64-bit) (November 20, 2012)
# Date: Sat 20 Sep 2014 16:11:37


from object_library import all_couplings, Coupling

from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot



GC_1 = Coupling(name = 'GC_1',
                value = '-(AH*ch*complex(0,1))',
                order = {'HIW':1})

GC_2 = Coupling(name = 'GC_2',
                value = '-(ee*complex(0,1))/3.',
                order = {'QED':1})

GC_3 = Coupling(name = 'GC_3',
                value = '(2*ee*complex(0,1))/3.',
                order = {'QED':1})

GC_4 = Coupling(name = 'GC_4',
                value = '-(ee*complex(0,1))',
                order = {'QED':1})

GC_5 = Coupling(name = 'GC_5',
                value = '-G',
                order = {'QCD':1})

GC_6 = Coupling(name = 'GC_6',
                value = 'complex(0,1)*G',
                order = {'QCD':1})

GC_7 = Coupling(name = 'GC_7',
                value = 'complex(0,1)*G**2',
                order = {'QCD':2})

GC_8 = Coupling(name = 'GC_8',
                value = '-(ch*complex(0,1)*GH)',
                order = {'HIG':1})

GC_9 = Coupling(name = 'GC_9',
                value = '-(ch*G*GH)',
                order = {'HIG':1,'QCD':1})

GC_10 = Coupling(name = 'GC_10',
                 value = 'ch*complex(0,1)*G**2*GH',
                 order = {'HIG':1,'QCD':2})

GC_11 = Coupling(name = 'GC_11',
                 value = 'ca*cw*complex(0,1)*gw',
                 order = {'QED':1})

GC_12 = Coupling(name = 'GC_12',
                 value = '-(complex(0,1)*gw**2)',
                 order = {'QED':2})

GC_13 = Coupling(name = 'GC_13',
                 value = 'ca**2*cw**2*complex(0,1)*gw**2',
                 order = {'QED':2})

GC_14 = Coupling(name = 'GC_14',
                 value = '-(cw*complex(0,1)*gw*sa)',
                 order = {'QED':1})

GC_15 = Coupling(name = 'GC_15',
                 value = '-(ca*cw**2*complex(0,1)*gw**2*sa)',
                 order = {'QED':2})

GC_16 = Coupling(name = 'GC_16',
                 value = 'cw**2*complex(0,1)*gw**2*sa**2',
                 order = {'QED':2})

GC_17 = Coupling(name = 'GC_17',
                 value = '-(AH*complex(0,1)*sh)',
                 order = {'HIW':1})

GC_18 = Coupling(name = 'GC_18',
                 value = '-(complex(0,1)*GH*sh)',
                 order = {'HIG':1})

GC_19 = Coupling(name = 'GC_19',
                 value = '-(G*GH*sh)',
                 order = {'HIG':1,'QCD':1})

GC_20 = Coupling(name = 'GC_20',
                 value = 'complex(0,1)*G**2*GH*sh',
                 order = {'HIG':1,'QCD':2})

GC_21 = Coupling(name = 'GC_21',
                 value = '-3*ch**3*complex(0,1)*kap*sh + 6*ch**3*complex(0,1)*rho*sh + 3*ch*complex(0,1)*kap*sh**3 - 6*ch*complex(0,1)*lam*sh**3',
                 order = {'QED':2})

GC_22 = Coupling(name = 'GC_22',
                 value = '3*ch**3*complex(0,1)*kap*sh - 6*ch**3*complex(0,1)*lam*sh - 3*ch*complex(0,1)*kap*sh**3 + 6*ch*complex(0,1)*rho*sh**3',
                 order = {'QED':2})

GC_23 = Coupling(name = 'GC_23',
                 value = '-(ch**4*complex(0,1)*kap) + 4*ch**2*complex(0,1)*kap*sh**2 - 6*ch**2*complex(0,1)*lam*sh**2 - 6*ch**2*complex(0,1)*rho*sh**2 - complex(0,1)*kap*sh**4',
                 order = {'QED':2})

GC_24 = Coupling(name = 'GC_24',
                 value = '-6*ch**4*complex(0,1)*rho - 6*ch**2*complex(0,1)*kap*sh**2 - 6*complex(0,1)*lam*sh**4',
                 order = {'QED':2})

GC_25 = Coupling(name = 'GC_25',
                 value = '-6*ch**4*complex(0,1)*lam - 6*ch**2*complex(0,1)*kap*sh**2 - 6*complex(0,1)*rho*sh**4',
                 order = {'QED':2})

GC_26 = Coupling(name = 'GC_26',
                 value = '(ch**2*ee**2*complex(0,1))/(2.*sw**2)',
                 order = {'QED':2})

GC_27 = Coupling(name = 'GC_27',
                 value = '(ch*ee**2*complex(0,1)*sh)/(2.*sw**2)',
                 order = {'QED':2})

GC_28 = Coupling(name = 'GC_28',
                 value = '(ee**2*complex(0,1)*sh**2)/(2.*sw**2)',
                 order = {'QED':2})

GC_29 = Coupling(name = 'GC_29',
                 value = '(ee*complex(0,1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_30 = Coupling(name = 'GC_30',
                 value = '(CKM1x1*ee*complex(0,1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_31 = Coupling(name = 'GC_31',
                 value = '(CKM1x2*ee*complex(0,1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_32 = Coupling(name = 'GC_32',
                 value = '(CKM2x1*ee*complex(0,1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_33 = Coupling(name = 'GC_33',
                 value = '(CKM2x2*ee*complex(0,1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_34 = Coupling(name = 'GC_34',
                 value = '-(ca*cw*ee*complex(0,1))/(2.*sw)',
                 order = {'QED':1})

GC_35 = Coupling(name = 'GC_35',
                 value = '(ca*cw*ee*complex(0,1))/(2.*sw)',
                 order = {'QED':1})

GC_36 = Coupling(name = 'GC_36',
                 value = '-(cw*ee*complex(0,1)*sa)/(2.*sw)',
                 order = {'QED':1})

GC_37 = Coupling(name = 'GC_37',
                 value = '(cw*ee*complex(0,1)*sa)/(2.*sw)',
                 order = {'QED':1})

GC_38 = Coupling(name = 'GC_38',
                 value = 'complex(0,1)*gw*sw',
                 order = {'QED':1})

GC_39 = Coupling(name = 'GC_39',
                 value = '-2*ca*cw*complex(0,1)*gw**2*sw',
                 order = {'QED':2})

GC_40 = Coupling(name = 'GC_40',
                 value = '2*cw*complex(0,1)*gw**2*sa*sw',
                 order = {'QED':2})

GC_41 = Coupling(name = 'GC_41',
                 value = 'complex(0,1)*gw**2*sw**2',
                 order = {'QED':2})

GC_42 = Coupling(name = 'GC_42',
                 value = '(ee*eta*complex(0,1)*sa)/(6.*cw) - (ca*ee*complex(0,1)*sw)/(6.*cw)',
                 order = {'QED':1})

GC_43 = Coupling(name = 'GC_43',
                 value = '-(ee*eta*complex(0,1)*sa)/(2.*cw) + (ca*ee*complex(0,1)*sw)/(2.*cw)',
                 order = {'QED':1})

GC_44 = Coupling(name = 'GC_44',
                 value = '-(ee*eta*complex(0,1)*sa)/(2.*cw) + (ca*cw*ee*complex(0,1))/(2.*sw) + (ca*ee*complex(0,1)*sw)/(2.*cw)',
                 order = {'QED':1})

GC_45 = Coupling(name = 'GC_45',
                 value = '(ca*ee*eta*complex(0,1))/(6.*cw) + (ee*complex(0,1)*sa*sw)/(6.*cw)',
                 order = {'QED':1})

GC_46 = Coupling(name = 'GC_46',
                 value = '-(ca*ee*eta*complex(0,1))/(2.*cw) - (ee*complex(0,1)*sa*sw)/(2.*cw)',
                 order = {'QED':1})

GC_47 = Coupling(name = 'GC_47',
                 value = '-(ca*ee*eta*complex(0,1))/(2.*cw) - (cw*ee*complex(0,1)*sa)/(2.*sw) - (ee*complex(0,1)*sa*sw)/(2.*cw)',
                 order = {'QED':1})

GC_48 = Coupling(name = 'GC_48',
                 value = 'ca**2*ch**2*ee**2*complex(0,1) + (ch**2*ee**2*eta**2*complex(0,1)*sa**2)/(2.*cw**2) + (2*eta**2*complex(0,1)*gX**2*sa**2*sh**2)/chi**2 + (ca**2*ch**2*cw**2*ee**2*complex(0,1))/(2.*sw**2) - (ca*ch**2*ee**2*eta*complex(0,1)*sa)/sw - (ca*ch**2*ee**2*eta*complex(0,1)*sa*sw)/cw**2 + (ca**2*ch**2*ee**2*complex(0,1)*sw**2)/(2.*cw**2)',
                 order = {'QED':2})

GC_49 = Coupling(name = 'GC_49',
                 value = '-(ca*ch**2*ee**2*complex(0,1)*sa) + (ca*ch**2*ee**2*eta**2*complex(0,1)*sa)/(2.*cw**2) + (2*ca*eta**2*complex(0,1)*gX**2*sa*sh**2)/chi**2 - (ca*ch**2*cw**2*ee**2*complex(0,1)*sa)/(2.*sw**2) - (ca**2*ch**2*ee**2*eta*complex(0,1))/(2.*sw) + (ch**2*ee**2*eta*complex(0,1)*sa**2)/(2.*sw) - (ca**2*ch**2*ee**2*eta*complex(0,1)*sw)/(2.*cw**2) + (ch**2*ee**2*eta*complex(0,1)*sa**2*sw)/(2.*cw**2) - (ca*ch**2*ee**2*complex(0,1)*sa*sw**2)/(2.*cw**2)',
                 order = {'QED':2})

GC_50 = Coupling(name = 'GC_50',
                 value = '(ca**2*ch**2*ee**2*eta**2*complex(0,1))/(2.*cw**2) + ch**2*ee**2*complex(0,1)*sa**2 + (2*ca**2*eta**2*complex(0,1)*gX**2*sh**2)/chi**2 + (ch**2*cw**2*ee**2*complex(0,1)*sa**2)/(2.*sw**2) + (ca*ch**2*ee**2*eta*complex(0,1)*sa)/sw + (ca*ch**2*ee**2*eta*complex(0,1)*sa*sw)/cw**2 + (ch**2*ee**2*complex(0,1)*sa**2*sw**2)/(2.*cw**2)',
                 order = {'QED':2})

GC_51 = Coupling(name = 'GC_51',
                 value = 'ca**2*ch*ee**2*complex(0,1)*sh + (ch*ee**2*eta**2*complex(0,1)*sa**2*sh)/(2.*cw**2) - (2*ch*eta**2*complex(0,1)*gX**2*sa**2*sh)/chi**2 + (ca**2*ch*cw**2*ee**2*complex(0,1)*sh)/(2.*sw**2) - (ca*ch*ee**2*eta*complex(0,1)*sa*sh)/sw - (ca*ch*ee**2*eta*complex(0,1)*sa*sh*sw)/cw**2 + (ca**2*ch*ee**2*complex(0,1)*sh*sw**2)/(2.*cw**2)',
                 order = {'QED':2})

GC_52 = Coupling(name = 'GC_52',
                 value = '-(ca*ch*ee**2*complex(0,1)*sa*sh) + (ca*ch*ee**2*eta**2*complex(0,1)*sa*sh)/(2.*cw**2) - (2*ca*ch*eta**2*complex(0,1)*gX**2*sa*sh)/chi**2 - (ca*ch*cw**2*ee**2*complex(0,1)*sa*sh)/(2.*sw**2) - (ca**2*ch*ee**2*eta*complex(0,1)*sh)/(2.*sw) + (ch*ee**2*eta*complex(0,1)*sa**2*sh)/(2.*sw) - (ca**2*ch*ee**2*eta*complex(0,1)*sh*sw)/(2.*cw**2) + (ch*ee**2*eta*complex(0,1)*sa**2*sh*sw)/(2.*cw**2) - (ca*ch*ee**2*complex(0,1)*sa*sh*sw**2)/(2.*cw**2)',
                 order = {'QED':2})

GC_53 = Coupling(name = 'GC_53',
                 value = '(ca**2*ch*ee**2*eta**2*complex(0,1)*sh)/(2.*cw**2) - (2*ca**2*ch*eta**2*complex(0,1)*gX**2*sh)/chi**2 + ch*ee**2*complex(0,1)*sa**2*sh + (ch*cw**2*ee**2*complex(0,1)*sa**2*sh)/(2.*sw**2) + (ca*ch*ee**2*eta*complex(0,1)*sa*sh)/sw + (ca*ch*ee**2*eta*complex(0,1)*sa*sh*sw)/cw**2 + (ch*ee**2*complex(0,1)*sa**2*sh*sw**2)/(2.*cw**2)',
                 order = {'QED':2})

GC_54 = Coupling(name = 'GC_54',
                 value = '(2*ch**2*eta**2*complex(0,1)*gX**2*sa**2)/chi**2 + ca**2*ee**2*complex(0,1)*sh**2 + (ee**2*eta**2*complex(0,1)*sa**2*sh**2)/(2.*cw**2) + (ca**2*cw**2*ee**2*complex(0,1)*sh**2)/(2.*sw**2) - (ca*ee**2*eta*complex(0,1)*sa*sh**2)/sw - (ca*ee**2*eta*complex(0,1)*sa*sh**2*sw)/cw**2 + (ca**2*ee**2*complex(0,1)*sh**2*sw**2)/(2.*cw**2)',
                 order = {'QED':2})

GC_55 = Coupling(name = 'GC_55',
                 value = '(2*ca*ch**2*eta**2*complex(0,1)*gX**2*sa)/chi**2 - ca*ee**2*complex(0,1)*sa*sh**2 + (ca*ee**2*eta**2*complex(0,1)*sa*sh**2)/(2.*cw**2) - (ca*cw**2*ee**2*complex(0,1)*sa*sh**2)/(2.*sw**2) - (ca**2*ee**2*eta*complex(0,1)*sh**2)/(2.*sw) + (ee**2*eta*complex(0,1)*sa**2*sh**2)/(2.*sw) - (ca**2*ee**2*eta*complex(0,1)*sh**2*sw)/(2.*cw**2) + (ee**2*eta*complex(0,1)*sa**2*sh**2*sw)/(2.*cw**2) - (ca*ee**2*complex(0,1)*sa*sh**2*sw**2)/(2.*cw**2)',
                 order = {'QED':2})

GC_56 = Coupling(name = 'GC_56',
                 value = '(2*ca**2*ch**2*eta**2*complex(0,1)*gX**2)/chi**2 + (ca**2*ee**2*eta**2*complex(0,1)*sh**2)/(2.*cw**2) + ee**2*complex(0,1)*sa**2*sh**2 + (cw**2*ee**2*complex(0,1)*sa**2*sh**2)/(2.*sw**2) + (ca*ee**2*eta*complex(0,1)*sa*sh**2)/sw + (ca*ee**2*eta*complex(0,1)*sa*sh**2*sw)/cw**2 + (ee**2*complex(0,1)*sa**2*sh**2*sw**2)/(2.*cw**2)',
                 order = {'QED':2})

GC_57 = Coupling(name = 'GC_57',
                 value = '(ch*ee**2*complex(0,1)*v)/(2.*sw**2)',
                 order = {'QED':1})

GC_58 = Coupling(name = 'GC_58',
                 value = '(ee**2*complex(0,1)*sh*v)/(2.*sw**2)',
                 order = {'QED':1})

GC_59 = Coupling(name = 'GC_59',
                 value = '(ca**2*ee**2*eta**2*complex(0,1)*sh*v)/(2.*cw**2) + ee**2*complex(0,1)*sa**2*sh*v + (cw**2*ee**2*complex(0,1)*sa**2*sh*v)/(2.*sw**2) + (ca*ee**2*eta*complex(0,1)*sa*sh*v)/sw + (ca*ee**2*eta*complex(0,1)*sa*sh*sw*v)/cw**2 + (ee**2*complex(0,1)*sa**2*sh*sw**2*v)/(2.*cw**2) + (2*ca**2*ch*eta**2*complex(0,1)*gX**2*xi)/chi**2',
                 order = {'QED':1})

GC_60 = Coupling(name = 'GC_60',
                 value = '-(ca*ee**2*complex(0,1)*sa*sh*v) + (ca*ee**2*eta**2*complex(0,1)*sa*sh*v)/(2.*cw**2) - (ca*cw**2*ee**2*complex(0,1)*sa*sh*v)/(2.*sw**2) - (ca**2*ee**2*eta*complex(0,1)*sh*v)/(2.*sw) + (ee**2*eta*complex(0,1)*sa**2*sh*v)/(2.*sw) - (ca**2*ee**2*eta*complex(0,1)*sh*sw*v)/(2.*cw**2) + (ee**2*eta*complex(0,1)*sa**2*sh*sw*v)/(2.*cw**2) - (ca*ee**2*complex(0,1)*sa*sh*sw**2*v)/(2.*cw**2) + (2*ca*ch*eta**2*complex(0,1)*gX**2*sa*xi)/chi**2',
                 order = {'QED':1})

GC_61 = Coupling(name = 'GC_61',
                 value = 'ca**2*ee**2*complex(0,1)*sh*v + (ee**2*eta**2*complex(0,1)*sa**2*sh*v)/(2.*cw**2) + (ca**2*cw**2*ee**2*complex(0,1)*sh*v)/(2.*sw**2) - (ca*ee**2*eta*complex(0,1)*sa*sh*v)/sw - (ca*ee**2*eta*complex(0,1)*sa*sh*sw*v)/cw**2 + (ca**2*ee**2*complex(0,1)*sh*sw**2*v)/(2.*cw**2) + (2*ch*eta**2*complex(0,1)*gX**2*sa**2*xi)/chi**2',
                 order = {'QED':1})

GC_62 = Coupling(name = 'GC_62',
                 value = '(ca**2*ch*ee**2*eta**2*complex(0,1)*v)/(2.*cw**2) + ch*ee**2*complex(0,1)*sa**2*v + (ch*cw**2*ee**2*complex(0,1)*sa**2*v)/(2.*sw**2) + (ca*ch*ee**2*eta*complex(0,1)*sa*v)/sw + (ca*ch*ee**2*eta*complex(0,1)*sa*sw*v)/cw**2 + (ch*ee**2*complex(0,1)*sa**2*sw**2*v)/(2.*cw**2) - (2*ca**2*eta**2*complex(0,1)*gX**2*sh*xi)/chi**2',
                 order = {'QED':1})

GC_63 = Coupling(name = 'GC_63',
                 value = '-(ca*ch*ee**2*complex(0,1)*sa*v) + (ca*ch*ee**2*eta**2*complex(0,1)*sa*v)/(2.*cw**2) - (ca*ch*cw**2*ee**2*complex(0,1)*sa*v)/(2.*sw**2) - (ca**2*ch*ee**2*eta*complex(0,1)*v)/(2.*sw) + (ch*ee**2*eta*complex(0,1)*sa**2*v)/(2.*sw) - (ca**2*ch*ee**2*eta*complex(0,1)*sw*v)/(2.*cw**2) + (ch*ee**2*eta*complex(0,1)*sa**2*sw*v)/(2.*cw**2) - (ca*ch*ee**2*complex(0,1)*sa*sw**2*v)/(2.*cw**2) - (2*ca*eta**2*complex(0,1)*gX**2*sa*sh*xi)/chi**2',
                 order = {'QED':1})

GC_64 = Coupling(name = 'GC_64',
                 value = 'ca**2*ch*ee**2*complex(0,1)*v + (ch*ee**2*eta**2*complex(0,1)*sa**2*v)/(2.*cw**2) + (ca**2*ch*cw**2*ee**2*complex(0,1)*v)/(2.*sw**2) - (ca*ch*ee**2*eta*complex(0,1)*sa*v)/sw - (ca*ch*ee**2*eta*complex(0,1)*sa*sw*v)/cw**2 + (ca**2*ch*ee**2*complex(0,1)*sw**2*v)/(2.*cw**2) - (2*eta**2*complex(0,1)*gX**2*sa**2*sh*xi)/chi**2',
                 order = {'QED':1})

GC_65 = Coupling(name = 'GC_65',
                 value = '-3*ch**2*complex(0,1)*kap*sh*v - 6*complex(0,1)*lam*sh**3*v - 6*ch**3*complex(0,1)*rho*xi - 3*ch*complex(0,1)*kap*sh**2*xi',
                 order = {'QED':1})

GC_66 = Coupling(name = 'GC_66',
                 value = '2*ch**2*complex(0,1)*kap*sh*v - 6*ch**2*complex(0,1)*lam*sh*v - complex(0,1)*kap*sh**3*v - ch**3*complex(0,1)*kap*xi + 2*ch*complex(0,1)*kap*sh**2*xi - 6*ch*complex(0,1)*rho*sh**2*xi',
                 order = {'QED':1})

GC_67 = Coupling(name = 'GC_67',
                 value = '-(ch**3*complex(0,1)*kap*v) + 2*ch*complex(0,1)*kap*sh**2*v - 6*ch*complex(0,1)*lam*sh**2*v - 2*ch**2*complex(0,1)*kap*sh*xi + 6*ch**2*complex(0,1)*rho*sh*xi + complex(0,1)*kap*sh**3*xi',
                 order = {'QED':1})

GC_68 = Coupling(name = 'GC_68',
                 value = '-6*ch**3*complex(0,1)*lam*v - 3*ch*complex(0,1)*kap*sh**2*v + 3*ch**2*complex(0,1)*kap*sh*xi + 6*complex(0,1)*rho*sh**3*xi',
                 order = {'QED':1})

GC_69 = Coupling(name = 'GC_69',
                 value = '-((ch*complex(0,1)*yb)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_70 = Coupling(name = 'GC_70',
                 value = '-((complex(0,1)*sh*yb)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_71 = Coupling(name = 'GC_71',
                 value = '-((ch*complex(0,1)*yc)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_72 = Coupling(name = 'GC_72',
                 value = '-((complex(0,1)*sh*yc)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_73 = Coupling(name = 'GC_73',
                 value = '-((ch*complex(0,1)*ye)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_74 = Coupling(name = 'GC_74',
                 value = '-((complex(0,1)*sh*ye)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_75 = Coupling(name = 'GC_75',
                 value = '-((ch*complex(0,1)*ym)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_76 = Coupling(name = 'GC_76',
                 value = '-((complex(0,1)*sh*ym)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_77 = Coupling(name = 'GC_77',
                 value = '-((ch*complex(0,1)*yt)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_78 = Coupling(name = 'GC_78',
                 value = '-((complex(0,1)*sh*yt)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_79 = Coupling(name = 'GC_79',
                 value = '-((ch*complex(0,1)*ytau)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_80 = Coupling(name = 'GC_80',
                 value = '-((complex(0,1)*sh*ytau)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_81 = Coupling(name = 'GC_81',
                 value = '(ee*complex(0,1)*complexconjugate(CKM1x1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_82 = Coupling(name = 'GC_82',
                 value = '(ee*complex(0,1)*complexconjugate(CKM1x2))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_83 = Coupling(name = 'GC_83',
                 value = '(ee*complex(0,1)*complexconjugate(CKM2x1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_84 = Coupling(name = 'GC_84',
                 value = '(ee*complex(0,1)*complexconjugate(CKM2x2))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

