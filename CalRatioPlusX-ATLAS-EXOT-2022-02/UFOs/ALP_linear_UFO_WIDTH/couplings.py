# This file was automatically created by FeynRules =.2.4
# Mathematica version: 10.0 for Linux x86 (64-bit) (December 4, 2014)
# Date: Fri 9 Mar 2018 12:00:41


from object_library import all_couplings, Coupling

from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot



GC_1 = Coupling(name = 'GC_1',
                value = '-(ee*complex(0,1))/3.',
                order = {'QED':1})

GC_2 = Coupling(name = 'GC_2',
                value = '(2*ee*complex(0,1))/3.',
                order = {'QED':1})

GC_3 = Coupling(name = 'GC_3',
                value = '-(ee*complex(0,1))',
                order = {'QED':1})

GC_4 = Coupling(name = 'GC_4',
                value = 'ee*complex(0,1)',
                order = {'QED':1})

GC_5 = Coupling(name = 'GC_5',
                value = 'ee**2*complex(0,1)',
                order = {'QED':2})

GC_6 = Coupling(name = 'GC_6',
                value = '(CGtil*complex(0,1))/(2.*fa)',
                order = {'NP':1})

GC_7 = Coupling(name = 'GC_7',
                value = '(2*CWtil*complex(0,1))/fa',
                order = {'NP':1})

GC_8 = Coupling(name = 'GC_8',
                value = '(4*CWtil*ee*complex(0,1))/fa',
                order = {'NP':1,'QED':1})

GC_9 = Coupling(name = 'GC_9',
                value = '-G',
                order = {'QCD':1})

GC_10 = Coupling(name = 'GC_10',
                 value = 'complex(0,1)*G',
                 order = {'QCD':1})

GC_11 = Coupling(name = 'GC_11',
                 value = '(CGtil*G)/fa',
                 order = {'NP':1,'QCD':1})

GC_12 = Coupling(name = 'GC_12',
                 value = 'complex(0,1)*G**2',
                 order = {'QCD':2})

GC_13 = Coupling(name = 'GC_13',
                 value = '-6*complex(0,1)*lam',
                 order = {'QED':2})

GC_14 = Coupling(name = 'GC_14',
                 value = '(ee**2*complex(0,1))/(2.*sw**2)',
                 order = {'QED':2})

GC_15 = Coupling(name = 'GC_15',
                 value = '-((ee**2*complex(0,1))/sw**2)',
                 order = {'QED':2})

GC_16 = Coupling(name = 'GC_16',
                 value = '(cw**2*ee**2*complex(0,1))/sw**2',
                 order = {'QED':2})

GC_17 = Coupling(name = 'GC_17',
                 value = '(ee*complex(0,1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_18 = Coupling(name = 'GC_18',
                 value = '(CKM1x1*ee*complex(0,1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_19 = Coupling(name = 'GC_19',
                 value = '(CKM1x2*ee*complex(0,1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_20 = Coupling(name = 'GC_20',
                 value = '(CKM1x3*ee*complex(0,1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_21 = Coupling(name = 'GC_21',
                 value = '(CKM2x1*ee*complex(0,1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_22 = Coupling(name = 'GC_22',
                 value = '(CKM2x2*ee*complex(0,1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_23 = Coupling(name = 'GC_23',
                 value = '(CKM2x3*ee*complex(0,1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_24 = Coupling(name = 'GC_24',
                 value = '(CKM3x1*ee*complex(0,1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_25 = Coupling(name = 'GC_25',
                 value = '(CKM3x2*ee*complex(0,1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_26 = Coupling(name = 'GC_26',
                 value = '(CKM3x3*ee*complex(0,1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_27 = Coupling(name = 'GC_27',
                 value = '(cw*ee*complex(0,1))/sw',
                 order = {'QED':1})

GC_28 = Coupling(name = 'GC_28',
                 value = '(-2*cw*ee**2*complex(0,1))/sw',
                 order = {'QED':2})

GC_29 = Coupling(name = 'GC_29',
                 value = '(4*cw*CWtil*ee*complex(0,1))/(fa*sw)',
                 order = {'NP':1,'QED':1})

GC_30 = Coupling(name = 'GC_30',
                 value = '(ee*complex(0,1)*sw)/(3.*cw)',
                 order = {'QED':1})

GC_31 = Coupling(name = 'GC_31',
                 value = '(-2*ee*complex(0,1)*sw)/(3.*cw)',
                 order = {'QED':1})

GC_32 = Coupling(name = 'GC_32',
                 value = '(ee*complex(0,1)*sw)/cw',
                 order = {'QED':1})

GC_33 = Coupling(name = 'GC_33',
                 value = '-(cw*ee*complex(0,1))/(2.*sw) - (ee*complex(0,1)*sw)/(6.*cw)',
                 order = {'QED':1})

GC_34 = Coupling(name = 'GC_34',
                 value = '(cw*ee*complex(0,1))/(2.*sw) - (ee*complex(0,1)*sw)/(6.*cw)',
                 order = {'QED':1})

GC_35 = Coupling(name = 'GC_35',
                 value = '-(cw*ee*complex(0,1))/(2.*sw) + (ee*complex(0,1)*sw)/(2.*cw)',
                 order = {'QED':1})

GC_36 = Coupling(name = 'GC_36',
                 value = '(cw*ee*complex(0,1))/(2.*sw) + (ee*complex(0,1)*sw)/(2.*cw)',
                 order = {'QED':1})

GC_37 = Coupling(name = 'GC_37',
                 value = '(-2*CBtil*cw*complex(0,1)*sw)/fa + (2*cw*CWtil*complex(0,1)*sw)/fa',
                 order = {'NP':1})

GC_38 = Coupling(name = 'GC_38',
                 value = 'ee**2*complex(0,1) + (cw**2*ee**2*complex(0,1))/(2.*sw**2) + (ee**2*complex(0,1)*sw**2)/(2.*cw**2)',
                 order = {'QED':2})

GC_39 = Coupling(name = 'GC_39',
                 value = '(2*cw**2*CWtil*complex(0,1))/fa + (2*CBtil*complex(0,1)*sw**2)/fa',
                 order = {'NP':1})

GC_40 = Coupling(name = 'GC_40',
                 value = '(2*CBtil*cw**2*complex(0,1))/fa + (2*CWtil*complex(0,1)*sw**2)/fa',
                 order = {'NP':1})

GC_41 = Coupling(name = 'GC_41',
                 value = '-6*complex(0,1)*lam*vev',
                 order = {'QED':1})

GC_42 = Coupling(name = 'GC_42',
                 value = '(ee**2*complex(0,1)*vev)/(2.*sw**2)',
                 order = {'QED':1})

GC_43 = Coupling(name = 'GC_43',
                 value = 'ee**2*complex(0,1)*vev + (cw**2*ee**2*complex(0,1)*vev)/(2.*sw**2) + (ee**2*complex(0,1)*sw**2*vev)/(2.*cw**2)',
                 order = {'QED':1})

GC_44 = Coupling(name = 'GC_44',
                 value = '-((complex(0,1)*yb)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_45 = Coupling(name = 'GC_45',
                 value = '-((CaPhi*yb)/(fa*cmath.sqrt(2)))',
                 order = {'NP':1,'QED':1})

GC_46 = Coupling(name = 'GC_46',
                 value = '-((CaPhi*vev*yb)/(fa*cmath.sqrt(2)))',
                 order = {'NP':1})

GC_47 = Coupling(name = 'GC_47',
                 value = '-((complex(0,1)*yc)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_48 = Coupling(name = 'GC_48',
                 value = '(CaPhi*yc)/(fa*cmath.sqrt(2))',
                 order = {'NP':1,'QED':1})

GC_49 = Coupling(name = 'GC_49',
                 value = '(CaPhi*vev*yc)/(fa*cmath.sqrt(2))',
                 order = {'NP':1})

GC_50 = Coupling(name = 'GC_50',
                 value = '-((complex(0,1)*ydo)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_51 = Coupling(name = 'GC_51',
                 value = '-((CaPhi*ydo)/(fa*cmath.sqrt(2)))',
                 order = {'NP':1,'QED':1})

GC_52 = Coupling(name = 'GC_52',
                 value = '-((CaPhi*vev*ydo)/(fa*cmath.sqrt(2)))',
                 order = {'NP':1})

GC_53 = Coupling(name = 'GC_53',
                 value = '-((complex(0,1)*ye)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_54 = Coupling(name = 'GC_54',
                 value = '-((CaPhi*ye)/(fa*cmath.sqrt(2)))',
                 order = {'NP':1,'QED':1})

GC_55 = Coupling(name = 'GC_55',
                 value = '-((CaPhi*vev*ye)/(fa*cmath.sqrt(2)))',
                 order = {'NP':1})

GC_56 = Coupling(name = 'GC_56',
                 value = '-((complex(0,1)*ym)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_57 = Coupling(name = 'GC_57',
                 value = '-((CaPhi*ym)/(fa*cmath.sqrt(2)))',
                 order = {'NP':1,'QED':1})

GC_58 = Coupling(name = 'GC_58',
                 value = '-((CaPhi*vev*ym)/(fa*cmath.sqrt(2)))',
                 order = {'NP':1})

GC_59 = Coupling(name = 'GC_59',
                 value = '-((complex(0,1)*ys)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_60 = Coupling(name = 'GC_60',
                 value = '-((CaPhi*ys)/(fa*cmath.sqrt(2)))',
                 order = {'NP':1,'QED':1})

GC_61 = Coupling(name = 'GC_61',
                 value = '-((CaPhi*vev*ys)/(fa*cmath.sqrt(2)))',
                 order = {'NP':1})

GC_62 = Coupling(name = 'GC_62',
                 value = '-((complex(0,1)*yt)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_63 = Coupling(name = 'GC_63',
                 value = '(CaPhi*yt)/(fa*cmath.sqrt(2))',
                 order = {'NP':1,'QED':1})

GC_64 = Coupling(name = 'GC_64',
                 value = '(CaPhi*vev*yt)/(fa*cmath.sqrt(2))',
                 order = {'NP':1})

GC_65 = Coupling(name = 'GC_65',
                 value = '-((complex(0,1)*ytau)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_66 = Coupling(name = 'GC_66',
                 value = '-((CaPhi*ytau)/(fa*cmath.sqrt(2)))',
                 order = {'NP':1,'QED':1})

GC_67 = Coupling(name = 'GC_67',
                 value = '-((CaPhi*vev*ytau)/(fa*cmath.sqrt(2)))',
                 order = {'NP':1})

GC_68 = Coupling(name = 'GC_68',
                 value = '-((complex(0,1)*yup)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_69 = Coupling(name = 'GC_69',
                 value = '(CaPhi*yup)/(fa*cmath.sqrt(2))',
                 order = {'NP':1,'QED':1})

GC_70 = Coupling(name = 'GC_70',
                 value = '(CaPhi*vev*yup)/(fa*cmath.sqrt(2))',
                 order = {'NP':1})

GC_71 = Coupling(name = 'GC_71',
                 value = '(ee*complex(0,1)*complexconjugate(CKM1x1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_72 = Coupling(name = 'GC_72',
                 value = '(ee*complex(0,1)*complexconjugate(CKM1x2))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_73 = Coupling(name = 'GC_73',
                 value = '(ee*complex(0,1)*complexconjugate(CKM1x3))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_74 = Coupling(name = 'GC_74',
                 value = '(ee*complex(0,1)*complexconjugate(CKM2x1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_75 = Coupling(name = 'GC_75',
                 value = '(ee*complex(0,1)*complexconjugate(CKM2x2))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_76 = Coupling(name = 'GC_76',
                 value = '(ee*complex(0,1)*complexconjugate(CKM2x3))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_77 = Coupling(name = 'GC_77',
                 value = '(ee*complex(0,1)*complexconjugate(CKM3x1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_78 = Coupling(name = 'GC_78',
                 value = '(ee*complex(0,1)*complexconjugate(CKM3x2))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_79 = Coupling(name = 'GC_79',
                 value = '(ee*complex(0,1)*complexconjugate(CKM3x3))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

