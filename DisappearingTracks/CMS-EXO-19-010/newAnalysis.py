import os
import sys

import regenerate

if len(sys.argv) < 2:
    print("Need to give an analysis name!")
    raise SystemExit

ANAME=sys.argv[1]

ADIR='Analyses/'

HNAME=ADIR+ANAME+'.h'
CCNAME=ADIR+ANAME+'.cc'

HHNAME=ANAME+'_h'


if os.path.exists(HNAME) or os.path.exists(CCNAME):
    print('Analysis already exists! If you want to create a new blank version, just delete the .h and .cc files')
    raise SystemExit

############### Write HEADER!

try:
    OF=open(HNAME,'w')
except:
    print('Could not create header file')
    raise SystemExit

#OF.write('#pragma once\n\n')


OF.write('#ifndef '+HHNAME+'\n')
OF.write('#define '+HHNAME+'\n')

OF.write('#include "include/BaseAnalysis.h"\n')
OF.write('#include "include/HEPData.h"\n\n')

OF.write('using namespace std;\n\n')

OF.write('class '+ANAME+' : public BaseAnalysis {\n')
OF.write('    public:\n\n')
OF.write('       ~'+ANAME+'()\n')
OF.write('       void init();\n\n')
OF.write('       void Execute(std::mt19937 &engine);\n\n')
OF.write('       void Finalise() { return; }\n\n')
OF.write('       '+ANAME+'() { this->setup(); this->analysisname="'+ANAME+'"; this->init(); };\n\n')
OF.write('};\n\n')

OF.write('#endif\n')
OF.close()




############### Write CPP file!

try:
    OF=open(CCNAME,'w')
except:
    print('Could not create source file')
    raise SystemExit

#OF.write('#pragma once\n\n')

OF.write('#include "'+ANAME+'.h"\n\n')

OF.write('void '+ANAME+'::init() {\n\n')
OF.write('cout << "---------------------------------------------------" << endl;\n')
OF.write('cout << "-- Analysis '+ANAME+' --" << endl;\n')
OF.write('cout << "---------------------------------------------------" << endl;\n')
OF.write('}\n\n')

OF.write('void '+ANAME+'::~'+ANAME+'() {\n\n')
OF.write('}\n\n')

OF.write('void '+ANAME+'::Execute(std::mt19937 &engine) {\n\n')
OF.write('}\n\n')


OF.close()


###### Now to regenerate the analysis list

regenerate.regenerate_analysislist()






