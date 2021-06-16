import os
import sys

#import re
import shutil

mass=float(sys.argv[1])
qcut=mass/4.0




try:
    IF=open('lhe_MLM.cfg','r')
except:
    raise SystemExit


try:
    os.remove('MLM.cfg')
except:
    pass

try:
    OF=open('MLM.cfg','w')
except:
    raise SystemExit

foundline=False

for line in IF:
    if 'JetMatching:qCut' in line:
        OF.write('JetMatching:qCut = %.1f  ! minimum kt jet measure between partons\n' %qcut)
        foundline=True
    else:
        OF.write(line)
if not foundline:
    OF.write('JetMatching:qCut = %.1f  ! minimum kt jet measure between partons\n' %qcut)
IF.close()
OF.close()





try:
    IF=open('lhe_CKKWL.cfg','r')
except:
    raise SystemExit


try:
    os.remove('CKKWL.cfg')
except:
    pass

try:
    OF=open('CKKWL.cfg','w')
except:
    raise SystemExit

foundline=False

for line in IF:
    if 'Merging:TMS' in line:
        OF.write('Merging:TMS = %.1f  !merging scale value\n' %qcut)
        foundline=True
    else:
        OF.write(line)
if not foundline:
    OF.write('Merging:TMS = = %.1f  !merging scale value \n' %qcut)
IF.close()
OF.close()


