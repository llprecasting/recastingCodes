def regenerate_analysislist():
    import os
    import sys
    import fnmatch
    try:
        OF=open('analysislist.h','w')
    except:
        print('Could not recreate the analysis list!')
        raise SystemExit


    ANALYSES=[]
    AHEADS=[]
    for f in os.listdir('Analyses'):
        if not fnmatch.fnmatch(f,'*.h'):
            continue

        try:
            ANAME=os.path.splitext(f)[0]
        except:
            continue
        ANALYSES.append(ANAME)
        AHEADS.append(f)

    for f in AHEADS:
        OF.write('#include "Analyses/'+f+'"\n')

    OF.write('''\n\n// Here we use macros to set up the list of all analyses
// and functions to check whether to run them. Sadly we can't
// also generate the include lines above. 
// However, deleting an entry from here will completely disable the analysis;
// we also have settings at runtime to do that on a per-run basis.\n\n''')

    OF.write('#define MAP_ANALYSES(F) ')
    for f in ANALYSES:
        OF.write('\\'+'\nF('+f+')')

    OF.write('\n')
    OF.close()

if __name__ == "__main__":
    regenerate_analysislist()
