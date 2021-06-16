#include "Analyses/HSCP_ATLAS.h"
#include "Analyses/DT_CMS.h"


// Here we use macros to set up the list of all analyses
// and functions to check whether to run them. Sadly we can't
// also generate the include lines above. 
// However, deleting an entry from here will completely disable the analysis;
// we also have settings at runtime to do that on a per-run basis.

#define MAP_ANALYSES(F) \
F(HSCP_ATLAS)\
F(DT_CMS)
