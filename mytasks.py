#
# User defined tasks setup.
# Generated from buildmytask.
#

import sys
from casa_stack_manip import stack_frame_find

if sys.path[1] != '/data/SHARED/WORKAREA/ARC_TOOLS/CASA-PolTools/trunk':
  sys.path.insert(1, '/data/SHARED/WORKAREA/ARC_TOOLS/CASA-PolTools/trunk')
from odict import odict
if not globals().has_key('mytasks') :
  mytasks = odict()

mytasks['polsimulate'] = 'Version 1.3.2 - Basic simulator of ALMA/J-VLA (and VLBI) full-polarization observations. The output should be imaged with CLEAN (with stokes=IQUV) and the polarization vectors should be computed with immath (with options poli and pola). See the ALMA Polarization CASA Guide for more information.\n\n'
mytasks['polsolve'] = 'Version 1.0.1b - Leakage solver for circular polarizers and extended polarization calibrators.\n\n'

if not globals().has_key('task_location') :
  task_location = odict()

task_location['polsimulate'] = '/data/SHARED/WORKAREA/ARC_TOOLS/CASA-PolTools/trunk'
task_location['polsolve'] = '/data/SHARED/WORKAREA/ARC_TOOLS/CASA-PolTools/trunk'
myglobals = stack_frame_find( )
tasksum = myglobals['tasksum'] 
for key in mytasks.keys() :
  tasksum[key] = mytasks[key]

from polsimulate_cli import  polsimulate_cli as polsimulate
from polsolve_cli import  polsolve_cli as polsolve
