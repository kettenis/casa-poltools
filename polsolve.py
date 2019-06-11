#
# This file was generated using xslt from its XML file
#
# Copyright 2009, Associated Universities Inc., Washington DC
#
import sys
import os
from  casac import *
import string
from taskinit import casalog
from taskinit import xmlpath
#from taskmanager import tm
import task_polsolve
def polsolve(vis='input.ms', spw=0, field='0', mounts=[], feed_rotation=[], DR=[], DL=[], DRSolve=[], DLSolve=[], CLEAN_models=[1.0], Pfrac=[0.0], EVPA=[0.0], PolSolve=[], parang_corrected=True, target_field='', plot_parang=False, min_elev_plot=10.0, wgt_power=1.0):

        """Version 1.0.1b - Leakage solver for circular polarizers and extended polarization calibrators.\n\n

  TODO: Combine IFs (i.e., spws) and parameterize the frequency dependence of the Dterms!	  

  POLSOLVE EXAMPLES:

  Let's suppose that we have a CLEAN image in AIPS with two CC tables from IMAGR (one per source component).
  We produce two ascii files (e.g., CC1.dat and CC2.dat) with the output of PRTAB:

  PRTAB(BPRINT=1, EPRINT=0, NDIG = 8, OUTPR = \'PWD:CC1.dat\', DOCRT = -1, INVER = 1, INEXT = \'CC\')
  PRTAB(BPRINT=1, EPRINT=0, NDIG = 8, OUTPR = \'PWD:CC2.dat\', DOCRT = -1, INVER = 2, INEXT = \'CC\')

  (if polsolve fails using these files, try to remove all unnecessary lines from them, and only keep the 
  model information).

  Then, if we want to solve for all Dterms and the Stokes parameters of these two components, the keywords 
  to use would be:

  field = \'Name of Calibrator\'
  CLEAN_models = [\'CC1.dat\' , \'CC2.dat\']
  Pfrac = [0., 0.]
  EVPA = [0., 0.]

  -If there are 8 antennas, and one of them (e.g., the fifth one) has a Nasmyth-Left mount:

  mounts = [\'AZ\' for i in range(8)] ; mounts[4] = \'NL\'

  
  -If we want to fix (i.e., to NOT fit) the Dterms for L of the first antenna:

  DLSolve = [True for i in range(8)] ; DLSolve[0] = False

  
  -If we want to fix the fractional polarization of the first component to 10% (with EVPA of 20 deg.):

  Pfrac    = [0.10,   0.0]
  EVPA     = [  20.,  0.0]
  PolSolve = [False, True]


  After running, the task will create a Dterms table, with the name equal to that of the input
  measurement set plus the suffix \'.spwI.Dterms\' (where \'I\' is the number of the spectral 
  window used in the fit).




        """

#
#    The following is work around to avoid a bug with current python translation
#
        mytmp = {}

        mytmp['vis'] = vis
        mytmp['spw'] = spw
        mytmp['field'] = field
        mytmp['mounts'] = mounts
        mytmp['feed_rotation'] = feed_rotation
        mytmp['DR'] = DR
        mytmp['DL'] = DL
        mytmp['DRSolve'] = DRSolve
        mytmp['DLSolve'] = DLSolve
        mytmp['CLEAN_models'] = CLEAN_models
        mytmp['Pfrac'] = Pfrac
        mytmp['EVPA'] = EVPA
        mytmp['PolSolve'] = PolSolve
        mytmp['parang_corrected'] = parang_corrected
        mytmp['target_field'] = target_field
        mytmp['plot_parang'] = plot_parang
        mytmp['min_elev_plot'] = min_elev_plot
        mytmp['wgt_power'] = wgt_power
	pathname="file:///home/marti/WORKAREA/ARC_TOOLS/CASA-PolTools/casa-poltools/"
	trec = casac.utils().torecord(pathname+'polsolve.xml')

        casalog.origin('polsolve')
        if trec.has_key('polsolve') and casac.utils().verify(mytmp, trec['polsolve']) :
	    result = task_polsolve.polsolve(vis, spw, field, mounts, feed_rotation, DR, DL, DRSolve, DLSolve, CLEAN_models, Pfrac, EVPA, PolSolve, parang_corrected, target_field, plot_parang, min_elev_plot, wgt_power)

	else :
	  result = False
        return result
