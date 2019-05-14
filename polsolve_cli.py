#
# This file was generated using xslt from its XML file
#
# Copyright 2014, Associated Universities Inc., Washington DC
#
import sys
import os
#from casac import *
import casac
import string
import time
import inspect
import gc
import numpy
from casa_stack_manip import stack_frame_find
from odict import odict
from types import *
from task_polsolve import polsolve
class polsolve_cli_:
    __name__ = "polsolve"
    rkey = None
    i_am_a_casapy_task = None
    # The existence of the i_am_a_casapy_task attribute allows help()
    # (and other) to treat casapy tasks as a special case.

    def __init__(self) :
       self.__bases__ = (polsolve_cli_,)
       self.__doc__ = self.__call__.__doc__

       self.parameters={'vis':None, 'spw':None, 'field':None, 'mounts':None, 'feed_rotation':None, 'DR':None, 'DL':None, 'DRSolve':None, 'DLSolve':None, 'CLEAN_models':None, 'Pfrac':None, 'EVPA':None, 'PolSolve':None, 'parang_corrected':None, 'target_field':None, 'plot_parang':None, 'min_elev_plot':None, 'wgt_power':None, }


    def result(self, key=None):
	    #### and add any that have completed...
	    return None


    def __call__(self, vis=None, spw=None, field=None, mounts=None, feed_rotation=None, DR=None, DL=None, DRSolve=None, DLSolve=None, CLEAN_models=None, Pfrac=None, EVPA=None, PolSolve=None, parang_corrected=None, target_field=None, plot_parang=None, min_elev_plot=None, wgt_power=None, ):

        """Version 1.0.1b - Leakage solver for circular polarizers and extended polarization calibrators.\n\n

	Detailed Description:
Version 1.0.1b - Leakage solver for circular polarizers and extended polarization calibrators.\n\n
	Arguments :
		vis:	Name of input measurement set. The data\nshould already be calibrated \n(in bandpass and gains).
		   Default Value: input.ms

		spw:	Spectral window to fit for.
		   Default Value: 0

		field:	Field name (or id) to use as calibrator.
		   Default Value: 0

		mounts:	List of the antenna mounts (must be given\n in the same order as the ANTENNA table \nof the ms). Use this in case the mounts \nwere not properly imported into the ms.\nA mount type is specified with two characters.\nSupported mounts: \nalt-az (\'AZ\'), \nequatorial (\'EQ\'), \nX-Y (\'XY\'), \nNasmyth Right (\'NR\') and \nNasmyth Left (\'NL\'). \nDefault means all antennas are alt-az.
		   Default Value: []

		feed_rotation:	Rotation of the feed of each antenna with respect to the local horizontal-vertical frame. One value per antenna. Empty list assumes a null feed angle for all antennas.
		   Default Value: []

		DR:	List of complex numbers (length equal\n to the number of antennas). If not empty, \n these are the a-priori values of \n the DR leakage terms to use \n in the fit.
		   Default Value: []

		DL:	List of complex numbers (length equal\n to the number of antennas). If not empty, \n these are the a-priori values of \n the DL leakage terms to use \n in the fit.
		   Default Value: []

		DRSolve:	List of booleans (length equal\nto the number of antennas). If not empty,\n it will tell which DR terms are fitted. \n The DR[i] term will be fitted if \n DRSolve[i] is True. Otherwise, DR[i] \n will be fixed to its a-priori value. \n Default (i.e., empty list) means \n to fit all DRs. 
		   Default Value: []

		DLSolve:	Just as DRSolve, but for the DL terms.
		   Default Value: []

		CLEAN_models:	List of CLEAN model files (CCs given \n in PRTAB format). Each file will \n correspond to a source component with \n the same polarization state. If one number \n is given (instead of a list of filenames),\n a centered point source (with that \n flux density) will be used.
		   Default Value: [1.0]

		Pfrac:	List of fractional polarizations (one number \n per source component). Pfrac \n values must fall between 0 and 1.
		   Default Value: [0.0]

		EVPA:	List of EVPAs in degrees (one number \n per source component). Angles \n are measured from North to East.
		   Default Value: [0.0]

		PolSolve:	List of booleans (one per source component) \n that tell which source components \n are to be fitted in polarization. \n If PolSolve[i] is True, the fractional \n polarization and EVPA of the ith source \n component will be fitted, together with \n the antenna Dterms. If False, all \n Stokes parameters of the ith component will \n be fixed in the fit. \n Empty list means to fit the polarization \n of all the source components.
		   Default Value: []

		parang_corrected:	If True, the data are assumed to be \n already corrected for parallactic angle. \n This is usually the case, unless \n you are working with data generated with \n polsimulate with no parang correction.
		   Default Value: True

		target_field:	List of sources to which apply the Dterm (and parangle) correction. It must follow the CASA syntax if a range of field ids is given. Empy list means NOT to apply the Dterms (i.e., just save them in a calibration table). If you want to apply the calibration, DO NOT FORGET TO *ALWAYS* RUN CLEARCAL BEFORE POLSOLVE!!
		   Default Value: 

		plot_parang:	If True, plot the time evolution of the antenna feed angles (i.e., parallactic angle plus correction from the antenna mounts).
		   Default Value: False

		min_elev_plot:	 In degrees. If plot_parang is True, points with elevations lower than this limit will be plotted in red. THIS DOES NOT FLAG THE DATA. If you want to flag them, run the flagdata task.
		   Default Value: 10.0

		wgt_power:	 Power for the visibility weights. Unity means to leave the weights untouched (i.e., equivalent to natural weighting, but for the fit). Zero means equal weights for all visibilities (i.e., equivalent to uniform weighting for the fit).
		   Default Value: 1.0

	Returns: bool

	Example :


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
	if not hasattr(self, "__globals__") or self.__globals__ == None :
           self.__globals__=stack_frame_find( )
	#casac = self.__globals__['casac']
	casalog = self.__globals__['casalog']
	casa = self.__globals__['casa']
	#casalog = casac.casac.logsink()
        self.__globals__['__last_task'] = 'polsolve'
        self.__globals__['taskname'] = 'polsolve'
        ###
        self.__globals__['update_params'](func=self.__globals__['taskname'],printtext=False,ipython_globals=self.__globals__)
        ###
        ###
        #Handle globals or user over-ride of arguments
        #
        if type(self.__call__.func_defaults) is NoneType:
            function_signature_defaults={}
	else:
	    function_signature_defaults=dict(zip(self.__call__.func_code.co_varnames[1:],self.__call__.func_defaults))
	useLocalDefaults = False

        for item in function_signature_defaults.iteritems():
                key,val = item
                keyVal = eval(key)
                if (keyVal == None):
                        #user hasn't set it - use global/default
                        pass
                else:
                        #user has set it - use over-ride
			if (key != 'self') :
			   useLocalDefaults = True

	myparams = {}
	if useLocalDefaults :
	   for item in function_signature_defaults.iteritems():
	       key,val = item
	       keyVal = eval(key)
	       exec('myparams[key] = keyVal')
	       self.parameters[key] = keyVal
	       if (keyVal == None):
	           exec('myparams[key] = '+ key + ' = self.itsdefault(key)')
		   keyVal = eval(key)
		   if(type(keyVal) == dict) :
                      if len(keyVal) > 0 :
		         exec('myparams[key] = ' + key + ' = keyVal[len(keyVal)-1][\'value\']')
		      else :
		         exec('myparams[key] = ' + key + ' = {}')

        else :
            print ''

            myparams['vis'] = vis = self.parameters['vis']
            myparams['spw'] = spw = self.parameters['spw']
            myparams['field'] = field = self.parameters['field']
            myparams['mounts'] = mounts = self.parameters['mounts']
            myparams['feed_rotation'] = feed_rotation = self.parameters['feed_rotation']
            myparams['DR'] = DR = self.parameters['DR']
            myparams['DL'] = DL = self.parameters['DL']
            myparams['DRSolve'] = DRSolve = self.parameters['DRSolve']
            myparams['DLSolve'] = DLSolve = self.parameters['DLSolve']
            myparams['CLEAN_models'] = CLEAN_models = self.parameters['CLEAN_models']
            myparams['Pfrac'] = Pfrac = self.parameters['Pfrac']
            myparams['EVPA'] = EVPA = self.parameters['EVPA']
            myparams['PolSolve'] = PolSolve = self.parameters['PolSolve']
            myparams['parang_corrected'] = parang_corrected = self.parameters['parang_corrected']
            myparams['target_field'] = target_field = self.parameters['target_field']
            myparams['plot_parang'] = plot_parang = self.parameters['plot_parang']
            myparams['min_elev_plot'] = min_elev_plot = self.parameters['min_elev_plot']
            myparams['wgt_power'] = wgt_power = self.parameters['wgt_power']


	result = None

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
	pathname="file:///data/SHARED/WORKAREA/ARC_TOOLS/CASA-PolTools/trunk/"
	trec = casac.casac.utils().torecord(pathname+'polsolve.xml')

        casalog.origin('polsolve')
	try :
          #if not trec.has_key('polsolve') or not casac.casac.utils().verify(mytmp, trec['polsolve']) :
	    #return False

          casac.casac.utils().verify(mytmp, trec['polsolve'], True)
          scriptstr=['']
          saveinputs = self.__globals__['saveinputs']
          if type(self.__call__.func_defaults) is NoneType:
              saveinputs=''
          else:
              saveinputs('polsolve', 'polsolve.last', myparams, self.__globals__,scriptstr=scriptstr)
          tname = 'polsolve'
          spaces = ' '*(18-len(tname))
          casalog.post('\n##########################################'+
                       '\n##### Begin Task: ' + tname + spaces + ' #####')
          if (casa['state']['telemetry-enabled']):
              casalog.poststat('Begin Task: ' + tname)
          if type(self.__call__.func_defaults) is NoneType:
              casalog.post(scriptstr[0]+'\n', 'INFO')
          else :
              casalog.post(scriptstr[1][1:]+'\n', 'INFO')
          result = polsolve(vis, spw, field, mounts, feed_rotation, DR, DL, DRSolve, DLSolve, CLEAN_models, Pfrac, EVPA, PolSolve, parang_corrected, target_field, plot_parang, min_elev_plot, wgt_power)
          if (casa['state']['telemetry-enabled']):
              casalog.poststat('End Task: ' + tname)
          casalog.post('##### End Task: ' + tname + '  ' + spaces + ' #####'+
                       '\n##########################################')

	except Exception, instance:
          if(self.__globals__.has_key('__rethrow_casa_exceptions') and self.__globals__['__rethrow_casa_exceptions']) :
             raise
          else :
             #print '**** Error **** ',instance
	     tname = 'polsolve'
             casalog.post('An error occurred running task '+tname+'.', 'ERROR')
             pass
	casalog.origin('')

        gc.collect()
        return result
#
#
#
#    def paramgui(self, useGlobals=True, ipython_globals=None):
#        """
#        Opens a parameter GUI for this task.  If useGlobals is true, then any relevant global parameter settings are used.
#        """
#        import paramgui
#	if not hasattr(self, "__globals__") or self.__globals__ == None :
#           self.__globals__=stack_frame_find( )
#
#        if useGlobals:
#	    if ipython_globals == None:
#                myf=self.__globals__
#            else:
#                myf=ipython_globals
#
#            paramgui.setGlobals(myf)
#        else:
#            paramgui.setGlobals({})
#
#        paramgui.runTask('polsolve', myf['_ip'])
#        paramgui.setGlobals({})
#
#
#
#
    def defaults(self, param=None, ipython_globals=None, paramvalue=None, subparam=None):
	if not hasattr(self, "__globals__") or self.__globals__ == None :
           self.__globals__=stack_frame_find( )
        if ipython_globals == None:
            myf=self.__globals__
        else:
            myf=ipython_globals

        a = odict()
        a['vis']  = 'input.ms'
        a['spw']  = 0
        a['field']  = '0'
        a['mounts']  = []
        a['feed_rotation']  = []
        a['DR']  = []
        a['DL']  = []
        a['DRSolve']  = []
        a['DLSolve']  = []
        a['CLEAN_models']  = [1.0]
        a['Pfrac']  = [0.0]
        a['EVPA']  = [0.0]
        a['PolSolve']  = []
        a['parang_corrected']  = True
        a['target_field']  = ''
        a['plot_parang']  = False
        a['min_elev_plot']  = 10.0
        a['wgt_power']  = 1.0


### This function sets the default values but also will return the list of
### parameters or the default value of a given parameter
        if(param == None):
                myf['__set_default_parameters'](a)
        elif(param == 'paramkeys'):
                return a.keys()
        else:
            if(paramvalue==None and subparam==None):
               if(a.has_key(param)):
                  return a[param]
               else:
                  return self.itsdefault(param)
            else:
               retval=a[param]
               if(type(a[param])==dict):
                  for k in range(len(a[param])):
                     valornotval='value'
                     if(a[param][k].has_key('notvalue')):
                        valornotval='notvalue'
                     if((a[param][k][valornotval])==paramvalue):
                        retval=a[param][k].copy()
                        retval.pop(valornotval)
                        if(subparam != None):
                           if(retval.has_key(subparam)):
                              retval=retval[subparam]
                           else:
                              retval=self.itsdefault(subparam)
		     else:
                        retval=self.itsdefault(subparam)
               return retval


#
#
    def check_params(self, param=None, value=None, ipython_globals=None):
      if ipython_globals == None:
          myf=self.__globals__
      else:
          myf=ipython_globals
#      print 'param:', param, 'value:', value
      try :
         if str(type(value)) != "<type 'instance'>" :
            value0 = value
            value = myf['cu'].expandparam(param, value)
            matchtype = False
            if(type(value) == numpy.ndarray):
               if(type(value) == type(value0)):
                  myf[param] = value.tolist()
               else:
                  #print 'value:', value, 'value0:', value0
                  #print 'type(value):', type(value), 'type(value0):', type(value0)
                  myf[param] = value0
                  if type(value0) != list :
                     matchtype = True
            else :
               myf[param] = value
            value = myf['cu'].verifyparam({param:value})
            if matchtype:
               value = False
      except Exception, instance:
         #ignore the exception and just return it unchecked
         myf[param] = value
      return value
#
#
    def description(self, key='polsolve', subkey=None):
        desc={'polsolve': 'Version 1.0.1b - Leakage solver for circular polarizers and extended polarization calibrators.\n\n',
               'vis': 'Name of input measurement set. The data\nshould already be calibrated \n(in bandpass and gains).',
               'spw': 'Spectral window to fit for.',
               'field': 'Field name (or id) to use as calibrator.',
               'mounts': 'List of the antenna mounts (must be given\n in the same order as the ANTENNA table \nof the ms). Use this in case the mounts \nwere not properly imported into the ms.\nA mount type is specified with two characters.\nSupported mounts: \nalt-az (\'AZ\'), \nequatorial (\'EQ\'), \nX-Y (\'XY\'), \nNasmyth Right (\'NR\') and \nNasmyth Left (\'NL\'). \nDefault means all antennas are alt-az.',
               'feed_rotation': 'Rotation of the feed of each antenna with respect to the local horizontal-vertical frame. One value per antenna. Empty list assumes a null feed angle for all antennas.',
               'DR': 'List of complex numbers (length equal\n to the number of antennas). If not empty, \n these are the a-priori values of \n the DR leakage terms to use \n in the fit.',
               'DL': 'List of complex numbers (length equal\n to the number of antennas). If not empty, \n these are the a-priori values of \n the DL leakage terms to use \n in the fit.',
               'DRSolve': 'List of booleans (length equal\nto the number of antennas). If not empty,\n it will tell which DR terms are fitted. \n The DR[i] term will be fitted if \n DRSolve[i] is True. Otherwise, DR[i] \n will be fixed to its a-priori value. \n Default (i.e., empty list) means \n to fit all DRs. ',
               'DLSolve': 'Just as DRSolve, but for the DL terms.',
               'CLEAN_models': 'List of CLEAN model files (CCs given \n in PRTAB format). Each file will \n correspond to a source component with \n the same polarization state. If one number \n is given (instead of a list of filenames),\n a centered point source (with that \n flux density) will be used.',
               'Pfrac': 'List of fractional polarizations (one number \n per source component). Pfrac \n values must fall between 0 and 1.',
               'EVPA': 'List of EVPAs in degrees (one number \n per source component). Angles \n are measured from North to East.',
               'PolSolve': 'List of booleans (one per source component) \n that tell which source components \n are to be fitted in polarization. \n If PolSolve[i] is True, the fractional \n polarization and EVPA of the ith source \n component will be fitted, together with \n the antenna Dterms. If False, all \n Stokes parameters of the ith component will \n be fixed in the fit. \n Empty list means to fit the polarization \n of all the source components.',
               'parang_corrected': 'If True, the data are assumed to be \n already corrected for parallactic angle. \n This is usually the case, unless \n you are working with data generated with \n polsimulate with no parang correction.',
               'target_field': 'List of sources to which apply the Dterm (and parangle) correction. It must follow the CASA syntax if a range of field ids is given. Empy list means NOT to apply the Dterms (i.e., just save them in a calibration table). If you want to apply the calibration, DO NOT FORGET TO *ALWAYS* RUN CLEARCAL BEFORE POLSOLVE!!',
               'plot_parang': 'If True, plot the time evolution of the antenna feed angles (i.e., parallactic angle plus correction from the antenna mounts).',
               'min_elev_plot': ' In degrees. If plot_parang is True, points with elevations lower than this limit will be plotted in red. THIS DOES NOT FLAG THE DATA. If you want to flag them, run the flagdata task.',
               'wgt_power': ' Power for the visibility weights. Unity means to leave the weights untouched (i.e., equivalent to natural weighting, but for the fit). Zero means equal weights for all visibilities (i.e., equivalent to uniform weighting for the fit).',

              }

        if(desc.has_key(key)) :
           return desc[key]

    def itsdefault(self, paramname) :
        a = {}
        a['vis']  = 'input.ms'
        a['spw']  = 0
        a['field']  = '0'
        a['mounts']  = []
        a['feed_rotation']  = []
        a['DR']  = []
        a['DL']  = []
        a['DRSolve']  = []
        a['DLSolve']  = []
        a['CLEAN_models']  = [1.0]
        a['Pfrac']  = [0.0]
        a['EVPA']  = [0.0]
        a['PolSolve']  = []
        a['parang_corrected']  = True
        a['target_field']  = ''
        a['plot_parang']  = False
        a['min_elev_plot']  = 10.0
        a['wgt_power']  = 1.0

        #a = sys._getframe(len(inspect.stack())-1).f_globals

        if a.has_key(paramname) :
	      return a[paramname]
polsolve_cli = polsolve_cli_()
