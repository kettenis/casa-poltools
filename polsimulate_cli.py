#
# This file was generated using xslt from its XML file
#
# Copyright 2014, Associated Universities Inc., Washington DC
#
import sys
import os
import datetime
#from casac import *
import casac
import string
import time
import inspect
import numpy
from casa_stack_manip import stack_frame_find
from odict import odict
from types import *
from task_polsimulate import polsimulate
class polsimulate_cli_:
    __name__ = "polsimulate"
    rkey = None
    i_am_a_casapy_task = None
    # The existence of the i_am_a_casapy_task attribute allows help()
    # (and other) to treat casapy tasks as a special case.

    def __init__(self) :
       self.__bases__ = (polsimulate_cli_,)
       self.__doc__ = self.__call__.__doc__

       self.parameters={'vis':None, 'reuse':None, 'array_configuration':None, 'elevation_cutoff':None, 'feed':None, 'mounts':None, 'ConstDt0':None, 'ConstDt1':None, 'LO':None, 'BBs':None, 'spw_width':None, 'nchan':None, 'model_image':None, 'I':None, 'Q':None, 'U':None, 'V':None, 'RM':None, 'spec_index':None, 'RAoffset':None, 'Decoffset':None, 'spectrum_file':None, 'phase_center':None, 'incell':None, 'inbright':None, 'inwidth':None, 'innu0':None, 'H0':None, 'onsource_time':None, 'observe_time':None, 'visib_time':None, 'nscan':None, 'apply_parang':None, 'export_uvf':None, 'corrupt':None, 'seed':None, 'Dt_amp':None, 'Dt_noise':None, 'tau0':None, 't_sky':None, 't_ground':None, 't_receiver':None, }


    def result(self, key=None):
	    #### and add any that have completed...
	    return None


    def __call__(self, vis=None, reuse=None, array_configuration=None, elevation_cutoff=None, feed=None, mounts=None, ConstDt0=None, ConstDt1=None, LO=None, BBs=None, spw_width=None, nchan=None, model_image=None, I=None, Q=None, U=None, V=None, RM=None, spec_index=None, RAoffset=None, Decoffset=None, spectrum_file=None, phase_center=None, incell=None, inbright=None, inwidth=None, innu0=None, H0=None, onsource_time=None, observe_time=None, visib_time=None, nscan=None, apply_parang=None, export_uvf=None, corrupt=None, seed=None, Dt_amp=None, Dt_noise=None, tau0=None, t_sky=None, t_ground=None, t_receiver=None, ):

        """Version 1.3.2 - Basic simulator of ALMA/J-VLA (and VLBI) full-polarization observations. The output should be imaged with CLEAN (with stokes=IQUV) and the polarization vectors should be computed with immath (with options poli and pola). See the ALMA Polarization CASA Guide for more information.\n\n

	Detailed Description:
Version 1.3.2 - Basic simulator of ALMA/J-VLA (and VLBI) full-polarization observations. The output should be imaged with CLEAN (with stokes=IQUV) and the polarization vectors should be computed with immath (with options poli and pola). See the ALMA Polarization CASA Guide for more information.\n\n
	Arguments :
		vis:	Name of output measurement set.
		   Default Value: polsimulate_output.ms

		reuse:	If True, and the measurement set exists, will reuse it (so the antenna names and source coordinates will be reused as well).
		   Default Value: False

		array_configuration:	Array configuration file, where the antenna coordinates and diameters are set (see the files in data/alma/simmos of your CASA path). Default is an ALMA configuration. The name of the array is taken as the part of the name of this file before the first dot (i.e., it is \'alma\' for the default filename).
		   Default Value: alma.out04.cfg

		elevation_cutoff:	Minimum allowed elevation for the antennas (degrees).
		   Default Value: 5.0

		feed:	Polarization basis for the measurement set. Can be linear (e.g., for ALMA) or circular (e.g., for VLA). Default is linear.
		   Default Value: linear

		mounts:	For VLBI observations, this is a list of the antenna mounts (given in the same order as in the array configuration file). A mount type is specified with two characters. Supported mounts: alt-az (\'AZ\'), equatorial (\'EQ\'), X-Y (\'XY\'), Nasmyth Right (\'NR\') and Nasmyth Left (\'NL\'). Default means all antennas are alt-az.
		   Default Value: []

		ConstDt0:	List of complex numbers (length equal to the number of antennas). If not empty, the first polarizer (i.e., either R or X, depending on the value of \'feed\') will be contamined with a leakage given by these Dterms.
		   Default Value: []

		ConstDt1:	List of complex numbers (length equal to the number of antennas). If not empty, the second polarizer will be contamined with a leakage given by these Dterms.
		   Default Value: []

		LO:	Frequency of the first LO in Hz (this will define the ALMA band of observation). Default is 100 GHz (i.e., ALMA Band 3)
		   Default Value: 100.e9

		BBs:	List with the baseband frequency offsets (in Hz). There will be one spectral window per baseband. 
		   Default Value: [-7.e9, -5.e9, 5.e9, 7.e9]

		spw_width:	Width of the spectral windows in Hz (the spws will be centered at each BB).
		   Default Value: 2.e9

		nchan:	Number of channels in each spw (all spws will have the same number of channels).
		   Default Value: 8

		model_image:	List of four images (Stokes I, Q, U, and V, respectively) to be used as observed extended sources. Image cubes are allowed. Default is to NOT simulate extended sources. BEWARE OF THE CURRENT ALMA LIMITATION for extended sources (should fall within the inner 1/3 of the primary beam FWHM).
		   Default Value: []

		I:	List of Stokes I (in Jy) for a set of point sources to simulate. These are added to the \'model_image\' (if it is provided). Default is to add NO source. The flux density is referred to the LO frequency. Example: [1.0] for a 1 Jy source.
		   Default Value: []

		Q:	List of Stokes Q (in Jy) for the sources defined above. Default is no source. Example: [0.0] for no Q signal. These values are referred to the LO frequency.
		   Default Value: []

		U:	>List of Stokes U (in Jy) for the sources defined above. Default is no source. Example: [0.0] for no U signal. These values are referred to the LO frequency.
		   Default Value: []

		V:	>List of Stokes V (in Jy) for the sources defined above. Default is no source. Example: [0.0] for no V signal. These values are referred to the LO frequency.
		   Default Value: []

		RM:	List of Rotation Measures (RM, in rad/m**2.) for the sources defined above. Default is no source. Example: [0.0] for no RM.
		   Default Value: []

		spec_index:	List of spectral indices for the sources defined above. Default is no source. Example: [0.0] for a flat spectrum.
		   Default Value: []

		RAoffset:	List of right-ascension offsets (in arcsec) for the sources defined above. The first source is assumed to be at the phase center, so all sources will be shifted RAoffset[0] arcsec (so that the first source in the list is at the phase center).
		   Default Value: []

		Decoffset:	List of declination offsets (in arcsec) for the sources defined above. The first source is assumed to be at the phase center, so all sources will be shifted Decoffset[0] arcsec.
		   Default Value: []

		spectrum_file:	File with user-defined spectra of I, Q, U, and V. See help for details about the file format. This source WILL BE ADDED TO THE PHASE CENTER, together with the source defined in the model_image model and all those given in the I, Q, U, V, lists.
		   Default Value: 

		phase_center:	Coordinates of the observed source (will override the coordinates defined in model_image, if an image is being used). This keyword MUST BE defined.
		   Default Value: J2000 00h00m00.00 -00d00m00.00

		incell:	Pixel size of the model_image. If not empty, will override the original value stored in the image. Example: \'0.1arcsec\'. All the Stokes images (I, Q, U, and V) will be set the same way. USE WITH CARE.
		   Default Value: 

		inbright:	Peak intensity of the I Stokes model_image (will override the original value stored in the image). Default is to use the original brightness unit. All the Stokes images (I, Q, U, and V) will be set the same way. Default (i.e., 0.0) means to NOT scale the image birghtness. USE WITH CARE.
		   Default Value: 0.0

		inwidth:	Width of the frequency channels in the model_image. If not empty, will override the original value stored in the image. Example: \'10MHz\'. All the Stokes images (I, Q, U, and V) will be set the same way.
		   Default Value: 

		innu0:	Frequency of the first image channel (e.g., \'1GHz\'). Default (empty) means to use the value in the image.
		   Default Value: 

		H0:	If a float is given, it is the Hour Angle at the start of the observations, as observed from the array center (in hr). If a string is given, it is the exact start of the observations (UT time) in the format \'2017/01/01/00:00:00\'.
		   Default Value: -1.5

		onsource_time:	Total effective integration time on the source (in hr). Default is 1 hr.
		   Default Value: 1.0

		observe_time:	Total observing time (i.e., including overheads) in hr. Default is 3h, so there will be an observing efficiency of 0.33 if \'onsource_time\' is set to one hour. 
		   Default Value: 3.0

		visib_time:	Integration time per visibility. This is a string, where \'s\' stands for \'seconds\'.
		   Default Value: 6s

		nscan:	Number of scans. Can be provided as a \'listobs\' file \n(in that case, observe_time, onsource_time \n and H0 are not used, but taken from the listobs).\n If a list in \'listobs\' format is not provided \nthen all scans will be set to equal length. \n If just an integer is given, the scans will be \n homogeneously distributed across the \n total observing time. If a list is given, \n the values will be taken as the starting times \n of the scans, relative to the \n duration of the experiment. For instance, if \n \'observe_time = 6.0\' then \'nscan = [0., 0.5, 0.75]\' \n will make three scans, one at the start of the \n observations, the second one 3 hours later \n and the third one 4.5 hours after \n the start of the observations.
		   Default Value: 50

		apply_parang:	If True, applies the parallactic-angle correction. \n If False, the data polarization will be given \n in the antenna frame (i.e., just as true raw data).
		   Default Value: False

		export_uvf:	If True, exports the measurement into uvfits format (for its use in e.g., AIPS/Difmap).
		   Default Value: True

		corrupt:	Whether to add random noise to the visibilities.
		   Default Value: True

		seed:	Seed of the random number generator in the sm tool.
		   Default Value: 42

		Dt_amp:	Will add random Dterms (antenna-wise). Dt_amp is the typical absolute value of the Dterms (real and imag). The actual values will be computed from a random Gaussian distribution.
		   Default Value: 0.0

		Dt_noise:	Will add random channel-dependent contribution to the Dterms. Dt_noise is the typical residual channel noise in the Dterms (real and imag). The actual values for each frequency channel will be those of Dt_amp PLUS a random Gaussian distribution of width Dt_noise. Default is 0.001, similar to the spectral spread of Dterms seen in the SV ALMA polarization data (see the CASA Guide).
		   Default Value: 0.001

		tau0:	Atmospheric opacity at zenith.
		   Default Value: 0.0

		t_sky:	Sky temperature (in K).
		   Default Value: 250.0

		t_ground:	Ground temperature (in K).
		   Default Value: 270.0

		t_receiver:	Receiver temperature (in K).
		   Default Value: 50.0

	Returns: bool

	Example :


  EXAMPLES OF USE:

  1.- Simulate only a point source with constant fractional polarization of 10%,
  an EVPA of zero degrees, no V Stokes, and a steep spectrum:
  
  model_image=[]
  I = [1.0], Q = [0.1], U = [0.], V = [0.], RM = [0.], spec_index = [-1.0]
  spectrum_file = ''


  
  2.- Simulate an extended source, defined by a set of images
  (one image for each Stokes parameter):

  model_image=['I.image', 'Q.image', 'U.image', 'V.image']
  I = [], Q = [], U = [], RM = [], spec_index = []
  spectrum_file=''

  
  
  3.- Simulate the TWO SOURCES of the previous examples
  TOGETHER (i.e., one OVER THE OTHER):

  model_image=['I.image', 'Q.image', 'U.image']
  I = [1.0], Q = [0.1], U = [0.], V = [0.], RM = [0.], spec_index = [-1.0]
  spectrum_file = ''



  4.- Simulate a point source with a user-defined 
  polarization spectrum for I, Q, U, and V:

  model_image=[]
  I = [], Q = [], U = [], V = [], RM = [], spec_index = []
  spectrum_file='my_spectrum.dat'

  The format of my_spectrum is an ASCII file with several rows 
  (one row per frequency channel). The format of each row is:

  FREQUENCY (HZ)    I (Jy)    Q (Jy)    U (Jy)  V (Jy)
  Some example rows (for a Band 6 simulation) could be:

  246.e9   1.0   0.07   0.05   0.00
  247.e9   1.2   0.06   0.05   0.00 
  248.e9   1.3   0.05   0.04   0.00
  ...

  
  The spectrum will be INTERPOLATED to the frequencies of the 
  spectral windows in the measurement set. BEWARE and check that 
  the spw frequencies are properly covered by the spectrum defined 
  in the spectrum file!



  5.- Simulate a source with two compact component of different 
  polarization. The first component will have: I=1, Q=0.03, U=0.02;
  the second component will have I=0.4, Q=0.02, U=-0.02, with an
  offset of 0.5 mas in RA and 0.2 mas in Dec:

  I = [1.00,  0.40]
  Q = [0.03,  0.02]
  U = [0.02, -0.02]
  V = [0.00,  0.00]
  RAoffset =  [0.0, 0.0005]
  Decoffset = [0.0, 0.0002]
  RM = [0.0, 0.0]
  spec_index = [0.0, 0.0]

  The source will be observed with the VLBA at 5GHz, using two spectral
  windows (IFs) of 64 channels and 512MHz each (one centered at 4744 MHz 
  and the other centered at 5256 MHz):

  LO = 5.e9
  BBs = [4.744e9, 5.256e9]
  spw_width = 512.e6
  nchan = 64


  If you want to mimick the exact uv coverage of a real observation
  (or if you want to have a full control over the observing schedule) 
  you can set \'nscan\' to the name of an ascii file, expected to have 
  the time and scan information in \'listobs\' format. 
  For instance, the contents of that file could be:

  ###################################

  Observed from   11-Apr-2017/01:02:30.0   to   11-Apr-2017/08:45:00.0 (UTC)

  11-Apr-2017/01:02:30.0 - 01:08:30.0     
              01:38:26.8 - 01:44:00.0    
              06:51:00.0 - 07:06:00.0    


  ###################################

  and so on.

  Notice that, in this case, the values of \'observe_time\', \'onsource_time\' 
  and \'H0\' will be ignored.



        """
	if not hasattr(self, "__globals__") or self.__globals__ == None :
           self.__globals__=stack_frame_find( )
	#casac = self.__globals__['casac']
	casalog = self.__globals__['casalog']
	casa = self.__globals__['casa']
	#casalog = casac.casac.logsink()
        self.__globals__['__last_task'] = 'polsimulate'
        self.__globals__['taskname'] = 'polsimulate'
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
            myparams['reuse'] = reuse = self.parameters['reuse']
            myparams['array_configuration'] = array_configuration = self.parameters['array_configuration']
            myparams['elevation_cutoff'] = elevation_cutoff = self.parameters['elevation_cutoff']
            myparams['feed'] = feed = self.parameters['feed']
            myparams['mounts'] = mounts = self.parameters['mounts']
            myparams['ConstDt0'] = ConstDt0 = self.parameters['ConstDt0']
            myparams['ConstDt1'] = ConstDt1 = self.parameters['ConstDt1']
            myparams['LO'] = LO = self.parameters['LO']
            myparams['BBs'] = BBs = self.parameters['BBs']
            myparams['spw_width'] = spw_width = self.parameters['spw_width']
            myparams['nchan'] = nchan = self.parameters['nchan']
            myparams['model_image'] = model_image = self.parameters['model_image']
            myparams['I'] = I = self.parameters['I']
            myparams['Q'] = Q = self.parameters['Q']
            myparams['U'] = U = self.parameters['U']
            myparams['V'] = V = self.parameters['V']
            myparams['RM'] = RM = self.parameters['RM']
            myparams['spec_index'] = spec_index = self.parameters['spec_index']
            myparams['RAoffset'] = RAoffset = self.parameters['RAoffset']
            myparams['Decoffset'] = Decoffset = self.parameters['Decoffset']
            myparams['spectrum_file'] = spectrum_file = self.parameters['spectrum_file']
            myparams['phase_center'] = phase_center = self.parameters['phase_center']
            myparams['incell'] = incell = self.parameters['incell']
            myparams['inbright'] = inbright = self.parameters['inbright']
            myparams['inwidth'] = inwidth = self.parameters['inwidth']
            myparams['innu0'] = innu0 = self.parameters['innu0']
            myparams['H0'] = H0 = self.parameters['H0']
            myparams['onsource_time'] = onsource_time = self.parameters['onsource_time']
            myparams['observe_time'] = observe_time = self.parameters['observe_time']
            myparams['visib_time'] = visib_time = self.parameters['visib_time']
            myparams['nscan'] = nscan = self.parameters['nscan']
            myparams['apply_parang'] = apply_parang = self.parameters['apply_parang']
            myparams['export_uvf'] = export_uvf = self.parameters['export_uvf']
            myparams['corrupt'] = corrupt = self.parameters['corrupt']
            myparams['seed'] = seed = self.parameters['seed']
            myparams['Dt_amp'] = Dt_amp = self.parameters['Dt_amp']
            myparams['Dt_noise'] = Dt_noise = self.parameters['Dt_noise']
            myparams['tau0'] = tau0 = self.parameters['tau0']
            myparams['t_sky'] = t_sky = self.parameters['t_sky']
            myparams['t_ground'] = t_ground = self.parameters['t_ground']
            myparams['t_receiver'] = t_receiver = self.parameters['t_receiver']


	result = None

#
#    The following is work around to avoid a bug with current python translation
#
        mytmp = {}

        mytmp['vis'] = vis
        mytmp['reuse'] = reuse
        mytmp['array_configuration'] = array_configuration
        mytmp['elevation_cutoff'] = elevation_cutoff
        mytmp['feed'] = feed
        mytmp['mounts'] = mounts
        mytmp['ConstDt0'] = ConstDt0
        mytmp['ConstDt1'] = ConstDt1
        mytmp['LO'] = LO
        mytmp['BBs'] = BBs
        mytmp['spw_width'] = spw_width
        mytmp['nchan'] = nchan
        mytmp['model_image'] = model_image
        mytmp['I'] = I
        mytmp['Q'] = Q
        mytmp['U'] = U
        mytmp['V'] = V
        mytmp['RM'] = RM
        mytmp['spec_index'] = spec_index
        mytmp['RAoffset'] = RAoffset
        mytmp['Decoffset'] = Decoffset
        mytmp['spectrum_file'] = spectrum_file
        mytmp['phase_center'] = phase_center
        mytmp['incell'] = incell
        mytmp['inbright'] = inbright
        mytmp['inwidth'] = inwidth
        mytmp['innu0'] = innu0
        mytmp['H0'] = H0
        mytmp['onsource_time'] = onsource_time
        mytmp['observe_time'] = observe_time
        mytmp['visib_time'] = visib_time
        mytmp['nscan'] = nscan
        mytmp['apply_parang'] = apply_parang
        mytmp['export_uvf'] = export_uvf
        mytmp['corrupt'] = corrupt
        mytmp['seed'] = seed
        mytmp['Dt_amp'] = Dt_amp
        mytmp['Dt_noise'] = Dt_noise
        mytmp['tau0'] = tau0
        mytmp['t_sky'] = t_sky
        mytmp['t_ground'] = t_ground
        mytmp['t_receiver'] = t_receiver
	pathname="file:///home/marti/WORKAREA/ARC_TOOLS/CASA-PolTools/casa-poltools/"
	trec = casac.casac.utils().torecord(pathname+'polsimulate.xml')

        casalog.origin('polsimulate')
	try :
          #if not trec.has_key('polsimulate') or not casac.casac.utils().verify(mytmp, trec['polsimulate']) :
	    #return False

          casac.casac.utils().verify(mytmp, trec['polsimulate'], True)
          scriptstr=['']
          saveinputs = self.__globals__['saveinputs']

          # Save .last file for this task execution. MPI servers don't write it (CASR-329).
          from mpi4casa.MPIEnvironment import MPIEnvironment
          do_full_logging = MPIEnvironment.is_mpi_disabled_or_client()
          if type(self.__call__.func_defaults) is NoneType:
              saveinputs=''
          else:
              saveinputs('polsimulate', 'polsimulate.last', myparams, self.__globals__,scriptstr=scriptstr, do_save_inputs=do_full_logging)

          tname = 'polsimulate'
          spaces = ' '*(18-len(tname))
          casalog.post('\n##########################################'+
                       '\n##### Begin Task: ' + tname + spaces + ' #####')
          # Don't do telemetry from MPI servers (CASR-329)
          if do_full_logging and casa['state']['telemetry-enabled']:
              #casalog.poststat('Begin Task: ' + tname)
              task_starttime = str(datetime.datetime.now())
          if type(self.__call__.func_defaults) is NoneType:
              casalog.post(scriptstr[0]+'\n', 'INFO')
          else:
              casalog.post(scriptstr[1][1:]+'\n', 'INFO')

          # Effective call to the task as defined in gcwrap/python/scripts/task_*
          result = polsimulate(vis, reuse, array_configuration, elevation_cutoff, feed, mounts, ConstDt0, ConstDt1, LO, BBs, spw_width, nchan, model_image, I, Q, U, V, RM, spec_index, RAoffset, Decoffset, spectrum_file, phase_center, incell, inbright, inwidth, innu0, H0, onsource_time, observe_time, visib_time, nscan, apply_parang, export_uvf, corrupt, seed, Dt_amp, Dt_noise, tau0, t_sky, t_ground, t_receiver)

          if do_full_logging and casa['state']['telemetry-enabled']:
              task_endtime = str(datetime.datetime.now())
              casalog.poststat( 'Task ' + tname + ' complete. Start time: ' + task_starttime + ' End time: ' + task_endtime )
          casalog.post('##### End Task: ' + tname + '  ' + spaces + ' #####'+
                       '\n##########################################')

	except Exception, instance:
          if(self.__globals__.has_key('__rethrow_casa_exceptions') and self.__globals__['__rethrow_casa_exceptions']) :
             raise
          else :
             #print '**** Error **** ',instance
	     tname = 'polsimulate'
             casalog.post('An error occurred running task '+tname+'.', 'ERROR')
             pass
	casalog.origin('')

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
#        paramgui.runTask('polsimulate', myf['_ip'])
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
        a['vis']  = 'polsimulate_output.ms'
        a['reuse']  = False
        a['array_configuration']  = 'alma.out04.cfg'
        a['elevation_cutoff']  = 5.0
        a['feed']  = 'linear'
        a['mounts']  = []
        a['ConstDt0']  = []
        a['ConstDt1']  = []
        a['LO']  = 100.e9
        a['BBs']  = [-7.e9, -5.e9, 5.e9, 7.e9]
        a['spw_width']  = 2.e9
        a['nchan']  = 8
        a['model_image']  = []
        a['I']  = []
        a['Q']  = []
        a['U']  = []
        a['V']  = []
        a['RM']  = []
        a['spec_index']  = []
        a['RAoffset']  = []
        a['Decoffset']  = []
        a['spectrum_file']  = ''
        a['phase_center']  = 'J2000 00h00m00.00 -00d00m00.00'
        a['incell']  = ''
        a['inbright']  = 0.0
        a['inwidth']  = ''
        a['innu0']  = ''
        a['H0']  = -1.5
        a['onsource_time']  = 1.0
        a['observe_time']  = 3.0
        a['visib_time']  = '6s'
        a['nscan']  = 50
        a['apply_parang']  = False
        a['export_uvf']  = True
        a['corrupt']  = True
        a['seed']  = 42
        a['Dt_amp']  = 0.0
        a['Dt_noise']  = 0.001
        a['tau0']  = 0.0
        a['t_sky']  = 250.0
        a['t_ground']  = 270.0
        a['t_receiver']  = 50.0


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
    def description(self, key='polsimulate', subkey=None):
        desc={'polsimulate': 'Version 1.3.2 - Basic simulator of ALMA/J-VLA (and VLBI) full-polarization observations. The output should be imaged with CLEAN (with stokes=IQUV) and the polarization vectors should be computed with immath (with options poli and pola). See the ALMA Polarization CASA Guide for more information.\n\n',
               'vis': 'Name of output measurement set.',
               'reuse': 'If True, and the measurement set exists, will reuse it (so the antenna names and source coordinates will be reused as well).',
               'array_configuration': 'Array configuration file, where the antenna coordinates and diameters are set (see the files in data/alma/simmos of your CASA path). Default is an ALMA configuration. The name of the array is taken as the part of the name of this file before the first dot (i.e., it is \'alma\' for the default filename).',
               'elevation_cutoff': 'Minimum allowed elevation for the antennas (degrees).',
               'feed': 'Polarization basis for the measurement set. Can be linear (e.g., for ALMA) or circular (e.g., for VLA). Default is linear.',
               'mounts': 'For VLBI observations, this is a list of the antenna mounts (given in the same order as in the array configuration file). A mount type is specified with two characters. Supported mounts: alt-az (\'AZ\'), equatorial (\'EQ\'), X-Y (\'XY\'), Nasmyth Right (\'NR\') and Nasmyth Left (\'NL\'). Default means all antennas are alt-az.',
               'ConstDt0': 'List of complex numbers (length equal to the number of antennas). If not empty, the first polarizer (i.e., either R or X, depending on the value of \'feed\') will be contamined with a leakage given by these Dterms.',
               'ConstDt1': 'List of complex numbers (length equal to the number of antennas). If not empty, the second polarizer will be contamined with a leakage given by these Dterms.',
               'LO': 'Frequency of the first LO in Hz (this will define the ALMA band of observation). Default is 100 GHz (i.e., ALMA Band 3)',
               'BBs': 'List with the baseband frequency offsets (in Hz). There will be one spectral window per baseband. ',
               'spw_width': 'Width of the spectral windows in Hz (the spws will be centered at each BB).',
               'nchan': 'Number of channels in each spw (all spws will have the same number of channels).',
               'model_image': 'List of four images (Stokes I, Q, U, and V, respectively) to be used as observed extended sources. Image cubes are allowed. Default is to NOT simulate extended sources. BEWARE OF THE CURRENT ALMA LIMITATION for extended sources (should fall within the inner 1/3 of the primary beam FWHM).',
               'I': 'List of Stokes I (in Jy) for a set of point sources to simulate. These are added to the \'model_image\' (if it is provided). Default is to add NO source. The flux density is referred to the LO frequency. Example: [1.0] for a 1 Jy source.',
               'Q': 'List of Stokes Q (in Jy) for the sources defined above. Default is no source. Example: [0.0] for no Q signal. These values are referred to the LO frequency.',
               'U': '>List of Stokes U (in Jy) for the sources defined above. Default is no source. Example: [0.0] for no U signal. These values are referred to the LO frequency.',
               'V': '>List of Stokes V (in Jy) for the sources defined above. Default is no source. Example: [0.0] for no V signal. These values are referred to the LO frequency.',
               'RM': 'List of Rotation Measures (RM, in rad/m**2.) for the sources defined above. Default is no source. Example: [0.0] for no RM.',
               'spec_index': 'List of spectral indices for the sources defined above. Default is no source. Example: [0.0] for a flat spectrum.',
               'RAoffset': 'List of right-ascension offsets (in arcsec) for the sources defined above. The first source is assumed to be at the phase center, so all sources will be shifted RAoffset[0] arcsec (so that the first source in the list is at the phase center).',
               'Decoffset': 'List of declination offsets (in arcsec) for the sources defined above. The first source is assumed to be at the phase center, so all sources will be shifted Decoffset[0] arcsec.',
               'spectrum_file': 'File with user-defined spectra of I, Q, U, and V. See help for details about the file format. This source WILL BE ADDED TO THE PHASE CENTER, together with the source defined in the model_image model and all those given in the I, Q, U, V, lists.',
               'phase_center': 'Coordinates of the observed source (will override the coordinates defined in model_image, if an image is being used). This keyword MUST BE defined.',
               'incell': 'Pixel size of the model_image. If not empty, will override the original value stored in the image. Example: \'0.1arcsec\'. All the Stokes images (I, Q, U, and V) will be set the same way. USE WITH CARE.',
               'inbright': 'Peak intensity of the I Stokes model_image (will override the original value stored in the image). Default is to use the original brightness unit. All the Stokes images (I, Q, U, and V) will be set the same way. Default (i.e., 0.0) means to NOT scale the image birghtness. USE WITH CARE.',
               'inwidth': 'Width of the frequency channels in the model_image. If not empty, will override the original value stored in the image. Example: \'10MHz\'. All the Stokes images (I, Q, U, and V) will be set the same way.',
               'innu0': 'Frequency of the first image channel (e.g., \'1GHz\'). Default (empty) means to use the value in the image.',
               'H0': 'If a float is given, it is the Hour Angle at the start of the observations, as observed from the array center (in hr). If a string is given, it is the exact start of the observations (UT time) in the format \'2017/01/01/00:00:00\'.',
               'onsource_time': 'Total effective integration time on the source (in hr). Default is 1 hr.',
               'observe_time': 'Total observing time (i.e., including overheads) in hr. Default is 3h, so there will be an observing efficiency of 0.33 if \'onsource_time\' is set to one hour. ',
               'visib_time': 'Integration time per visibility. This is a string, where \'s\' stands for \'seconds\'.',
               'nscan': 'Number of scans. Can be provided as a \'listobs\' file \n(in that case, observe_time, onsource_time \n and H0 are not used, but taken from the listobs).\n If a list in \'listobs\' format is not provided \nthen all scans will be set to equal length. \n If just an integer is given, the scans will be \n homogeneously distributed across the \n total observing time. If a list is given, \n the values will be taken as the starting times \n of the scans, relative to the \n duration of the experiment. For instance, if \n \'observe_time = 6.0\' then \'nscan = [0., 0.5, 0.75]\' \n will make three scans, one at the start of the \n observations, the second one 3 hours later \n and the third one 4.5 hours after \n the start of the observations.',
               'apply_parang': 'If True, applies the parallactic-angle correction. \n If False, the data polarization will be given \n in the antenna frame (i.e., just as true raw data).',
               'export_uvf': 'If True, exports the measurement into uvfits format (for its use in e.g., AIPS/Difmap).',
               'corrupt': 'Whether to add random noise to the visibilities.',
               'seed': 'Seed of the random number generator in the sm tool.',
               'Dt_amp': 'Will add random Dterms (antenna-wise). Dt_amp is the typical absolute value of the Dterms (real and imag). The actual values will be computed from a random Gaussian distribution.',
               'Dt_noise': 'Will add random channel-dependent contribution to the Dterms. Dt_noise is the typical residual channel noise in the Dterms (real and imag). The actual values for each frequency channel will be those of Dt_amp PLUS a random Gaussian distribution of width Dt_noise. Default is 0.001, similar to the spectral spread of Dterms seen in the SV ALMA polarization data (see the CASA Guide).',
               'tau0': 'Atmospheric opacity at zenith.',
               't_sky': 'Sky temperature (in K).',
               't_ground': 'Ground temperature (in K).',
               't_receiver': 'Receiver temperature (in K).',

              }

        if(desc.has_key(key)) :
           return desc[key]

    def itsdefault(self, paramname) :
        a = {}
        a['vis']  = 'polsimulate_output.ms'
        a['reuse']  = False
        a['array_configuration']  = 'alma.out04.cfg'
        a['elevation_cutoff']  = 5.0
        a['feed']  = 'linear'
        a['mounts']  = []
        a['ConstDt0']  = []
        a['ConstDt1']  = []
        a['LO']  = 100.e9
        a['BBs']  = [-7.e9, -5.e9, 5.e9, 7.e9]
        a['spw_width']  = 2.e9
        a['nchan']  = 8
        a['model_image']  = []
        a['I']  = []
        a['Q']  = []
        a['U']  = []
        a['V']  = []
        a['RM']  = []
        a['spec_index']  = []
        a['RAoffset']  = []
        a['Decoffset']  = []
        a['spectrum_file']  = ''
        a['phase_center']  = 'J2000 00h00m00.00 -00d00m00.00'
        a['incell']  = ''
        a['inbright']  = 0.0
        a['inwidth']  = ''
        a['innu0']  = ''
        a['H0']  = -1.5
        a['onsource_time']  = 1.0
        a['observe_time']  = 3.0
        a['visib_time']  = '6s'
        a['nscan']  = 50
        a['apply_parang']  = False
        a['export_uvf']  = True
        a['corrupt']  = True
        a['seed']  = 42
        a['Dt_amp']  = 0.0
        a['Dt_noise']  = 0.001
        a['tau0']  = 0.0
        a['t_sky']  = 250.0
        a['t_ground']  = 270.0
        a['t_receiver']  = 50.0

        #a = sys._getframe(len(inspect.stack())-1).f_globals

        if a.has_key(paramname) :
	      return a[paramname]
polsimulate_cli = polsimulate_cli_()
