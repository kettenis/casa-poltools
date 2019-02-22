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
import task_polsimulate
def polsimulate(vis='polsimulate_output.ms', reuse=False, array_configuration='alma.out04.cfg', elevation_cutoff=5.0, feed='linear', mounts=[], ConstDt0=[], ConstDt1=[], LO=100.e9, BBs=[-7.e9, -5.e9, 5.e9, 7.e9], spw_width=2.e9, nchan=8, model_image=[], I=[], Q=[], U=[], V=[], RM=[], spec_index=[], RAoffset=[], Decoffset=[], spectrum_file='', phase_center='J2000 00h00m00.00 -00d00m00.00', incell='', inbright=0.0, inwidth='', innu0='', H0=-1.5, onsource_time=1.0, observe_time=3.0, visib_time='6s', nscan=50, apply_parang=False, export_uvf=True, corrupt=True, seed=42, Dt_amp=0.0, Dt_noise=0.001, tau0=0.0, t_sky=250.0, t_ground=270.0, t_receiver=50.0):

        """Version 1.3.2 - Basic simulator of ALMA/J-VLA (and VLBI) full-polarization observations. The output should be imaged with CLEAN (with stokes=IQUV) and the polarization vectors should be computed with immath (with options poli and pola). See the ALMA Polarization CASA Guide for more information.\n\n

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
	pathname="file:///data/SHARED/WORKAREA/ARC_TOOLS/CASA-PolTools/trunk/"
	trec = casac.utils().torecord(pathname+'polsimulate.xml')

        casalog.origin('polsimulate')
        if trec.has_key('polsimulate') and casac.utils().verify(mytmp, trec['polsimulate']) :
	    result = task_polsimulate.polsimulate(vis, reuse, array_configuration, elevation_cutoff, feed, mounts, ConstDt0, ConstDt1, LO, BBs, spw_width, nchan, model_image, I, Q, U, V, RM, spec_index, RAoffset, Decoffset, spectrum_file, phase_center, incell, inbright, inwidth, innu0, H0, onsource_time, observe_time, visib_time, nscan, apply_parang, export_uvf, corrupt, seed, Dt_amp, Dt_noise, tau0, t_sky, t_ground, t_receiver)

	else :
	  result = False
        return result
