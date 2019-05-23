# PolSolve - A task to solve fot Dterms with extended calibrator sources.
#
# Copyright (c) Ivan Marti-Vidal - Observatorio de Yebes (2018). 
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>,
# or write to the Free Software Foundation, Inc., 
# 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# a. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# b. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the
#    distribution.
# c. Neither the name of the author nor the names of contributors may 
#    be used to endorse or promote products derived from this software 
#    without specific prior written permission.
#
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#


# NOTES: Only supports circular feeds.
#        Antenna mounts always have to be given explicitely.
#        2nd order corrections (slower).
#        Parangle for orbital elements is not implemented.

# TODO: 
#       1. Combine IFs and parameterize frequency dependence.
#
#       2. Work directly on FITS-IDI (or uvfits) files.
#
#       3. Write a more complete documentation.
#
#       4. Solve for RM.
#
#       5. Apply the calibration into the corrected column.   DONE!
#       (to overcome mounts limitations in CASA).
#
#       6. Multi-source self-consistent fitting.
#
#

import gc
import os
import numpy as np
import pylab as pl
from taskinit import gentools
from scipy.optimize import minimize
import PolSolver as PS 
from taskinit import casalog
import sys

ms = gentools(['ms'])[0]
tb = gentools(['tb'])[0]

__version__ = '1.1b'


#####################
# UNIT TEST LINES:
if __name__=='__main__':

  vis                =  "3C279_3601_UVFIX2.ms"
  spw                =  0
  field              =  "0"
  mounts             =  ['AZ', 'NR', 'NR', 'AZ', 'NL', 'NL', 'NL', 'AZ', 'NL']
  feed_rotation      =  []
  DR                 =  []
  DL                 =  []
  DRSolve            =  [True, True, True, False, True, True, True, True, False]
  DLSolve            =  [True, True, True, False, True, True, True, True, False]
  CLEAN_models       =  ['3C279_3601_image.fits_CC00.dat', '3C279_3601_image.fits_CC01.dat']
  Pfrac              =  [0.0, 0.0]
  EVPA               =  [0.0, 0.0]
  PolSolve           =  [True, True]
  parang_corrected    =  True
  target_field       =  ""
  plot_parang        =  True
  min_elev_plot      =  0.0
  wgt_power          =  1.0
#
#
#
##################






def polsolve(vis = 'input.ms', spw=0, field = '0', mounts = [], feed_rotation = [], DR = [], DL = [],
                DRSolve = [], DLSolve = [], CLEAN_models = [1.0], Pfrac = [0.0], 
                EVPA = [0.0], PolSolve = [True], parang_corrected = True, target_field = '', plot_parang=False, min_elev_plot=10.0, wgt_power=1.0):



########
# Uncomment for testing/debugging
#  return True
#if __name__=='__main__':
########



  STOKES_CODES = {'I': 1,  'Q': 2,  'U': 3,  'V': 4, 
                 'RR': 5, 'RL': 6, 'LR': 7, 'LL': 8,
                 'XX': 9, 'XY':10, 'YX':11, 'YY':12}

  POL_ORDER = ['RR','RL','LR','LL']

###########################################
# HELPER FUNCTIONS
###########################################


# Helper functions to print info/errors in terminal + logger:

  def printError(msg):
    print '\n', msg, '\n'
    casalog.post('PolSolve: '+msg)
    raise Exception(msg)

  def printMsg(msg,dolog=True):
    print msg
    if dolog:
      casalog.post('PolSolve: '+msg)





# Helper function (return list of target fields, given in CASA format):

  def getFields(vis,field_str):
    """ Takes a string of field range selection and returns a list 
    with the field ids. Examples of strings:
    
    '0~3'
    '0,1,2,3~5'
    'M87,3C279'
    """
      
    tb.open(os.path.join(vis,'FIELD'))  
    NAMES = list(tb.getcol('NAME'))
    
    # Remove annoying spaces:
    fields = field_str.replace(' ','')
    
    # Field entries separated by commas:
    selFields = []
    for fld in fields.split(','):
      if '~' in fld: # Field range (must be integers)
        try:  
          f0,f1 = map(int,fld.split('~'))  
        except:
          printError('WRONG FORMAT FOR TARGET FIELDS!')
        if f0<0 or f1>len(NAMES): # Field ids must exist!
          printError('WRONG RANGE OF TARGET FIELDS!')  
        selFields += range(f0,f1+1)
      elif fld in NAMES: # Field name is given.
        selFields += [NAMES.index(fld)]
      else: # Field id is given (must be an integer).
        try:
          selFields += [int(fld)]
        except:
          printError('WRONG FIELD %s!'%fld)
          
    return selFields





# Return GMST from UT time (from NRAO webpage):
  def GMST(MJD):
    Days = MJD/86400.  
    t = (Days -51544.0)/36525.
    Hh = (Days - np.floor(Days))
    GMsec = 24110.54841 + 8640184.812866*t + 0.093104*t*t - 0.0000062*t*t*t
    return (GMsec/86400. + Hh)*2.*np.pi



# Helper function to compute feed (i.e., mount + parallactic) angles:

  def getParangle(vis,spw,field,mounts, feedAngles, scan = -1, doplot=False):
      
    """ Returns feed rotation (mount + parangle) for the visibilities
    corresponding to field (id) or to the scan number (if >= 0). """
      
      
    tb.open(os.path.join(vis,'FIELD'))
    scoord = tb.getcol('PHASE_DIR')
    tb.close()      

# NOTE: Revise this for multi-source measurement sets!!
    RA = scoord[0,0,field]
    Dec = scoord[1,0,field]

    CosDec = np.cos(Dec)
    SinDec = np.sin(Dec)
    
# Load antenna info:
    tb.open(os.path.join(vis,'ANTENNA'))
    apos = tb.getcol('POSITION')
    nant = len(apos)
    tb.close()

    Lat = np.arctan2(apos[2,:],np.sqrt(apos[0,:]**2. + apos[1,:]**2.))
    Tlat = np.tan(Lat)
    Lon = np.arctan2(apos[1,:],apos[0,:])


# Load data:
    ms.open(vis)
    ms.selectinit(datadescid=spw)
    if scan < 0:
      ms.select({'field_id':field})
    else:
      ms.select({'scan_number':scan})
      
    DATA = ms.getdata(['u','v','w','antenna1','antenna2', 'time'])
    Nvis  = len(DATA['u'])
    
    ms.close()

# Compute PAs using UV coordinates:

    PAs = np.zeros((Nvis,2))
    
    GOODS = np.ones(Nvis,dtype=np.bool)

 #   V2 = SinDec*DATA['v'] - CosDec*DATA['w']
    
 #   Bx = -(apos[0,DATA['antenna2']]-apos[0,DATA['antenna1']])
 #   By = -(apos[1,DATA['antenna2']]-apos[1,DATA['antenna1']])
 #   Bz = -(apos[2,DATA['antenna2']]-apos[2,DATA['antenna1']])

 #   CH = DATA['u']*By - V2*Bx
 #   SH = DATA['u']*Bx + V2*By
    
    CT1 = CosDec*Tlat[DATA['antenna1']]
    CT2 = CosDec*Tlat[DATA['antenna2']]
    
 #   HAng2 = np.arctan2(SH,CH)
    HAng = (GMST(DATA['time']) - RA)%(2.*np.pi)

 #   print '\n'
 #   print DATA['time'][10],GMST(DATA['time'])[10] , HAng[10], HAng2[10]
 #   print DATA['time'][50], GMST(DATA['time'])[50], HAng[50], HAng2[50]


    H1 = HAng + Lon[DATA['antenna1']]
    H2 = HAng + Lon[DATA['antenna2']]
    
    
 #   Autos = (CH==0.)*(SH==0.)
    Autos = DATA['antenna1']==DATA['antenna2']

    H1[Autos] = 0.0
    H2[Autos] = 0.0
    
    GOODS[Autos] = False

    E1 = np.arcsin(SinDec*np.sin(Lat[DATA['antenna1']])+np.cos(Lat[DATA['antenna1']])*CosDec*np.cos(H1))
    E2 = np.arcsin(SinDec*np.sin(Lat[DATA['antenna2']])+np.cos(Lat[DATA['antenna2']])*CosDec*np.cos(H2))

    GOODS[E1<min_elev_plot*np.pi/180.] = False
    GOODS[E2<min_elev_plot*np.pi/180.] = False
    GOODS[E1>np.pi/2.] = False
    GOODS[E2>np.pi/2.] = False


    PAZ1 = np.arctan2(np.sin(H1), CT1 - SinDec*np.cos(H1))
    PAZ2 = np.arctan2(np.sin(H2), CT2 - SinDec*np.cos(H2))


# Add mount rotations:

    for mti,mt in enumerate(mounts):
      filt = DATA['antenna1'] == mti 

      if mt=='AZ':
        PAs[filt,0] = -PAZ1[filt]
      elif mt=='EQ':
        PAs[filt,0] = 0.0
      elif mt=='XY':
        PAs[filt,0] = -np.arctan2(np.cos(H1[filt]),SinDec*np.sin(H1[filt]))
      elif mt=='NR':
        PAs[filt,0] = -PAZ1[filt] - E1[filt]
      elif mt=='NL':
        PAs[filt,0] = -PAZ1[filt] + E1[filt]

      PAs[filt,0] += feedAngles[mti]

      filt = DATA['antenna2'] == mti 
      
      if mt=='AZ':
        PAs[filt,1] = -PAZ2[filt]
      elif mt=='EQ':
        PAs[filt,1] = 0.0
      elif mt=='XY':
        PAs[filt,1] = -np.arctan2(np.cos(H2[filt]),SinDec*np.sin(H2[filt]))
      elif mt=='NR':
        PAs[filt,1] = -PAZ2[filt] - E2[filt]
      elif mt=='NL':
        PAs[filt,1] = -PAZ2[filt] + E2[filt]

      PAs[filt,1] += feedAngles[mti]

# Release memory:
    if doplot:
      TT = np.copy(DATA['time']); A1 = np.copy(DATA['antenna1']); A2 = np.copy(DATA['antenna2'])

    del DATA['antenna1'], DATA['antenna2'], DATA['u'], DATA['v'], DATA['w'], DATA['time']
#    del H1, H2, E1, E2, PAZ1, PAZ2, Bx, By, Bz, CH, SH, CT1, CT2, V2, filt
    del H1, H2, E1, E2, PAZ1, PAZ2, CT1, CT2, filt

# Finished!
    if doplot:
      return [PAs, TT, A1, A2, GOODS]
    else:
      return PAs






# Helper function to print matrices nicely:
  def printMatrix(a,f="% 6.3f"):
    print "Matrix["+("%d" %a.shape[0])+"]["+("%d" %a.shape[1])+"]"
    rows = a.shape[0]
    cols = a.shape[1]
    for i in range(0,rows):
      for j in range(0,cols):
         print(f%a[i,j]) 
      print
    print



# Helper function that fills the vector of variables,
# given the vector of fitting parameters:
  def setFitVal(p,NmodComp, nant, FITSOU, DRSolve, DLSolve):
    """ Fills vector of variables given the fitting parameters.
    First, Stokes Q and U for all fittable components. 
    Then, the fittable DR and DL (in real/imag base). 
    
    Fittable source components are given by the list of booleans FITSOU.
    Fittable Dterms are given by the lists of booleans DRSolve and DLSolve.
    NmodComp and nant is the total number (i.e., fixed and fittable) of 
    source components and antennas, respectively.
    """

    i = 0
    for si in range(NmodComp):
      if FITSOU[si]:
        FitVal[2*si  ]  = p[i] 
        FitVal[2*si+1]  = p[i+1] 
        i += 2
    for ni in range(nant):
      if DRSolve[ni]:
        FitVal[2*NmodComp+2*ni  ]  = p[i]
        FitVal[2*NmodComp+2*ni+1]  = p[i+1]
        i += 2
    for ni in range(nant):
      if DLSolve[ni]:
        FitVal[2*(NmodComp+nant)+2*ni  ]  = p[i]
        FitVal[2*(NmodComp+nant)+2*ni+1]  = p[i+1]
        i += 2


# Computes the Chi2, given a set of fitting parameter values (p)
# and also computes the Hessian and the residuals vector:
  def getChi2(p, NmodComp, nant, FITSOU, DRSolve, DLSolve, arg=0):
    """ Returns either the Chi square and Hessian (+ residuals) or
    the optimum model flux scaling, depending on the value of \"arg\" 
    (see documentation of C++ module). Same arguments as setFitVal."""
    
    setFitVal(p,NmodComp,nant, FITSOU, DRSolve, DLSolve)
    return PS.getHessian(0.0)[arg]









###########################################
# END OF HELPER FUNCTIONS
###########################################












# SCRIPT STARTS!


  printMsg( '\n\n  POLSOLVE - VERSION %s  - I. Marti-Vidal (Yebes Observatory, Spain)'%__version__)

  GoodMounts = ['AZ','EQ','NR','NL','XY']

  printMsg('WARNING! Currently, PolSolve only handles circular-feed receivers')

# Currently, we only work with circular feeds:
  doCirc = True


# Load the source model components:

# First, some sanity checks:
  if not (type(CLEAN_models) is list):
    CLEAN_models = [CLEAN_models]  

  if len(set([len(CLEAN_models),len(Pfrac),len(EVPA),len(PolSolve)]))>1:
    printError("ERROR! CLEAN_models, Pfrac, EVPA, and PolSolve should all have the same length!")

# How many components do we have?
  NmodComp = max([len(CLEAN_models),1])

# Flux density and offset (RA, Dec) of the deltas for each component:
  DELTAS = [[] for dd in CLEAN_models]
  for mi in range(NmodComp):

# Case of AIPS-like CC files:
    if type(CLEAN_models[mi]) is str:
      if not os.path.exists(CLEAN_models[mi]):
         printError('ERROR! File %s does not exist!'%CLEAN_models[mi])
      
      ifile = open(CLEAN_models[mi])
      for line in ifile.readlines():
        spl = line.split()
        if len(line)>0 and len(spl)==4:
          try:
            test = int(spl[0])  
            DELTAS[mi].append(map(float,line.split()[1:4]))
          except:
            pass
      ifile.close()
      if len(DELTAS[mi])==0:
        printError('ERROR! Problem with CLEAN model file %s'%CLEAN_models[mi])


# Case of a centered point source:
    else:
      DELTAS[mi].append([float(CLEAN_models[mi]),0.,0.])

# Polarization of each component:
  PFRAC = map(float,Pfrac)
  POLANG = map(float,EVPA)

# Which components have a fittable polarization?
  FITSOU = map(bool,PolSolve)

# Print for checking:
  printMsg('\nThere are %i model components'%NmodComp)
  TotFluxAll = 0.0
  for sou in range(NmodComp):
    TotFlux = np.sum([di[0] for di in DELTAS[sou]])
    printMsg('Comp. %i is made of %i deltas (%.3f Jy). A-priori pol: %.2e (%.1f deg.)'%(sou+1,len(DELTAS[sou]),TotFlux,PFRAC[sou], POLANG[sou]))
    TotFluxAll += TotFlux


# Sanity checks:
  if not ms.open(vis):
      printError('ERROR with input measurement set!')

  if not ms.selectinit(datadescid=int(spw)):
      printError('ERROR with spw %i!'%spw)
  ms.close()


# Load info about calibrator:
  tb.open(os.path.join(vis,'FIELD'))
  snam = list(tb.getcol('NAME'))
  tb.close()
  if field not in snam:
      try:
        fid = int(field)
      except:
        printError('ERROR with field %s'%field)
  else:
      fid = snam.index(field)




# Load antenna info:
  tb.open(os.path.join(vis,'ANTENNA'))
  anam = list(tb.getcol('NAME'))
  nant = len(anam)
  tb.close()



# Sanity checks:
  try:
      DR = list(DR)
      DL = list(DL)
      DRSolve = list(DRSolve)
      DLSolve = list(DLSolve)
      mounts = list(mounts)

  except:
      printError('ERROR: DR and DL must be lists of complex numbers!\n DRSolve and DLSolve must be lists of booleans.\n mounts must be a list of strings.')


  if len(DR) not in [0,nant]:
      printError('ERROR: DR should have %i elements!'%nant)
  else:
     for it in range(len(DR)):
          DR[it] = np.complex128(DR[it])
             
  if len(DL) not in [0,nant]:
      printError('ERROR: DL should have %i elements!'%nant)
  else:
     for it in range(len(DL)):
          DL[it] = np.complex128(DL[it])

  if len(DRSolve) not in [0,nant]:
      printError('ERROR: DRSolve should have %i elements!'%nant)
  else:
     for it in range(len(DRSolve)):
          DRSolve[it] = bool(DRSolve[it])
             
  if len(DLSolve) not in [0,nant]:
      printError('ERROR: DLSolve should have %i elements!'%nant)
  else:
     for it in range(len(DLSolve)):
          DLSolve[it] = bool(DLSolve[it])


  if len(mounts) not in [0,nant]:
      printError('ERROR: mounts should have %i elements!'%nant)
  else:
     for it in range(len(mounts)):
       mounts[it] = str(mounts[it])
       if mounts[it] not in GoodMounts:   
         printError('ERROR: unknown mount %s!'%mounts[it])
         
# Default values for Dterms and mounts:
  if len(DR)==0:
      DR = [0.+1.j*0. for i in range(nant)]
      
  if len(DL)==0:
      DL = [0.+1.j*0. for i in range(nant)]
      
  if len(DRSolve)==0:
      DRSolve = [True for i in range(nant)]
      
  if len(DLSolve)==0:
      DLSolve = [True for i in range(nant)]

  if len(mounts)==0:
      mounts = ['AZ' for i in range(nant)]



# Feed angles:

  FeedAngles = np.zeros(nant)

  try:
    feed_rotation = list(feed_rotation)

  except:
    printError('ERROR: feed_rotation should be a list of floats!')
   
   
  if len(feed_rotation) not in [0,nant]:
    printError('ERROR: feed_rotation should have %i elements!'%nant)

  for it in range(len(feed_rotation)):
    try:  
      FeedAngles[it] = np.pi/180.*float(feed_rotation[it]) 
    except:
      printError('ERROR: bad feed rotation for antenna #%i!'%it)




    

# Print info for sanity checking:
  printMsg('\nThere are %i antennas'%nant)
  F = {True: 'FITTABLE',False: 'FIXED   '}
  for ai in range(nant):
    printMsg('Ant %i (%s). Mount %s. Leakage: DR (%.3f, %.3f, %s); DL (%.3f, %.3f, %s)'%(ai,anam[ai],mounts[ai],DR[ai].real, DR[ai].imag,F[DRSolve[ai]], DL[ai].real, DL[ai].imag, F[DLSolve[ai]]))


# Get polarization info from MS:
  ms.open(vis)
  polprods = [dd[0] for dd in ms.range(['corr_names'])['corr_names']]
  spwFreqs = ms.range('chan_freq')['chan_freq'][:,0]


# More sanity checks:
  for pi in polprods: # By now, only circular feeds are allowed.
    if pi not in ['RR','RL','LR','LL']:
      ms.close()  
      printError("ERROR! Wrong pol. product %s"%pi)

  if len(polprods) != 4:
    printError('Measurement set is not a valid full-polarization dataset!')


# Load data:
  ms.selectinit(datadescid=int(spw))
  if not ms.select({'field_id':fid}):
    ms.close()
    printError('Problem reading data for source #%i'%fid)

  DATA = ms.getdata(['u','v','w','antenna1','antenna2','data','weight','flag'])
  ms.close()

# Metadata and flags:
  Nchan,Nvis = np.shape(DATA['data'])[1:]

  GoodData = np.logical_not(DATA['flag']) 
  GoodWgt =  np.copy(DATA['weight']) 


#################
# Lines for testing (same weight for all corr products)
#  GoodData[:] = np.prod(np.logical_not(DATA['flag']),axis=0)[np.newaxis,:]
#  GoodWgt[:] = np.sum(DATA['weight'],axis=0)[np.newaxis,:]
#################


# Weight (set to zero for bad data and autocorrs):
  Wgt = [np.zeros((Nvis,Nchan)) for i in range(4)]
  for i in range(4):
    j = polprods.index(POL_ORDER[i])  
    Wgt[i][:,:] = GoodWgt[j,:,np.newaxis]
    Wgt[i][np.logical_not(np.transpose(GoodData[j,...]))] = 0.0
    Wgt[i][DATA['antenna1']==DATA['antenna2'],:] = 0.0

# Correct weight power:
    Bads = Wgt[i] == 0.0
    Wgt[i][:] = np.power(Wgt[i],wgt_power)
    Wgt[i][Bads] = 0.0



# Arrange data optimally (i.e., [time,channel,polariz]):
  DATAPol = np.zeros((Nvis,Nchan,4),dtype=np.complex128)
  for i in range(4):
    DATAPol[:,:,i] = np.transpose(DATA['data'][polprods.index(POL_ORDER[i]),:,:])

# Free some memory:
  del GoodData, GoodWgt, DATA['data'], DATA['weight'], DATA['flag']





######################################
# Computing parallactic angles:

  printMsg('\nComputing parallactic (feed) angles')


  if plot_parang:

    NMAX = 8  # Maximum number of antennas in the plot
    col = ['r','g','b','c','m','y','k','w'] # colors

    PAs, TT, A1, A2, FG = getParangle(vis,spw,fid,mounts,FeedAngles,doplot=True)

    MIN = PAs[:,0] < -np.pi
    PAs[MIN,0] += 2.*np.pi
    MIN = PAs[:,1] < -np.pi
    PAs[MIN,1] += 2.*np.pi

    MAX = PAs[:,0] > np.pi
    PAs[MAX,0] -= 2.*np.pi
    MAX = PAs[:,1] > np.pi
    PAs[MAX,1] -= 2.*np.pi


    UT = TT/86400.
    UT -= int(np.min(UT))
    UT *= 24.

  
    tb.open(vis+'/ANTENNA')
    NAMES = tb.getcol('NAME')
    tb.close()


    OFF = open('polsolve_feed-angles.dat','w')


# Subsets of plots of NMAX antennas:
    SUBSET = []
    j = 0
    i = 0
    while i< len(NAMES):
      if i+NMAX <= len(NAMES):
        SUBSET.append(NAMES[i:i+NMAX])
        i += NMAX
        j += 1
      else:
        SUBSET.append(NAMES[i:])
        j += 1
        break


    for k in range(j):
      fig = pl.figure()
      sub = fig.add_subplot(111)



      for i in range(len(SUBSET[k])):

        sub.set_ylabel('Feed angle (deg.)')

        if i + NMAX*k >0:  

          mask = (A2==i+NMAX*k)*(A1!=i+NMAX*k)*FG
          TOPLOT = np.copy(-PAs[mask,1]*180./np.pi)
          PLOTTED = np.sum(mask)>0

          #for mi,m in enumerate(mask):
           # if m:  
           #   print >> OFF,NAMES[i+NMAX*k],UT[mi], -PAs[mi,1]*180./np.pi

          sub.plot( UT[mask], TOPLOT, 'o%s'%col[i], label=NAMES[i+NMAX*k])
         
          mask = (A1==i+NMAX*k)*(A2!=i+NMAX*k)*FG
          TOPLOT = np.copy(-PAs[mask,0]*180./np.pi)

          sub.plot( UT[mask], TOPLOT, 'o%s'%col[i])

          #for mi,m in enumerate(mask):
           # if m:  
           #   print >> OFF,NAMES[i+NMAX*k],UT[mi], -PAs[mi,0]*180./np.pi


        else:
          mask = (A1==0)*(A2!=1)*(A2!=0)*FG  #+(A2==3))
          TOPLOT = np.copy(-PAs[mask,0]*180./np.pi)

          sub.plot( UT[mask], TOPLOT, 'o%s'%col[i], label=NAMES[i+NMAX*k])


          #for mi,m in enumerate(mask):
          #  if m:  
          #    print >> OFF,NAMES[i+NMAX*k],UT[mi], -PAs[mi,0]*180./np.pi
 
        pl.legend(numpoints=1)
        
      sub.set_xlim((np.min(UT),np.max(UT) + 0.3*(np.max(UT)-np.min(UT))))
      sub.set_ylim((-189.,189.))
      sub.set_xlabel('UT (h)')  

      pl.savefig('%s_Parang_plot_%i.png'%(vis,k))
     
    OFF.close()
    
  #  raw_input('Press ENTER to continue...')
  
  else:
 
    PAs = getParangle(vis,spw,fid,mounts,FeedAngles)
  


#######################################

  gc.collect()
  


################
# Compute Fourier transform of each source component:
  COMPS = np.zeros((Nvis,Nchan,NmodComp+1),dtype=np.complex128)

  FouFac = 1.j*2.*np.pi*(np.pi/180.)
  U = DATA['u'][:,np.newaxis]/(3.e8/spwFreqs[np.newaxis,:])
  V = DATA['v'][:,np.newaxis]/(3.e8/spwFreqs[np.newaxis,:])

  for i in range(NmodComp):
    printMsg( 'Computing source component %i of %i'%(i+1,NmodComp))
    gc.collect()
    
    for delt in DELTAS[i]:
      COMPS[:,:,i] += delt[0]*np.exp(FouFac*(delt[1]*U + delt[2]*V)) 
   #   print U[0], V[0], COMPS[100,0,i],DATAPol[100,0,0],EMA[100],EPA[100]

    COMPS[:,:,NmodComp] += COMPS[:,:,i]

# How much flux do we have in the DATA?
  DSUM = 0.5*np.sum(np.abs(np.average((Wgt[0]+Wgt[1])*(DATAPol[:,:,polprods.index('RR')] + DATAPol[:,:,polprods.index('LL')]),axis=1)))/(np.sum(np.average(Wgt[0]+Wgt[1],axis=1)))


#  print np.max(Wgt), DSUM

# How much flux do we have in the MODEL?
  MSUM = np.sum(np.abs(np.average((Wgt[0]+Wgt[1])*COMPS[:,:,NmodComp],axis=1)))/(np.sum(np.average((Wgt[0]+Wgt[1]),axis=1)))
  
#  print np.max(Wgt), MSUM
  
# Notice that the Chi2 will take out the effects of this FluxFactor!  
  FluxFactor = MSUM/DSUM
  printMsg('All model components describe up to %.2f percent of the signal\n'%(FluxFactor*100.))
  
 


####################
##### CODE FOR TESTING (WRITE MODEL INTO MS):
#  ms.open(vis,nomodify=False)
#  ms.selectinit(datadescid=int(spw))
#  ms.select({'field_id':fid})
#  MDATA = ms.getdata(['data','corrected_data','model_data'])
#  MDATA['model_data'][0,0,:] = COMPS[:,0,NmodComp]
#  MDATA['model_data'][1,0,:] = COMPS[:,0,NmodComp]
#  MDATA['model_data'][2,0,:] = COMPS[:,0,NmodComp]
#  MDATA['model_data'][3,0,:] = COMPS[:,0,NmodComp]
#
#  MDATA['corrected_data'][:,:,:] = MDATA['data'] - MDATA['model_data']
#  ms.putdata(MDATA)
#  ms.close()
#  del MDATA['model_data'],MDATA['data'],MDATA['corrected_data']
#  del MDATA
####################




# Sum and difference of PANGS (in complex form):
  EPA = np.exp(1.j*(PAs[:,0]+PAs[:,1])) ; EMA = np.exp(1.j*(PAs[:,0]-PAs[:,1]))

# Get data back into antenna frame (if it was in the sky frame):
  if parang_corrected:
    printMsg("Undoing parang correction")
    for i in range(4):
      if POL_ORDER[i][0]=='R':
        DATAPol[:,:,i] *= np.exp(1.j*(PAs[:,0]))[:,np.newaxis]
      else:
        DATAPol[:,:,i] /= np.exp(1.j*(PAs[:,0]))[:,np.newaxis]

      if POL_ORDER[i][1]=='R':
        DATAPol[:,:,i] /= np.exp(1.j*(PAs[:,1]))[:,np.newaxis]
      else:
        DATAPol[:,:,i] *= np.exp(1.j*(PAs[:,1]))[:,np.newaxis]












# Code for testing (turned off):
  if False:
    import pickle as pk
    OFF = open('FTMODEL.dat','w')
    pk.dump([U,V,COMPS,DATAPol,Wgt],OFF)
    OFF.close()

################


# Number of fittable parameters (i.e., sources + antennas):   
  Npar = 2*np.sum(DRSolve) + 2*np.sum(DLSolve) + 2*np.sum(FITSOU)

# Number of model variables (i.e., fittable + fixed parameters):
  Nvar = 4*nant + 2*NmodComp
  VarSou = np.sum(FITSOU) # Number of variables related to source components.

# Vectors to store parameter indices in the covariance matrix:
  PAntL = np.zeros(nant,dtype=np.int32)
  PAntR = np.zeros(nant,dtype=np.int32)
  PSou = np.zeros(NmodComp,dtype=np.int32)

# All model variables (both fittable and fixed):
  VAntL = np.zeros(nant,dtype=np.int32)
  VAntR = np.zeros(nant,dtype=np.int32)
  VSou = np.zeros(NmodComp,dtype=np.int32)

# Figure out the indices of each parameter of the cov. matrix:
# First elements are the source Q and U.
# Then, DRs (Re and Im). Then, DLs (Re and Im):
  k = 0
  for i in range(NmodComp):
    if FITSOU[i]:
      PSou[i] = k
      k += 2
    else:
      PSou[i] = -1

  for i in range(nant):
    if DRSolve[i]:
      PAntR[i] = k
      k += 2
    else:
      PAntR[i] = -1

  for i in range(nant):
    if DLSolve[i]:
      PAntL[i] = k
      k += 2
    else:
      PAntL[i] = -1


# Initial parameter values:
  InitVal = np.zeros(Nvar)
  i = 0
  for fi in range(NmodComp):
    InitVal[i] = Pfrac[fi]*np.cos(EVPA[fi]*np.pi/90.)
    InitVal[i+1] = Pfrac[fi]*np.sin(EVPA[fi]*np.pi/90.)
    VSou[fi] = i
    i += 2
  for ni in range(nant):
    InitVal[i] = DR[ni].real
    InitVal[i+1] = DR[ni].imag
    VAntR[ni] = i
    i += 2
  for ni in range(nant):
    InitVal[i] = DL[ni].real
    InitVal[i+1] = DL[ni].imag
    VAntL[ni] = i
    i += 2


# Model variables (fittable and fixed):
  FitVal = np.copy(InitVal)

# Memory location for Hessian and residuals vector:
  Hessian = np.zeros((Npar,Npar),dtype=np.float)
  ResVec = np.zeros(Npar,dtype=np.float)

# Set memory and send data to C++:
  success = PS.setData(DATAPol,Wgt,DATA['antenna1'],DATA['antenna2'],COMPS,EPA,EMA,PSou,PAntR,PAntL,VSou,VAntR,VAntL,FitVal,Hessian,ResVec,FluxFactor,doCirc)

  if success != 0:
    printError('ERROR READING DATA INTO C++: CODE %i'%success)



# Just code for testing (turned off):
  if False:
    myfit = minimize(getChi2,[0. for i in range(Npar)], args= (NmodComp, nant, FITSOU, DRSolve, DLSolve), method='nelder-mead')
    iff = open('polsolve.dterms','w')
    import pickle as pk
    pk.dump(myfit,iff)
    iff.close()
    print myfit








##################################
################################################
##############################################################
####################################################################
# Levenberg-Marquardt in-house implementation:


  def LMFit(pini):

    NITER = 3*Npar; Lambda = 1.e-6; kfac = 5.0; functol = 1.e-9


# For testing:
#    NITER = 1
#    Lambda = 0.0
##############

    Chi2 = 0

    p = np.array(pini)
    backupP = np.copy(p)

    HessianDiag = np.zeros(np.shape(Hessian),dtype=np.float64)
    backupHess = np.zeros(np.shape(Hessian),dtype=np.float64)
    backupGrad = np.zeros(np.shape(ResVec),dtype=np.float64)
    Inverse = np.zeros(np.shape(Hessian),dtype=np.float64)

    Hessian[:,:] = 0.0
    ResVec[:] = 0.0

    CurrChi = getChi2(pini, NmodComp, nant, FITSOU, DRSolve, DLSolve)
    backupHess[:,:] = Hessian
    backupGrad[:] = ResVec

    controlIter = 0

    for i in range(NITER):
      controlIter += 1
      for n in range(len(p)):
        HessianDiag[n,n] = Hessian[n,n]
      try:
        goodsol = True
        Inverse[:] = np.linalg.pinv(Hessian+Lambda*HessianDiag)
        Dpar = np.dot(Inverse,ResVec)
        DirDer = sum([Hessian[n,n]*Dpar[n]*Dpar[n] for n in range(len(p))])
        DirDer2 = np.sum(ResVec*Dpar)
        TheorImpr = DirDer-2.*DirDer2 

      except:
        goodsol=False
        Dpar = 0.0
        TheorImpr = -10.0
      p += Dpar
      if controlIter==NITER:
        break
      if goodsol:
        Hessian[:,:] = 0.0
        ResVec[:] = 0.0
        Chi2 = getChi2(p, NmodComp, nant, FITSOU, DRSolve, DLSolve)
        RealImpr = Chi2 - CurrChi
      else:
        RealImpr = 1.0

      if TheorImpr != 0.0:
        Ratio = RealImpr/TheorImpr
      else:
        Ratio = 0.0


      if Ratio < 0.25:

        if RealImpr<0.0:
           temp = np.sqrt(kfac)
        else:
           temp = kfac

      elif Ratio > 0.75:
        temp = 1./np.sqrt(kfac)

      elif not goodsol:
        temp = kfac

      else:
        temp = 1.0

      Lambda *= temp
      if Chi2 == 0.0 or CurrChi==0.0:
        break

      relchi = Chi2/CurrChi
      if relchi<1: relchi = 1./relchi 
      todivide = np.copy(backupP); todivide[todivide==0.0] = 1.0
      relpar = max([{True:1./pb,False:pb}[pb<1] for pb in p/todivide])

      if relchi-1.<functol and relpar-1.<functol: 
        p[:] = backupP
        break

      if CurrChi<Chi2:
        p[:] = backupP
        Hessian[:,:] = backupHess
        ResVec[:] = backupGrad
      else:
        CurrChi = Chi2
        backupHess[:,:] = Hessian
        backupGrad[:] = ResVec
        backupP[:] = p

 #   if controlIter == NITER:
 #     sys.stdout.write("\n\n REACHED MAXIMUM NUMBER OF ITERATIONS!\nThe algorithm may not have converged!\nPlease, check if the parameter values are meaningful.\n")
 #     sys.stdout.flush()


    try:
      return [p[:], np.linalg.pinv(Hessian), Chi2, getChi2(p, NmodComp, nant, FITSOU, DRSolve, DLSolve, arg=1)]
    except:
      return False


####################################################################
##############################################################
################################################
##################################








# Prepare vector of initial parameters:
# First elements are Q and U of the fittable source components.
# Then, the fittable Dterms (Re and Im) for R.
# Then, the fittable Dterms (Re and Im) for L.
  i = 0
  pini = np.zeros(Npar)
  for si in range(NmodComp):
    if FITSOU[si]:
      Q = Pfrac[si]*np.cos(EVPA[si]*np.pi/90.)
      U = Pfrac[si]*np.sin(EVPA[si]*np.pi/90.)
      pini[i] = Q ; pini[i+1] = U
      i += 2
  for si in range(nant):
    if DRSolve[si]:
      pini[i] = DR[si].real ; pini[i+1] = DR[si].imag
      i += 2

  for si in range(nant):
    if DLSolve[si]:
      pini[i] = DR[si].real ; pini[i+1] = DR[si].imag
      i += 2

# FIT!
  LM = LMFit(pini)

  if LM == False:
    printError('AN ERROR OCCURRED DURING THE MINIMIZATION PROCESS!')

  Errors = [np.sqrt(np.abs(LM[1][i,i])*LM[2]) for i in range(Npar)]

  ErrR = []
  ErrL = []


  ifile = open('%s.spw%i.PolSolve.CovMatrix'%(vis,int(spw)),'w')
  print >> ifile,'!  PARAMETERS:'

# Print final estimates of source polarization and Dterms:
  i = 0
  printMsg('\nFitting results:')
  for si in range(NmodComp):
    if FITSOU[si]:
      Q = LM[0][i]
      U = LM[0][i+1]
      ErrQ = Errors[i] ; ErrU = Errors[i+1]
      print >> ifile,' %i -> Component %i; Stokes Q'%(i,si)
      print >> ifile,' %i -> Component %i; Stokes U'%(i+1,si)
      i += 2
      # Unbiased estimates and errors (quick&dirty MonteCarlo):
      Nmc = 1000
      Qmc = np.random.normal(Q,ErrQ,Nmc) ; Umc = np.random.normal(U,ErrU,Nmc)
      Pmc = np.sqrt(Qmc*Qmc + Umc*Umc) ; EVmc = np.arctan2(Umc,Qmc)*90./np.pi

      for mi in range(Nmc):
        dAng = EVmc[mi] - EVmc[0]
        if dAng > 180.:
          EVmc[mi] -= 180.
        elif dAng < -180.:
          EVmc[mi] += 180.
          

      Pr = np.average(Pmc)
      ErrPr = np.std(Pmc)

      EV = np.average(EVmc)
      ErrEV = np.std(EVmc)

    else:
      Pr = Pfrac[si]
      EV = EVPA[si]
      ErrPr = 0.0 
      ErrEV = 0.0
    printMsg('Source component #%i: Pfrac = % 6.4f +/- %.3f; EVPA = % 6.2f +/- %.2f deg.'%(si,Pr,ErrPr,EV,ErrEV))

  printMsg('\n Dterms (Right):')
  for si in range(nant):
    if DRSolve[si]:
      Re = LM[0][i] ; Im = LM[0][i+1]
      ErrRe = Errors[i]; ErrIm = Errors[i+1]
      DR[si] = Re + 1.j*Im
      print >> ifile,' %i -> Antenna %i; Dterm R (real)'%(i,si)
      print >> ifile,' %i -> Antenna %i; Dterm R (imag)'%(i+1,si)
      i += 2
    else:
      Re = DR[si].real ; Im = DR[si].imag
      ErrRe = 0.0 ; ErrIm = 0.0
    ErrR.append(np.sqrt(ErrRe*ErrIm))
    printMsg('Antenna #%i (%s):  Real = % .2e +/- %.1e; Imag = % .2e +/- %.1e | Amp = %.4f | Phase: %.2f deg.'%(si,anam[si], Re,ErrRe,Im,ErrIm, np.sqrt(Re*Re+Im*Im),np.arctan2(Im,Re)*180./np.pi))

  printMsg('\n Dterms (Left):')
  for si in range(nant):
    if DLSolve[si]:
      Re = LM[0][i] ; Im = LM[0][i+1]
      ErrRe = Errors[i]; ErrIm = Errors[i+1]
      DL[si] = Re + 1.j*Im
      print >> ifile,' %i -> Antenna %i; Dterm L (real)'%(i,si)
      print >> ifile,' %i -> Antenna %i; Dterm L (imag)'%(i+1,si)      
      i += 2
    else:
      Re = DL[si].real ; Im = DL[si].imag
      ErrRe = 0.0 ; ErrIm = 0.0
    ErrL.append(np.sqrt(ErrRe*ErrIm))
    printMsg('Antenna #%i (%s):  Real = % .2e +/- %.1e; Imag = % .2e +/- %.1e | Amp = %.4f | Phase: %.2f deg.'%(si,anam[si], Re,ErrRe,Im,ErrIm, np.sqrt(Re*Re+Im*Im),np.arctan2(Im,Re)*180./np.pi))


# Write CovMatrix:
  print >> ifile,'!  POST-FIT COVARIANCE MATRIX:'

  fmt = '%.2e '*Npar
  for i in range(Npar):
    print >> ifile,fmt%tuple(LM[1][i,:])

  ifile.close()
  
  

# Write calibration table:
  form = '%.1f  % i  %i  %i  % i  %.1f  % i  % i  % .8f  % .8f  % .8f  % .8f  % .8f  % .8f  %i  %i  %.3f  %.3f\n'
  if os.path.exists('%s.spw%i.Dterms'%(vis,int(spw))):
    os.system('rm -rf %s.spw%i.Dterms'%(vis,int(spw)))
  ascf = open('%s.Dterms.ascii'%vis,'w')
  ascf.write('TIME FIELD_ID SPECTRAL_WINDOW_ID ANTENNA1 ANTENNA2 INTERVAL SCAN_NUMBER OBSERVATION_ID CPARAM PARAMERR FLAG SNR\n')
  ascf.write('D I I I I D I I X2,1 R2,1 B2,1 R2,1\n')


  for j in range(nant):
      ascf.write(form%(0,-1,int(spw),j,-1,0,-1,-1, DR[j].real, DR[j].imag, DL[j].real, DL[j].imag,ErrR[j], ErrL[j],0 ,0, 100., 100.))


  ascf.close()
  tb.fromascii(tablename='%s.spw%i.Dterms'%(vis,int(spw)),asciifile='%s.Dterms.ascii'%vis,sep=' ')
  tb.close()
  tb.open('%s.spw%i.Dterms'%(vis,int(spw)),nomodify=False)
  tb.putinfo({'readme': '', 'subType': 'D Jones', 'type': 'Calibration'})
  tb.close()






  targets = str(target_field)
  if len(targets)>0:
    selFields = getFields(vis,targets)
  #  if not parang_corrected:
  #    printError('IN THE CURRENT VERSION, DTERMS ARE ONLY APPLIED IF parang_corrected is True!')
  else:
    selFields = []
  
  
  DRa = -np.array(DR,dtype=np.complex128)
  DLa = -np.array(DL,dtype=np.complex128)
  SFac = 1./(1.-DRa*DLa)    
      
  for target in selFields:
    printMsg('\nApplying calibration to field id %i'%target)  
    ms.open(vis)
    ms.selectinit(datadescid=int(spw))
    scans = ms.range('scan_number')['scan_number']
    ms.close()
    for sci in scans:
      sys.stdout.write('\r  SCAN %i'%sci)
      sys.stdout.flush()
      
      PAs = getParangle(vis,spw,target,mounts, FeedAngles, scan=sci)
        
      ms.open(vis,nomodify=False)
      ms.selectinit(datadescid=spw)
      ms.select({'scan_number':sci})
      try:
        DATA = ms.getdata(['data','corrected_data','antenna1','antenna2'])
      except:
        ms.close()
        printError('COULD NOT READ CORRECTED DATA! DID YOU RUN CLEARCAL??')
        


# Sum and difference of PANGS (in complex form):
      EPA = np.exp(1.j*(PAs[:,0]+PAs[:,1])) ; EMA = np.exp(1.j*(PAs[:,0]-PAs[:,1]))


# Get data back into antenna frame:
      if parang_corrected:
        DATA['corrected_data'][polprods.index('RR'),:,:] = DATA['data'][polprods.index('RR'),:,:]*EMA[np.newaxis,:]
        DATA['corrected_data'][polprods.index('RL'),:,:] = DATA['data'][polprods.index('RL'),:,:]*EPA[np.newaxis,:]
        DATA['corrected_data'][polprods.index('LR'),:,:] = DATA['data'][polprods.index('LR'),:,:]/EPA[np.newaxis,:]
        DATA['corrected_data'][polprods.index('LL'),:,:] = DATA['data'][polprods.index('LL'),:,:]/EMA[np.newaxis,:]

# Apply Dterms:
      SDt = SFac[DATA['antenna1']]*np.conjugate(SFac[DATA['antenna2']])
      BKP_RR = np.copy(DATA['corrected_data'][polprods.index('RR'),:,:])*SDt
      BKP_RL = np.copy(DATA['corrected_data'][polprods.index('RL'),:,:])*SDt
      BKP_LR = np.copy(DATA['corrected_data'][polprods.index('LR'),:,:])*SDt
      BKP_LL = np.copy(DATA['corrected_data'][polprods.index('LL'),:,:])*SDt

      DATA['corrected_data'][polprods.index('RR'),:,:] = BKP_RR + np.conjugate(DRa[DATA['antenna2']])*BKP_RL + DRa[DATA['antenna1']]*BKP_LR + DRa[DATA['antenna1']]*np.conjugate(DRa[DATA['antenna2']])*BKP_LL
      DATA['corrected_data'][polprods.index('RL'),:,:] = BKP_RL + DRa[DATA['antenna1']]*np.conjugate(DLa[DATA['antenna2']])*BKP_LR + DRa[DATA['antenna1']]*BKP_LL + np.conjugate(DLa[DATA['antenna2']])*BKP_RR
      DATA['corrected_data'][polprods.index('LR'),:,:] = BKP_LR + DLa[DATA['antenna1']]*np.conjugate(DRa[DATA['antenna2']])*BKP_RL + DLa[DATA['antenna1']]*BKP_RR + np.conjugate(DRa[DATA['antenna2']])*BKP_LL
      DATA['corrected_data'][polprods.index('LL'),:,:] = BKP_LL + np.conjugate(DLa[DATA['antenna2']])*BKP_LR + DLa[DATA['antenna1']]*BKP_RL + DLa[DATA['antenna1']]*np.conjugate(DLa[DATA['antenna2']])*BKP_RR

# Put data into sky frame:
      DATA['corrected_data'][polprods.index('RR'),:,:] /= EMA[np.newaxis,:]
      DATA['corrected_data'][polprods.index('RL'),:,:] /= EPA[np.newaxis,:]
      DATA['corrected_data'][polprods.index('LR'),:,:] *= EPA[np.newaxis,:]
      DATA['corrected_data'][polprods.index('LL'),:,:] *= EMA[np.newaxis,:]

# Save data:
      ms.putdata(DATA)
      ms.close()
   

# Release memory:
      del BKP_RR, BKP_LL, BKP_RL, BKP_LR, DATA['corrected_data'], DATA['antenna1'], DATA['antenna2'], SDt
      gc.collect()
    

if __name__=='__main__':

  polsolve(vis, spw, field, mounts, DR, DL, DRSolve, DLSolve, 
           CLEAN_models, Pfrac, EVPA, PolSolve, parang_corrected,target_field)


