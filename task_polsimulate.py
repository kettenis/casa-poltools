# PolSimulate - A task to simulate simple full-polarization ALMA data.
#
# Copyright (c) Ivan Marti-Vidal - Nordic ARC Node (2016). 
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

import gc
from simutil import *
import os
import numpy as np
import scipy.interpolate as spint
from taskinit import gentools
from taskinit import casalog
from clearcal_cli import clearcal_cli as clearcal
from exportuvfits_cli import exportuvfits_cli as exportuvfits
from importuvfits_cli import importuvfits_cli as importuvfits

from ft_cli import ft_cli as ft
from simutil import *
ms = gentools(['ms'])[0]
sm = gentools(['sm'])[0]
me = gentools(['me'])[0]
tb = gentools(['tb'])[0]

__version__ = '1.3.2'

#####################
# UNIT TEST LINES:
if __name__=='__main__':

  vis                =  "EHT_TEST6"
  array_configuration    =  "EHT.cfg"
  feed               =  "circular"
  mounts             =  ['AZ', 'NR', 'NR', 'AZ', 'NL', 'NL', 'AZ', 'NL']
  ConstDt0           =  [(0.05+0.05j), (0.05+0.05j), (0.05+0.05j), (0.05+0.05j), (0.05+0.05j), (0.05+0.05j), (0.05+0.05j), (0.05+0.05j)]
  ConstDt1           =  [(0.05+0.05j), (0.05+0.05j), (0.05+0.05j), (0.05+0.05j), (0.05+0.05j), (0.05+0.05j), (0.05+0.05j), (0.05+0.05j)]
  LO                 =  2.3e+11
  BBs                =  [0.0]
  spw_width          =  2000000000.0
  nchan              =  4
  model_image        =  []
  I                  =  [1.4]
  Q                  =  [0.2]
  U                  =  [0.2]
  V                  =  [0.2]
  RM                 =  [0.0]
  spec_index         =  [0.0]
  RAoffset           =  [0.0]
  Decoffset          =  [0.0]
  spectrum_file      =  ""
  incenter           =  "J2000 12h30m49.4234 12d23m28.0439"
  incell             =  ""
  inbright           =  ""
  inwidth            =  ""
  H0                 =  "2017/04/01 00:00:00"
  onsource_time      =  4
  observe_time       =  5
  visib_time         =  "2s"
  nscan              =  [0., 0.2, 0.4, 0.6, 0.7]
  corrupt            =  False
  seed               =  42
  Dt_amp             =  0.0
  Dt_noise           =  0.0
  tau0               =  0.0
  t_sky              =  250.0
  t_ground           =  270.0
  t_receiver         =  50.0
  export_uvf = True
#
#
#
#
##################






def polsimulate(vis = 'polsimulate_output.ms', array_configuration='alma.out04.cfg', feed = 'linear',
                mounts = [], ConstDt0 = [], ConstDt1 = [], LO=200.e9, BBs = [-7.e9,-5.e9,5.e9,7.e9], spw_width = 2.e9, nchan = 128, 
                model_image=[], I = [], Q = [], U = [], V = [], RM = [], spec_index = [], RAoffset = [], Decoffset = [], spectrum_file = '',
                incenter = 'J2000 00h00m00.00 -00d00m00.00', incell='',inbright='', 
                inwidth='', H0 = -1.5, 
                onsource_time=1.5, observe_time = 3.0,visib_time='6s',nscan = 50, apply_parang = False, export_uvf=True,
                corrupt = True, seed=42, 
                Dt_amp = 0.00, Dt_noise = 0.001, tau0=0.0, t_sky=250.0, t_ground=270.0, t_receiver=50.0):



########
# Uncomment for testing/debugging
#  return True
#if __name__=='__main__':
########



  def printError(msg):
    print '\n', msg, '\n'
    casalog.post('PolSimulate: '+msg)
    raise Exception(msg)

  def printMsg(msg,dolog=True):
    print '\n',msg, '\n'
    if dolog:
      casalog.post('PolSimulate: '+msg)


  util = simutil('')


  printMsg( 'POLSIMULATE - VERSION %s  - Nordic ARC Node'%__version__)


  array = os.path.basename(array_configuration).split('.')[0].upper()


# ALMA bands:
  Bands = {'3':[84,119],'5':[163,211],'6':[211,275],'7':[275,370],'8':[385,500],'9':[602,720],'10':[787,950]}

# Receiver evector angles. Not used by now:
  Pangs = {'3':0, '5':0, '6':0, '7':0, '8':0, '9':0, '10':0}


  if array=='ALMA':
    found = False
    for band in Bands.keys():
      if LO/1.e9 > Bands[band][0] and LO/1.e9 < Bands[band][1]:
        found = True
        selb = band
        break

    if not found:
      printError("Frequency %.5fGHz does NOT correspond to any ALMA band!"%(LO/1.e9))
    else:
      printMsg( 'This is a Band %s ALMA observation.\n'%selb)

  if feed in ['linear','circular']:
    printMsg('Will simulate feeds in %s polarization basis'%feed)
  else:
    printError('Unknown feed %s'%feed)


# Load the different models:

# Point source(s):
  if len(set([len(I),len(Q),len(U),len(V),len(RM),len(spec_index),len(RAoffset),len(Decoffset)]))>1:
    printError("ERROR! I, Q, U, V, RM, RAoffset, Decoffset and spec_index should all have the same length!")

  NmodComp = max([len(I),1])

# Point source (user-defined spectrum):
  ismodel = False
  if type(spectrum_file) is str and len(spectrum_file)>0:
    if not os.path.exists(spectrum_file):
      printError("ERROR! spectrum_file is not found!")
    else:
     try:
      ismodel = True
      iff = open(spectrum_file)
      lines = iff.readlines() ; iff.close()
      model = np.zeros((5,len(lines)))
      for li,line in enumerate(lines):
        temp = map(float,line.split())
        model[:,li] = temp[:5]
      model2 = model[:,np.argsort(model[0,:])]
      interpI = spint.interp1d(model2[0,:],model2[1,:])
      interpQ = spint.interp1d(model2[0,:],model2[2,:])
      interpU = spint.interp1d(model2[0,:],model2[3,:])
      interpV = spint.interp1d(model2[0,:],model2[4,:])
     except:
      printError("ERROR! spectrum_file has an incorrect format!")



# Extended source (and cube!):
  if type(model_image) is list:
    if len(model_image)>0:
      new_mod = [m + '.polsim' for m in model_image]
      for i in range(4):
        os.system('rm -rf %s'%new_mod[i])
        returnpars = util.modifymodel(model_image[i], new_mod[i],
                               inbright,incenter,incell,
                               incenter,inwidth,0,
                               flatimage=False)
      Iim, Qim, Uim, Vim = new_mod
    else:
      Iim = '';Qim = '';Uim = ''; Vim = '';
  else:
    printError("ERROR! Unkown model_image!")



  if len(Iim)>0 and not (os.path.exists(Iim) and os.path.exists(Qim) and os.path.exists(Uim) and os.path.exists(Vim)):
    printError("ERROR! one or more model_image components does not exist!")


  if len(model_image)==0 and len(I)==0 and not ismodel:
    printError("ERROR! No model specified!")
  

  if not os.path.exists(array_configuration):
    antlist = os.getenv("CASAPATH").split(' ')[0] + "/data/alma/simmos/"+array_configuration
  else:
    antlist = array_configuration

  try:
    stnx, stny, stnz, stnd, antnames, arrname, arrpos = util.readantenna(antlist)
    nant = len(antnames)
  except:
    printError('ERROR with array configuration file!')

  if type(mounts) is not list:
     printError('ERROR: mounts must be a list of strings')
  elif len(mounts)>0:
      if len(mounts) != nant:
        printError('ERROR: length of \'mounts\' should be %i'%nant)
      else:
        for mi,mnt in enumerate(mounts):
          if mnt not in ['AZ','XY','EQ','NR','NL']:
              printError('ERROR: Unknown mount type %s'%mnt)
          else:
              printMsg('\n   ANTENNA %s WILL HAVE A MOUNT %s'%(antnames[mi],mnt))
  else:
      mounts = ['AZ' for i in range(nant)]
      printMsg('ALL ANTENNAS ARE ASSUMED TO HAVE ALT-AZ MOUNTS')
      
              

# Setting noise
  if corrupt:
   eta_p, eta_s, eta_b, eta_t, eta_q, t_rx = util.noisetemp(telescope=array,freq='%.9fHz'%(LO))
   eta_a = eta_p * eta_s * eta_b * eta_t
   if t_receiver!=0.0:
    t_rx = abs(t_receiver)
   tau0 = abs(tau0)
   t_sky = abs(t_sky)
   t_ground = abs(t_ground)
  else:
   printMsg('Since \'corrupt\' is False, noise will NOT be added!')
   Dt_noise = 0.0
   Dt_amp = 0.0


  try:
      ConstDt0 = list(ConstDt0)
      ConstDt1 = list(ConstDt1)
  except:
      printError('ERROR: ConstDt0 and ConstDt1 must be lists of copmlex numbers!')


  if len(ConstDt0) not in [0,nant]:
      printError('ERROR: ConstDt0 should have %i elements!'%nant)
  else:
     for it in range(len(ConstDt0)):
          ConstDt0[it] = np.complex128(ConstDt0[it])
             
  if len(ConstDt1) not in [0,nant]:
      printError('ERROR: ConstDt1 should have %i elements!'%nant)
  else:
     for it in range(len(ConstDt1)):
          ConstDt1[it] = np.complex128(ConstDt1[it])

  if len(ConstDt0)==0:
      ConstDt0 = [0.+1.j*0. for i in range(nant)]
      
  if len(ConstDt1)==0:
      ConstDt1 = [0.+1.j*0. for i in range(nant)]
      
  os.system('rm -rf '+vis)
  sm.open(vis)


# Setting the observatory and the observation:

  ALMA = me.observatory(array)
  if len(ALMA.keys())==0:
    ALMA = me.observatory('VLBA')

  usehourangle = False

  mount = 'alt-az'
  if type(H0) is str:
      refdate = H0
      H0 = 0.0
  else:
    H0 = float(H0)
    refdate='2017/01/01/00:00:00'
    usehourangle = True


  integ = visib_time


  sm.setconfig(telescopename=array, x=stnx, y=stny, z=stnz,
             dishdiameter=stnd.tolist(),
             mount=mount, antname=antnames, padname=antnames, 
             coordsystem='global', referencelocation=ALMA)
  spwnames = ['spw%i'%i for i in range(len(BBs))]

  dNu = spw_width/nchan/1.e9
  spwFreqs = []
  dtermsX = [] ; dtermsY = []
  ModI = [np.zeros((nchan,NmodComp)) for i in BBs]
  ModQ = [np.zeros((nchan,NmodComp)) for i in BBs]
  ModU = [np.zeros((nchan,NmodComp)) for i in BBs]
  ModV = [np.zeros((nchan,NmodComp)) for i in BBs]


# Spectral windows and D-terms:

  corrp = {'linear':'XX XY YX YY','circular':'RR RL LR LL'}

  for i in range(len(BBs)):
    Nu0 = (LO+BBs[i]-spw_width/2.)/1.e9
    sm.setspwindow(spwname=spwnames[i], freq='%.8fGHz'%(Nu0),
               deltafreq='%.9fGHz'%(dNu),
               freqresolution='%.9fGHz'%(dNu), 
               nchannels=nchan, refcode="BARY",
               stokes=corrp[feed])
    spwFreqs.append(1.e9*np.linspace(Nu0,Nu0+dNu*nchan,nchan))
        
    if Dt_amp>0.0:
      DtX = [[np.random.normal(0.,Dt_amp), np.random.normal(0.,Dt_noise)] for j in stnx]
      DtY = [[np.random.normal(0.,Dt_amp), np.random.normal(0.,Dt_noise)] for j in stnx]
    else:
      DtX = [[0.,0.] for j in stnx]
      DtY = [[0.,0.] for j in stnx]

    dtermsX.append([np.zeros(nchan,dtype=np.complex128) for j in range(nant)])
    dtermsY.append([np.zeros(nchan,dtype=np.complex128) for j in range(nant)])

    for j in range(nant):
      dtermsX[-1][j][:] = ConstDt0[j] 
      dtermsY[-1][j][:] = ConstDt1[j]
      if Dt_amp > 0.0:
        dtermsX[-1][j] += np.random.normal(0.,Dt_amp)*np.exp(1.j*2.*np.pi*np.random.random())
        dtermsY[-1][j] += np.random.normal(0.,Dt_amp)*np.exp(1.j*2.*np.pi*np.random.random())

      if Dt_noise > 0.0:
        dtermsX[-1][j] += np.random.normal(0.,Dt_noise,nchan)*np.exp(1.j*2.*np.pi*np.random.random(nchan))
        dtermsY[-1][j] += np.random.normal(0.,Dt_noise,nchan)*np.exp(1.j*2.*np.pi*np.random.random(nchan))



# Compute point models:

    if len(I)>0:
      Lam2 = np.power(299792458./spwFreqs[i],2.)
      LamLO2 = (299792458./LO)**2.

      for j in range(len(I)):
        p = (Q[j]**2. + U[j]**2.)**0.5*np.power(spwFreqs[i]/LO,spec_index[j])
        phi0 = np.arctan2(U[j],Q[j])
        ModI[i][:,j] += I[j]*np.power(spwFreqs[i]/LO,spec_index[j])
        ModQ[i][:,j] += p*np.cos(2.*(RM[j]*(Lam2-LamLO2)) + phi0)
        ModU[i][:,j] += p*np.sin(2.*(RM[j]*(Lam2-LamLO2)) + phi0)
        ModV[i][:,j] += V[j] 
    if ismodel:
        ModI[i,0] += interpI(spwFreqs[i])
        ModQ[i,0] += interpQ(spwFreqs[i])
        ModU[i,0] += interpU(spwFreqs[i])
        ModV[i,0] += interpV(spwFreqs[i])


# CASA sm tool FAILS with X Y receiver. Will change it later:
  #  sm.setfeed(mode='perfect R L',pol=[''])
  #  sm.setauto(0.0)


# Field name:
  if len(model_image)>0:

    source = '.'.join(os.path.basename(model_image[0]).split('.')[:-1])
  else:
    source = 'POLSIM'

  sm.setfield(sourcename=source, sourcedirection=incenter,
          calcode="TARGET", distance='0m')


  mereftime = me.epoch('TAI', refdate)

  if usehourangle:
    printMsg( ' Will shift the date of observation to match the Hour Angle range\n')

  sm.settimes(integrationtime=visib_time, usehourangle=usehourangle, 
            referencetime=mereftime)


# Set scans:

  starttimes = []
  stoptimes = []
  sources = []

  try:
    scop = len(nscan)
    scandur = float(onsource_time)/scop
    nscan.sort()
    T0s = [H0 + (float(observe_time)-scandur)/nscan[-1]*sci for sci in nscan]
  except:    
    scandur = float(onsource_time)/nscan
    scop = nscan
    if nscan>1:
      T0s = [H0 + (float(observe_time)-scandur)/(nscan-1)*i for i in range(nscan)]
    else:
      TOs = [H0]


  for i in range(scop):
    sttime = T0s[i]
    endtime = (sttime + scandur)
    if i< scop-1 and endtime > T0s[i+1]:
      printMsg('WARNING! There are overlapping scans! Will shift them ot avoid collisions!')
      for j in range(i+1,scop):
         T0s[j] += (endtime - T0s[i+1]) + 1.0

    starttimes.append(str(3600.*sttime)+'s')
    stoptimes.append(str(3600.*endtime)+'s')
    sources.append(source)

  for n in range(scop):
   for sp in spwnames: 
    sm.observemany(sourcenames=[sources[n]],spwname=sp,starttimes=[starttimes[n]],stoptimes=[stoptimes[n]],project='polsimulate')
  
  sm.close()


# Change feeds to XY:

  if feed == 'linear':

    printMsg( ' CHANGING FEEDS TO X-Y\n')
    tb.open(vis+'/FEED',nomodify = False)
    pols = tb.getcol('POLARIZATION_TYPE')
    pols[0][:] = 'X'
    pols[1][:] = 'Y'
    tb.putcol('POLARIZATION_TYPE',pols)
    tb.close()











##############################
# Create an auxiliary MS:
  printMsg( 'Creating the auxiliary single-pol datasets')
  if feed == 'linear':
    polprods = ['XX','XY','YX','YY','I','Q','U','V']
  elif feed == 'circular':
    polprods = ['RR','RL','LR','LL','I','Q','U','V']

  dvis = [vis +ss for ss in polprods]

  for dv in dvis:
    os.system('rm -rf '+dv)

  sm.open(dvis[0])
  sm.setconfig(telescopename=array, x=stnx, y=stny, z=stnz,
             dishdiameter=stnd.tolist(),
             mount=mount, antname=antnames, padname=antnames, 
             coordsystem='global', referencelocation=ALMA)
  spwnames = ['spw%i'%i for i in range(len(BBs))]
  for i in range(len(BBs)):
    sm.setspwindow(spwname=spwnames[i], freq='%.8fGHz'%((LO+BBs[i]-spw_width/2.)/1.e9),
               deltafreq='%.9fGHz'%(spw_width/nchan/1.e9),
               freqresolution='%.9fGHz'%(spw_width/nchan/1.e9), 
               nchannels=nchan, refcode="BARY",
               stokes=polprods[0])


  sm.setfield(sourcename=source, sourcedirection=incenter,
          calcode="TARGET", distance='0m')

  sm.settimes(integrationtime=visib_time, usehourangle=usehourangle, 
            referencetime=mereftime)

  for n in range(scop):
   for sp in spwnames: 
    sm.observemany(sourcenames=[sources[n]],spwname=sp,starttimes=[starttimes[n]],stoptimes=[stoptimes[n]],project='polsimulate')

  sm.close()
##############################

  if feed == 'linear':
    printMsg( ' VISIBS WILL BE COMPUTED IN X-Y BASIS')



# Simulate Stokes parameters:

  clearcal(vis=dvis[0],addmodel=True)
  clearcal(vis=vis,addmodel=True)

  for dv in dvis[1:]:
    os.system('cp -r %s %s'%(dvis[0],dv))


# Auxiliary arrays:
  spwscans = []
  for n in range(len(spwnames)):
    ms.open(dvis[0])
    ms.selectinit(datadescid=n)
    spwscans.append(np.copy(ms.range('scan_number')['scan_number']))
    ms.close()

  ms.open(dvis[0])
  ms.selectinit(datadescid=0)
  ms.select({'scan_number':int(spwscans[0][0])})
  dataI = np.copy(ms.getdata(['data'])['data']) 
  dataQ = np.copy(dataI) 
  dataU = np.copy(dataI) 
  dataV = np.copy(dataI) 
  ms.close()


  ntimes = np.shape(dataI)[-1]

  if feed=='linear':
    printMsg( ' Simulating X-Y feed observations')
  else:
    printMsg( ' Simulating R-L feed observations')


  PAs = []; ant1 = [{} for i in range(len(BBs))] ; ant2 = [{} for i in range(len(BBs))]
  Flags = []






######################################
# Computing parallactic angles:

  printMsg('Computing parallactic angles')
  
  dirst = incenter.split()
  csys = cs.newcoordsys(direction=True)
  csys.setdirection(refcode=dirst[0], refval=' '.join(dirst[1:]))
  Dec = csys.torecord()['direction0']['crval'][1]

  CosDec = np.cos(Dec)
  SinDec = np.sin(Dec)

  tb.open(dvis[4]+'/ANTENNA')
  apos = tb.getcol('POSITION')
  Lat = np.arctan2(apos[2,:],np.sqrt(apos[0,:]**2. + apos[1,:]**2.))
  Tlat = np.tan(Lat)
  Lon = np.arctan2(apos[1,:],apos[0,:])
  tb.close()


  for i in range(len(BBs)):
   PAs.append([])
   for ni,n in enumerate(spwscans[i]):

    ms.open(dvis[4],nomodify=False)
    ms.selectinit(datadescid=i)
    ms.select({'scan_number':int(n)})
    temp = ms.getdata(['antenna1','antenna2','u','v','w'])
    temp2 = ms.getdata(['flag'])
    
    ant1[i][n] = temp['antenna1']
    ant2[i][n] = temp['antenna2']
    
    Ndata = np.shape(temp['u'])[0]
    PAs[i].append(np.zeros((Ndata,2)))
    
    V2 = SinDec*temp['v'] - CosDec*temp['w']
    
    
    Bx = -(apos[0,temp['antenna2']]-apos[0,temp['antenna1']])
    By = -(apos[1,temp['antenna2']]-apos[1,temp['antenna1']])
    Bz = -(apos[2,temp['antenna2']]-apos[2,temp['antenna1']])

    CH = temp['u']*By - V2*Bx
    SH = temp['u']*Bx + V2*By
    
    CT1 = CosDec*Tlat[temp['antenna1']]
    CT2 = CosDec*Tlat[temp['antenna2']]
    
    HAng = np.arctan2(SH,CH)
    
    H1 = HAng + Lon[temp['antenna1']]
    H2 = HAng + Lon[temp['antenna2']]
    
    
    Autos = (CH==0.)*(SH==0.)
    H1[Autos] = 0.0
    H2[Autos] = 0.0
    
    E1 = np.arcsin(SinDec*np.sin(Lat[temp['antenna1']])+np.cos(Lat[temp['antenna1']])*CosDec*np.cos(H1))
    E2 = np.arcsin(SinDec*np.sin(Lat[temp['antenna2']])+np.cos(Lat[temp['antenna2']])*CosDec*np.cos(H2))

    if i == 0:
      temp2['flag'][...,np.logical_or(E1<0.0,E1>np.pi)] = True
      temp2['flag'][...,np.logical_or(E2<0.0,E2>np.pi)] = True
      NscF = np.sum(temp2['flag'][0,0,:])
      if NscF>0:
        printMsg('#%i visibs will be flagged for negative elevations in scan #%i'%(NscF, n))
  
    for j in range(Ndata):
      if mounts[temp['antenna1'][j]] == 'AZ':
          PAs[i][ni][j,0] = -np.arctan2(np.sin(H1[j]), CT1[j] - SinDec*np.cos(H1[j]))
      elif mounts[temp['antenna1'][j]] == 'EQ':
          PAs[i][ni][j,0] = 0.0
      elif mounts[temp['antenna1'][j]] == 'XY':
          PAs[i][ni][j,0] = -np.arctan2(np.cos(H1[j]),SinDec*np.sin(H1[j]))
      elif mounts[temp['antenna1'][j]] == 'NR':
          PAs[i][ni][j,0] = -np.arctan2(np.sin(H1[j]), CT1[j] - SinDec*np.cos(H1[j])) - E1[j]
      elif mounts[temp['antenna1'][j]] == 'NL':
          PAs[i][ni][j,0] = -np.arctan2(np.sin(H1[j]), CT1[j] - SinDec*np.cos(H1[j])) + E1[j]



      if mounts[temp['antenna2'][j]] == 'AZ':
          PAs[i][ni][j,1] = -np.arctan2(np.sin(H2[j]), CT2[j] - SinDec*np.cos(H2[j]))
      elif mounts[temp['antenna2'][j]] == 'EQ':
          PAs[i][ni][j,1] = 0.0
      elif mounts[temp['antenna2'][j]] == 'XY':
          PAs[i][ni][j,1] = -np.arctan2(np.cos(H2[j]),SinDec*np.sin(H2[j]))
      elif mounts[temp['antenna2'][j]] == 'NR':
          PAs[i][ni][j,1] = -np.arctan2(np.sin(H2[j]), CT2[j] - SinDec*np.cos(H2[j])) - E2[j]
      elif mounts[temp['antenna2'][j]] == 'NL':
          PAs[i][ni][j,1] = -np.arctan2(np.sin(H2[j]), CT2[j] - SinDec*np.cos(H2[j])) + E2[j]
    
    if i==0:
      Flags.append(np.copy(temp2['flag'][0,0,:]))      


    del E1, E2, H1, H2, HAng, CT1, CT2, CH, SH, V2, Bx, By, Bz

    ms.close()  
#######################################










  print "\nSimulating Stokes I"
  if len(Iim)>0:
    ft(vis=dvis[4], model = Iim, usescratch=True)


  print "\nSimulating Stokes Q"
  if len(Qim)>0:
    ft(vis=dvis[5], model = Qim, usescratch=True)


  print "\nSimulating Stokes U"
  if len(Uim)>0:
    ft(vis=dvis[6], model = Uim, usescratch=True)


  print "\nSimulating Stokes V"
  if len(Uim)>0:
    ft(vis=dvis[7], model = Uim, usescratch=True)


  printMsg( 'Computing the correlations')
  XX = np.zeros(np.shape(dataI),dtype=np.complex128)
  YY = np.zeros(np.shape(dataI),dtype=np.complex128)
  XY = np.zeros(np.shape(dataI),dtype=np.complex128)
  YX = np.zeros(np.shape(dataI),dtype=np.complex128)

#  if corrupt:
  XXa = np.zeros(nchan,dtype=np.complex128)
  YYa = np.zeros(nchan,dtype=np.complex128)
  XYa = np.zeros(nchan,dtype=np.complex128)
  YXa = np.zeros(nchan,dtype=np.complex128)
  XXb = np.zeros(nchan,dtype=np.complex128)
  YYb = np.zeros(nchan,dtype=np.complex128)
  XYb = np.zeros(nchan,dtype=np.complex128)
  YXb = np.zeros(nchan,dtype=np.complex128)


  FouFac = 1.j*2.*np.pi*(np.pi/180./3600.)
  for i in range(len(BBs)):
   printMsg( 'Doing spw %i'%i)
   gc.collect()
   for sci,sc in enumerate(spwscans[i]):
    print 'Scan %i of %i'%(sci+1,len(spwscans[i]))
    ms.open(dvis[4],nomodify=False)
    ms.selectinit(datadescid=i)
    ms.select({'scan_number':int(sc)})

    UVs = ms.getdata(['u','v'])
    U = UVs['u'][np.newaxis,:]/(3.e8/spwFreqs[i][:,np.newaxis])
    V = UVs['v'][np.newaxis,:]/(3.e8/spwFreqs[i][:,np.newaxis])

    if len(Iim)>0:
      dataI[:] = ms.getdata(['model_data'])['model_data']
    else:
      dataI[:] = 0.0
    if len(I)>0 or ismodel:
      for modi in range(NmodComp):
        dataI[:] += ModI[i][:,modi][:,np.newaxis]*np.exp(FouFac*((RAoffset[modi]-RAoffset[0])*U + (Decoffset[modi]-Decoffset[0])*V)) 

    ms.close()

    ms.open(dvis[5],nomodify=False)
    ms.selectinit(datadescid=i)
    ms.select({'scan_number':int(sc)})
 
    if len(Qim)>0:
      dataQ[:] = ms.getdata(['model_data'])['model_data']
    else:
      dataQ[:] = 0.0
    if len(I)>0 or ismodel:
      for modi in range(NmodComp):
        dataQ[:] += ModQ[i][:,modi][:,np.newaxis]*np.exp(FouFac*((RAoffset[modi]-RAoffset[0])*U + (Decoffset[modi]-Decoffset[0])*V))

    ms.close()

    ms.open(dvis[6],nomodify=False)
    ms.selectinit(datadescid=i)
    ms.select({'scan_number':int(sc)})

    if len(Uim)>0:
      dataU[:] = ms.getdata(['model_data'])['model_data']
    else:
      dataU[:] = 0.0
    if len(I)>0 or ismodel:
      for modi in range(NmodComp):
        dataU[:] += ModU[i][:,modi][:,np.newaxis]*np.exp(FouFac*((RAoffset[modi]-RAoffset[0])*U + (Decoffset[modi]-Decoffset[0])*V))

    ms.close()

    ms.open(dvis[7],nomodify=False)
    ms.selectinit(datadescid=i)
    ms.select({'scan_number':int(sc)})

    if len(Vim)>0:
      dataV[:] = ms.getdata(['model_data'])['model_data']
    else:
      dataV[:] = 0.0
    if len(I)>0 or ismodel:
      for modi in range(NmodComp):
        dataV[:] += ModV[i][:,modi][:,np.newaxis]*np.exp(FouFac*((RAoffset[modi]-RAoffset[0])*U + (Decoffset[modi]-Decoffset[0])*V))

    ms.close()

    del U, V, UVs['u'], UVs['v']

    for j in range(ntimes):
     PA = PAs[i][sci][j,:]
     C = np.cos(PA) ; S = np.sin(PA)
     EPA = np.exp(1.j*(PA[0]+PA[1])) ; EMA = np.exp(1.j*(PA[0]-PA[1]))


  # Visibilities in the antenna frame:
     if feed == 'linear':
         XXa[:] = ((dataQ[0,:,j]+dataI[0,:,j])*C[0]-(dataU[0,:,j]-1.j*dataV[0,:,j])*S[0])*C[1] - ((dataU[0,:,j]+1.j*dataV[0,:,j])*C[0] + (dataQ[0,:,j]-dataI[0,:,j])*S[0])*S[1]
         XYa[:] = ((dataU[0,:,j]+1.j*dataV[0,:,j])*C[0]+(dataQ[0,:,j]-dataI[0,:,j])*S[0])*C[1] + ((dataQ[0,:,j]+dataI[0,:,j])*C[0] - (dataU[0,:,j]-1.j*dataV[0,:,j])*S[0])*S[1]
         YXa[:] = ((dataU[0,:,j]-1.j*dataV[0,:,j])*C[0]+(dataQ[0,:,j]+dataI[0,:,j])*S[0])*C[1] + ((dataQ[0,:,j]-dataI[0,:,j])*C[0] - (dataU[0,:,j]+1.j*dataV[0,:,j])*S[0])*S[1]
         YYa[:] = -((dataQ[0,:,j]-dataI[0,:,j])*C[0]-(dataU[0,:,j]+1.j*dataV[0,:,j])*S[0])*C[1] + ((dataU[0,:,j]-1.j*dataV[0,:,j])*C[0] + (dataQ[0,:,j]+dataI[0,:,j])*S[0])*S[1]

     if feed == 'circular':
         XXa[:] = (dataI[0,:,j] + dataV[0,:,j])*EMA 
         YYa[:] = (dataI[0,:,j] - dataV[0,:,j])/EMA 
         XYa[:] = (dataQ[0,:,j] + 1.j*dataU[0,:,j])*EPA
         YXa[:] = (dataQ[0,:,j] - 1.j*dataU[0,:,j])/EPA


  # Apply leakage:
     XX[0,:,j] = XXa + YYa*dtermsX[i][ant1[i][sc][j]]*np.conjugate(dtermsX[i][ant2[i][sc][j]]) + XYa*np.conjugate(dtermsX[i][ant2[i][sc][j]]) + YXa*dtermsX[i][ant1[i][sc][j]]
     YY[0,:,j] = YYa + XXa*dtermsY[i][ant1[i][sc][j]]*np.conjugate(dtermsY[i][ant2[i][sc][j]]) + XYa*dtermsY[i][ant1[i][sc][j]] + YXa*np.conjugate(dtermsY[i][ant2[i][sc][j]])
     XY[0,:,j] = XYa + YYa*dtermsX[i][ant1[i][sc][j]] + XXa*np.conjugate(dtermsY[i][ant2[i][sc][j]]) + YXa*dtermsX[i][ant1[i][sc][j]]*np.conjugate(dtermsY[i][ant2[i][sc][j]])
     YX[0,:,j] = YXa + XXa*dtermsY[i][ant1[i][sc][j]] + YYa*np.conjugate(dtermsX[i][ant2[i][sc][j]]) + XYa*dtermsY[i][ant1[i][sc][j]]*np.conjugate(dtermsX[i][ant2[i][sc][j]])

     del PA, C, S, EPA, EMA

# Put back into sky frame:

     if apply_parang:

       if feed == 'linear':
         XXb[:]= (C[0]*XX[0,:,j] + S[0]*YX[0,:,j])*C[1] + (C[0]*XY[0,:,j]+S[0]*YY[0,:,j])*S[1]
         YYb[:] = -(S[0]*XY[0,:,j] - C[0]*YY[0,:,j])*C[1] + (S[0]*XX[0,:,j]-C[0]*YX[0,:,j])*S[1]
         XYb[:] = (C[0]*XY[0,:,j] + S[0]*YY[0,:,j])*C[1] - (C[0]*XX[0,:,j] + S[0]*YX[0,:,j])*S[1]
         YXb[:] = -(S[0]*XX[0,:,j] - C[0]*YX[0,:,j])*C[1] - (S[0]*XY[0,:,j] - C[0]*YY[0,:,j])*S[1]

         XX[0,:,j] = XXb         
         XY[0,:,j] = XYb
         YX[0,:,j] = YXb
         YY[0,:,j] = YYb

       if feed == 'circular':
         XX[0,:,j] /= EMA
         YY[0,:,j] *= EMA
         XY[0,:,j] /= EPA
         YX[0,:,j] *= EPA



  # Save:
    print polprods[0]
    ms.open(str(dvis[0]),nomodify=False)
    ms.selectinit(datadescid=i)
    ms.select({'scan_number':int(sc)})
    aux = ms.getdata(['data'])
    aux['data'][:] = XX[:]
    ms.putdata(aux)
    ms.close()
    del aux['data']


    print polprods[1]
    ms.open(str(dvis[1]),nomodify=False)
    ms.selectinit(datadescid=i)
    ms.select({'scan_number':int(sc)})
    aux = ms.getdata(['data'])
    aux['data'][:] = XY[:]
    ms.putdata(aux)
    ms.close()
    del aux['data']


    print polprods[2]
    ms.open(str(dvis[2]),nomodify=False)
    ms.selectinit(datadescid=i)
    ms.select({'scan_number':int(sc)})
    aux = ms.getdata(['data'])
    aux['data'][:] = YX[:]
    ms.putdata(aux)
    ms.close()
    del aux['data']


    print polprods[3]
    ms.open(str(dvis[3]),nomodify=False)
    ms.selectinit(datadescid=i)
    ms.select({'scan_number':int(sc)})
    aux = ms.getdata(['data'])
    aux['data'][:] = YY[:]
    ms.putdata(aux)
    ms.close()
    del aux['data']

  gc.collect()




# The sm tool IS BROKEN for full-polarization datasets!
# Write full-pol MS manually.

  if corrupt:
   printMsg( 'Corrupting')
   for i in range(len(BBs)):
    printMsg( 'Doing spw %i'%i)
    for pri,pr in enumerate(dvis[:4]):
     print 'Polprod %s'%polprods[pri]
     sm.openfromms(pr)
     sm.setseed(seed+4*i + 16*pri)
     sm.setdata(fieldid=[sources[0]],spwid=i)
     sm.setnoise(spillefficiency=eta_s,correfficiency=eta_q,
                 antefficiency=eta_a,trx=t_rx,
                 tau=tau0,tatmos=t_sky,tground=t_ground,tcmb=2.725,
                 mode="tsys-manual",senscoeff=-1)
     sm.corrupt()
     sm.done()

# Copy into full-pol ms:
  printMsg( 'Saving')

  for i in range(len(BBs)):
    printMsg( 'Doing spw %i'%i)
    for pri,pr in enumerate(dvis[:4]):
     print 'Polprod %s'%polprods[pri]
     for sc in spwscans[i]:
      ms.open(str(pr),nomodify=False)
      ms.selectinit(datadescid=i)
      ms.select({'scan_number':int(sc)})
      aux = ms.getdata(['data'])['data'][0,:]
      ms.close()
      ms.open(vis,nomodify=False)
      ms.selectinit(datadescid=i)
      ms.select({'scan_number':int(sc)})
      data = ms.getdata(['data'])
      data['data'][pri,:] = aux
      ms.putdata(data)
      ms.close()
      del data['data'], aux



# Flag out auto-correlations and negative elevations:
  printMsg('Flagging autocorrelations and negative elevations.')
  for i in range(len(BBs)):
    ms.open(vis,nomodify=False)
    ms.selectinit(datadescid=i)
    flags = ms.getdata(['antenna1','antenna2','flag'])
    autos = flags['antenna1'] == flags['antenna2']
    flags['flag'][:,:,autos] = True
    ms.putdata(flags)
    ms.close()
    del flags['antenna1'],flags['antenna2'],flags['flag'],autos
    for si,sc in enumerate(spwscans[i]):
      ms.open(vis,nomodify=False)
      ms.selectinit(datadescid=i)
      ms.select({'scan_number':int(sc)})
      flags = ms.getdata(['flag'])
      flags['flag'][...,Flags[si]] = True
      ms.putdata(flags)
      ms.close()
      del flags['flag']

# Update mounts in ANT table:
  tb.open(vis+'/ANTENNA',nomodify=False)
  MNT = tb.getcol('MOUNT')
  for mti, mnt in enumerate(mounts):
    if mnt == 'EQ':
        MNT[mti] = 'equatorial'
    elif mnt == 'XY':
        MNT[mti] = 'X-Y'
    elif mnt == 'NR':
        MNT[mti] = 'NASMYTH-R'
    elif mnt == 'NL':
        MNT[mti] = 'NASMYTH-L'

  tb.putcol('MOUNT',MNT)
  tb.close()
  


  printMsg( 'Clearing data')
  del XX, YY, XY, YX, dataI, dataQ, dataU, dataV
  gc.collect()
  clearcal(vis)

  for dv in dvis:
    os.system('rm -rf %s'%dv)


  if export_uvf:
    printMsg('Exporting to UVFITS')
    os.system('rm -rf %s.uvf'%vis)
    exportuvfits(vis=vis,fitsfile='%s.uvf'%vis,datacolumn='data')
    import pyfits as pf
    temp = pf.open('%s.uvf'%vis,mode='update')
    for mti,mt in enumerate(mounts):
        mtidx = {'AZ':0,'EQ':1,'OR':2,'XY':3,'NR':4,'NL':5}[mt]
        temp[2].data['MNTSTA'][mti] = mtidx
    temp.flush()
    temp.close()

  print '\n DONE!\n'




if __name__=='__main__':
 polsimulate(vis,array_configuration,feed,LO,BBs,spw_width,nchan,model_image,I,Q,U,V,RM,spec_index,
   spectrum_file,incenter,incell,inbright,inwidth,H0,onsource_time,observe_time,visib_time,nscan,
   corrupt, seed, Dt_amp, Dt_noise,tau0,t_sky,t_ground,t_receiver)


