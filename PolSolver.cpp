#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdio.h>  
#include <sys/types.h>
#include <new>
#include <complex>

/* Docstrings */
static char module_docstring[] =
    "Solver of leakage terms and source polarization.";
static char setData_docstring[] =
    "Reads data pointers and arranges data.";
static char getHessian_docstring[] =
    "Computes the Hessian matrix and residuals vector, given the parameter values";



/* Available functions */
static PyObject *setData(PyObject *self, PyObject *args);
static PyObject *getHessian(PyObject *self, PyObject *args);


// Global variables:
static bool DEBUG = false;
bool doCirc;  
typedef std::complex<double> cplx64d;

int Nchan, Nvis, NSou, NAnt, NPar;
int NITER = 0;
int *A1, *A2;
double *Wgt, *VARS, *Hessian, *DerVec;
double TotFlux;
cplx64d *DATA, *COMPS;
cplx64d *DR, *DL, *EPA, *EMA;

int *PSou, *PAntR, *PAntL, *VSou, *VAntR, *VAntL;


/* Module specification */
static PyMethodDef module_methods[] = {
    {"setData", setData, METH_VARARGS, setData_docstring},
    {"getHessian", getHessian, METH_VARARGS, getHessian_docstring},
    {NULL, NULL, 0, NULL}
};


/* Initialize the module */
PyMODINIT_FUNC initPolSolver(void)
{
    PyObject *m = Py_InitModule3("PolSolver", module_methods, module_docstring);
    if (m == NULL)
        return;

/* Load `numpy` functionality. */
    import_array();


}



//////////////////////////////////
// Gets the pointers to data, metadata, Hessian and residuals vector:
static PyObject *setData(PyObject *self, PyObject *args)
{

  PyObject *DATAPy, *A1Py, *A2Py, *COMPSPy, *EPAsPy, *EMAsPy, *PSouPy, *VARSPy; 
  PyObject *PAntRPy, *PAntLPy, *WgtPy, *VSouPy, *VAntRPy, *VAntLPy, *HessianPy, *ResVecPy;

  PyObject *Err;
  int i;


  if (!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOOdb", &DATAPy, &WgtPy, &A1Py, &A2Py, &COMPSPy, &EPAsPy, &EMAsPy, &PSouPy, &PAntRPy, &PAntLPy, &VSouPy, &VAntRPy, &VAntLPy, &VARSPy, &HessianPy, &ResVecPy,&TotFlux, &doCirc))
    {printf("Failed setData! Wrong arguments!\n"); fflush(stdout); Err = Py_BuildValue("i",-1); return Err;};


  Hessian = (double *)PyArray_DATA(HessianPy);
  DerVec = (double *)PyArray_DATA(ResVecPy);
  NPar =  PyArray_SHAPE(reinterpret_cast<PyArrayObject*>(ResVecPy))[0];

// Array must be of proper full-polarization data:
  int Ndims = PyArray_NDIM(DATAPy);
  if (Ndims != 3){printf("Array must have three dimenstions!\n"); fflush(stdout); Err = Py_BuildValue("i",-2); return Err;};

  Nchan = PyArray_SHAPE(reinterpret_cast<PyArrayObject*>(DATAPy))[1];
  Nvis = PyArray_SHAPE(reinterpret_cast<PyArrayObject*>(DATAPy))[0];

  NSou = PyArray_SHAPE(reinterpret_cast<PyArrayObject*>(PSouPy))[0];
  NAnt = PyArray_SHAPE(reinterpret_cast<PyArrayObject*>(PAntRPy))[0];

// get pointers to the data:
  DATA = (cplx64d *)PyArray_DATA(DATAPy);
  A1 = (int *)PyArray_DATA(A1Py);
  A2 = (int *)PyArray_DATA(A2Py);
  COMPS = (cplx64d *)PyArray_DATA(COMPSPy);
  EPA = (cplx64d *)PyArray_DATA(EPAsPy);
  EMA = (cplx64d *)PyArray_DATA(EMAsPy);
  Wgt = (double *)PyArray_DATA(WgtPy);
  PSou = (int *)PyArray_DATA(PSouPy);
  PAntR = (int *)PyArray_DATA(PAntRPy);
  PAntL = (int *)PyArray_DATA(PAntLPy);
  VSou = (int *) PyArray_DATA(VSouPy);
  VAntR = (int *)PyArray_DATA(VAntRPy);
  VAntL = (int *)PyArray_DATA(VAntLPy);
  VARS = (double *)PyArray_DATA(VARSPy);

  DR = reinterpret_cast<cplx64d*>(&VARS[VAntR[0]]);
  DL = reinterpret_cast<cplx64d*>(&VARS[VAntL[0]]);


if (DEBUG){

  printf("%i visibs with %i channels\n",Nvis,Nchan);

  printf("SOURCES: \n");
  for (i=0; i< NSou; i++){
    printf("#%i:  PSou = %i,  VSou = %i\n",i,PSou[i],VSou[i]);
    printf("  Starting values: %.3e, %.3e\n\n",VARS[VSou[i]],VARS[VSou[i]+1]);
  };

  printf("ANTENNAS: \n");
  for (i=0; i< NAnt; i++){
    printf("#%i:  PAntR = %i,  VAntR = %i\n",i,PAntR[i],VAntR[i]);
    printf("#%i:  PAntL = %i,  VAntL = %i\n",i,PAntL[i],VAntL[i]);
    printf("  Starting values: %.3e, %.3e | %.3e, %.3e\n\n",DR[i].real(),DR[i].imag(),DL[i].real(),DL[i].imag());
  };

  int chi=1;
  printf("FIRST %i VISIBILITIES (channel %i): \n",2*NAnt,chi);
  cplx64d RR, RL, LR, LL;
  double wgtaux;
  for (i=0; i < 2*NAnt; i++){
     RR = DATA[chi*4 + i*4*Nchan];
     LL = DATA[chi*4 + i*4*Nchan+1];
     RL = DATA[chi*4 + i*4*Nchan+2];
     LR = DATA[chi*4 + i*4*Nchan+3];
     wgtaux = Wgt[i*Nchan + chi];
     printf("#%i (wgt %.2e): RR = (%-.2e,%-.2e) | RL = (%-.2e,%-.2e) | LR = (%-.2e,%-.2e) | LL = (%-.2e,%-.2e)\n",i,wgtaux,RR.real(),RR.imag(), RL.real(),RL.imag(), LR.real(),LR.imag(), LL.real(),LL.imag());




  };


};

// Return success:
  NITER = 0;
  Err = Py_BuildValue("i",0);

  return Err;


};





/////////////////////////////////////////////////////////////
// Computes Hessian, residuals vector, and Chi Square.

static PyObject *getHessian(PyObject *self, PyObject *args)
{

  PyObject *Err;

  double Ifac, ChiSq, res, F[4];
  double IRatio, Iwt;
  double dArg;
  int i,j,k,l;
  long Ndata;



  if(doCirc){ 
      F[0] = 0.; F[1] = 1.; F[2] = 1.; F[3] = 0.;
  } else {
      F[0] = 1.; F[1] = 1.; F[2] = 1.; F[3] = 1.;
  }; // Relative weights for RR, RL, LR, LL / XX, XY, YX, YY
  
  
  double currWgt, ord2;
  
  
  // Turn on/off 2nd order corrections:
  ord2 = 1.0 ;
  

  int VisPar[NPar];
  int currPar;

  int iaux1, iaux2, iaux3, iaux4;

  if (!PyArg_ParseTuple(args,"d", &dArg))
    {printf("Failed getHessian! Wrong arguments!\n"); 
     Err = Py_BuildValue("i",-1); fflush(stdout); return Err;};

  cplx64d Itot, Qtot, Utot;
  cplx64d RRc, RLc, LRc, LLc;
  cplx64d Im = cplx64d(0.,1.);

  cplx64d RR, RL, LR, LL, auxC;
  cplx64d resid[4];
  cplx64d AllDer[NPar][4];

  ChiSq = 0.0;


// set matrix and vector to zero:
  for(i=0;i<NPar;i++){
    DerVec[i] = 0.;
    for(j=0;j<NPar;j++){
      Hessian[i*NPar + j] = 0.;
    };
  };

  IRatio = 0.; Iwt = 0.; Ifac = 0.;
  Ndata = 0;
  
// Loop over all data in proper order (time, channel, polariz):
  for(i=0; i<Nvis; i++){
    for(j=0; j<Nchan; j++){

// Indices in 1D for the multi-dimensional arrays:
     iaux1 = (i*Nchan + j)*(NSou+1);
     iaux2 = 4*(j + i*Nchan);
     iaux3 = i*Nchan + j;

// Index of current parameter in Hessian matrix:
     currPar = 0;

// Reset vector of derivatives to zero:
     for(k=0;k<NPar;k++){for(iaux4=0;iaux4<4;iaux4++){AllDer[k][iaux4] = cplx64d(0.,0.);};};

// Proceed only if data are good:
     if(Wgt[iaux3]>0.0){
    
      Ndata += 1;
      
// Add Fourier transforms of all source components:
      Itot = COMPS[iaux1+NSou]/TotFlux; // cplx64d(0.,0.);
      Qtot = cplx64d(0.,0.);
      Utot = cplx64d(0.,0.);

      for(k=0; k<NSou; k++){
    //    Itot += COMPS[iaux1 + k];
        Qtot += VARS[VSou[k]]*COMPS[iaux1 + k]/TotFlux;
        Utot += VARS[VSou[k]+1]*COMPS[iaux1 + k]/TotFlux;
        iaux4 = PSou[k];
        if(iaux4>=0){
          VisPar[currPar] = iaux4;
          VisPar[currPar+1] = iaux4+1;
          currPar += 2;
// Derivative w.r.t. Q for all corr products (RR, RL, LR, LL):
          AllDer[iaux4][0] = COMPS[iaux1+k]*(std::conj(DR[A2[i]])*EPA[i] + DR[A1[i]]/EPA[i]);
          AllDer[iaux4][1] = COMPS[iaux1+k]*(EPA[i] + ord2*std::conj(DL[A2[i]])*DR[A1[i]]/EPA[i]);
          AllDer[iaux4][2] = COMPS[iaux1+k]*(1./EPA[i] + ord2*std::conj(DR[A2[i]])*DL[A1[i]]*EPA[i]);
          AllDer[iaux4][3] = COMPS[iaux1+k]*(std::conj(DL[A2[i]])/EPA[i] + DL[A1[i]]*EPA[i]);
// Derivative w.r.t. U:
          AllDer[iaux4+1][0] = Im*COMPS[iaux1+k]*(std::conj(DR[A2[i]])*EPA[i] - DR[A1[i]]/EPA[i]);
          AllDer[iaux4+1][1] = Im*COMPS[iaux1+k]*(EPA[i] - ord2*std::conj(DL[A2[i]])*DR[A1[i]]/EPA[i]);
          AllDer[iaux4+1][2] = Im*COMPS[iaux1+k]*(-1./EPA[i] + ord2*std::conj(DR[A2[i]])*DL[A1[i]]*EPA[i]);
          AllDer[iaux4+1][3] = Im*COMPS[iaux1+k]*(-std::conj(DL[A2[i]])/EPA[i] + DL[A1[i]]*EPA[i]);
          
     ////////////////////////////    
     // FOR TESTING:     
     //     AllDer[iaux4][1] = cplx64d((COMPS[iaux1+k]*EPA[i]).real(),0.);
     //     AllDer[iaux4][2] = cplx64d((COMPS[iaux1+k]/EPA[i]).real(),0.);
     //     AllDer[iaux4+1][1] = cplx64d(0.,(COMPS[iaux1+k]*EPA[i]).imag());
     //     AllDer[iaux4+1][2] = cplx64d(0.,-(COMPS[iaux1+k]/EPA[i]).imag());
     ////////////////////////////    

          
          
        };
      };

      Ifac = std::abs(Itot);
      currWgt = Wgt[iaux3]*Ifac;

      // Stokes I (i.e., (RR + LL)/2 ) from the data, corrected by inverse of Dt matrices:
      Iwt += 0.5*std::abs((DATA[iaux2]*(1. + std::conj(DL[A2[i]])*DL[A1[i]])/EMA[i] + DATA[iaux2+1]*EMA[i]*(1. + std::conj(DR[A2[i]])*DR[A1[i]]) - DATA[iaux2+2]*(std::conj(DR[A2[i]])/EMA[i] + DL[A1[i]]*EMA[i]) - DATA[iaux2+3]*(DR[A1[i]]/EMA[i] + std::conj(DL[A2[i]])*EMA[i]))/(1.-DR[A1[i]]*DL[A1[i]])/(1.-std::conj(DR[A2[i]])*std::conj(DL[A2[i]])));
    
      // Stokes I from the model:
      IRatio += Ifac;
      
      
// Correlation products in antenna frame with no leakage:
      RRc = Itot*EMA[i];
      RLc = (Qtot + Im*Utot)*EPA[i];
      LRc = (Qtot - Im*Utot)/EPA[i];
      LLc = Itot/EMA[i];

// Apply leakage to model visibilities:
      RR = (RRc + RLc*std::conj(DR[A2[i]]) + LRc*DR[A1[i]] + ord2*LLc*std::conj(DR[A2[i]])*DR[A1[i]]);
      RL = (RLc + RRc*std::conj(DL[A2[i]]) + LLc*DR[A1[i]] + ord2*LRc*std::conj(DL[A2[i]])*DR[A1[i]]);
      LR = (LRc + RRc*DL[A1[i]] + LLc*std::conj(DR[A2[i]]) + ord2*RLc*std::conj(DR[A2[i]])*DL[A1[i]]);
      LL = (LLc + LRc*std::conj(DL[A2[i]]) + RLc*DL[A1[i]] + ord2*RRc*std::conj(DL[A2[i]])*DL[A1[i]]);

      if(DEBUG && NITER<2 && i<10 && j==0){
        printf("Perfect %i: RR = (%.2e, %.2e); RL = (%.2e, %.2e); LR = (%.2e, %.2e); LL = (%.2e, %.2e)\n",i,RRc.real(),RRc.imag(),RLc.real(),RLc.imag(),LRc.real(),LRc.imag(),LLc.real(),LLc.imag());
        printf("Corrupt %i: RR = (%.2e, %.2e); RL = (%.2e, %.2e); LR = (%.2e, %.2e); LL = (%.2e, %.2e)\n",i,RR.real(),RR.imag(),RL.real(),RL.imag(),LR.real(),LR.imag(),LL.real(),LL.imag());

      };

// Residuals for each correlation product:
      resid[0] =  (DATA[iaux2  ] - RR);
      resid[1] =  (DATA[iaux2+2] - RL); 
      resid[2] =  (DATA[iaux2+3] - LR);
      resid[3] =  (DATA[iaux2+1] - LL);

// Vector of derivative*residuals (first, for source components):
      for(k=0; k<NSou; k++){
        iaux4 = PSou[k];
        if(iaux4>=0){
          res = 0.0;
          for(l=0;l<4;l++){
            res += AllDer[iaux4][l].real()*resid[l].real()*F[l];
            res += AllDer[iaux4][l].imag()*resid[l].imag()*F[l];
          }; 
          DerVec[iaux4] += res*currWgt; res = 0.0;
          for(l=0;l<4;l++){
            res += AllDer[iaux4+1][l].real()*resid[l].real()*F[l];
            res += AllDer[iaux4+1][l].imag()*resid[l].imag()*F[l];
          }; 
          DerVec[iaux4+1] += res*currWgt;
        };
      };

// Derivatives w.r.t. Dterms:
      iaux4 = PAntR[A1[i]];
      if(iaux4>=0){
        VisPar[currPar] = iaux4;
        VisPar[currPar+1] = iaux4+1;
        currPar += 2;
        
        AllDer[iaux4][1] = LLc + ord2*LRc*std::conj(DL[A2[i]]);
        AllDer[iaux4+1][1] = Im*AllDer[iaux4][1];

        AllDer[iaux4][0] = LRc + ord2*LLc*std::conj(DR[A2[i]]);
        AllDer[iaux4+1][0] = Im*AllDer[iaux4][0];

        
        res  = AllDer[iaux4][0].real()*resid[0].real()*F[0];
        res += AllDer[iaux4][0].imag()*resid[0].imag()*F[0];
        res += AllDer[iaux4][1].real()*resid[1].real()*F[1];
        res += AllDer[iaux4][1].imag()*resid[1].imag()*F[1];
        DerVec[iaux4] += res*currWgt;

        res  = AllDer[iaux4+1][0].real()*resid[0].real()*F[0];
        res += AllDer[iaux4+1][0].imag()*resid[0].imag()*F[0];
        res += AllDer[iaux4+1][1].real()*resid[1].real()*F[1];
        res += AllDer[iaux4+1][1].imag()*resid[1].imag()*F[1];
        DerVec[iaux4+1] += res*currWgt;
      };

      iaux4 = PAntL[A1[i]];
      if(iaux4>=0){
        VisPar[currPar] = iaux4;
        VisPar[currPar+1] = iaux4+1;
        currPar += 2;
        
        AllDer[iaux4][2] = RRc + ord2*RLc*std::conj(DR[A2[i]]);
        AllDer[iaux4+1][2] = Im*AllDer[iaux4][2];

        AllDer[iaux4][3] = RLc + ord2*RRc*std::conj(DL[A2[i]]);
        AllDer[iaux4+1][3] = Im*AllDer[iaux4][3];

        res  = AllDer[iaux4][2].real()*resid[2].real()*F[2];
        res += AllDer[iaux4][2].imag()*resid[2].imag()*F[2];
        res += AllDer[iaux4][3].real()*resid[3].real()*F[3];
        res += AllDer[iaux4][3].imag()*resid[3].imag()*F[3];
        DerVec[iaux4] += res*currWgt;

        res  = AllDer[iaux4+1][2].real()*resid[2].real()*F[2];
        res += AllDer[iaux4+1][2].imag()*resid[2].imag()*F[2];
        res += AllDer[iaux4+1][3].real()*resid[3].real()*F[3];
        res += AllDer[iaux4+1][3].imag()*resid[3].imag()*F[3];
        DerVec[iaux4+1] += res*currWgt;
      };

      iaux4 = PAntR[A2[i]];
      if(iaux4>=0){
        VisPar[currPar] = iaux4;
        VisPar[currPar+1] = iaux4+1;
        currPar += 2;
        
        AllDer[iaux4][2] = LLc + ord2*RLc*DL[A1[i]];
        AllDer[iaux4+1][2] = -Im*AllDer[iaux4][2];

        
        AllDer[iaux4][0] = RLc + ord2*LLc*DR[A1[i]];
        AllDer[iaux4+1][0] = -Im*AllDer[iaux4][0];

        res  = AllDer[iaux4][0].real()*resid[0].real()*F[0];
        res += AllDer[iaux4][0].imag()*resid[0].imag()*F[0];
        res += AllDer[iaux4][2].real()*resid[2].real()*F[2];
        res += AllDer[iaux4][2].imag()*resid[2].imag()*F[2];
        DerVec[iaux4] += res*currWgt;

        res  = AllDer[iaux4+1][0].real()*resid[0].real()*F[0];
        res += AllDer[iaux4+1][0].imag()*resid[0].imag()*F[0];
        res += AllDer[iaux4+1][2].real()*resid[2].real()*F[2];
        res += AllDer[iaux4+1][2].imag()*resid[2].imag()*F[2];
        DerVec[iaux4+1] += res*currWgt;
      };

      iaux4 = PAntL[A2[i]];
      if(iaux4>=0){
        VisPar[currPar] = iaux4;
        VisPar[currPar+1] = iaux4+1;
        currPar += 2;

        AllDer[iaux4][1] = RRc + ord2*LRc*DR[A1[i]];
        AllDer[iaux4+1][1] = -Im*AllDer[iaux4][1];

        AllDer[iaux4][3] = LRc + ord2*RRc*DL[A1[i]];
        AllDer[iaux4+1][3] = -Im*AllDer[iaux4][3];

        res  = AllDer[iaux4][1].real()*resid[1].real()*F[1];
        res += AllDer[iaux4][1].imag()*resid[1].imag()*F[1];
        res += AllDer[iaux4][3].real()*resid[3].real()*F[3];
        res += AllDer[iaux4][3].imag()*resid[3].imag()*F[3];
        DerVec[iaux4] += res*currWgt;

        res  = AllDer[iaux4+1][1].real()*resid[1].real()*F[1];
        res += AllDer[iaux4+1][1].imag()*resid[1].imag()*F[1];
        res += AllDer[iaux4+1][3].real()*resid[3].real()*F[3];
        res += AllDer[iaux4+1][3].imag()*resid[3].imag()*F[3];
        DerVec[iaux4+1] += res*currWgt;
      };

      // Add up to the Chi Square:
      for(iaux4=0;iaux4<4;iaux4++){
        ChiSq += resid[iaux4].real()*resid[iaux4].real()*currWgt*F[iaux4];
        ChiSq += resid[iaux4].imag()*resid[iaux4].imag()*currWgt*F[iaux4];
      };

      
// Add up the product of derivatives to the Hessian matrix:
     for(k=0;k<currPar;k++){
       for(l=0;l<=k;l++){
         for(iaux4=0;iaux4<4;iaux4++){
           res  = AllDer[VisPar[k]][iaux4].real()*AllDer[VisPar[l]][iaux4].real()*F[iaux4];
           res += AllDer[VisPar[k]][iaux4].imag()*AllDer[VisPar[l]][iaux4].imag()*F[iaux4];
           Hessian[VisPar[k]*NPar + VisPar[l]] += res*currWgt;
           if(k!=l){Hessian[VisPar[l]*NPar + VisPar[k]] = Hessian[VisPar[k]*NPar + VisPar[l]];};
         };
       };
     };

     };

    };
  };

// Reduced Chi Square:
  ChiSq /= (double)(Ndata);

// Ratio of Stokes I between model and data:  
  IRatio /= Iwt;  

// Update ratio:  
  if (std::abs(IRatio - 1.0) > 1.e-3){TotFlux *= IRatio;}; 

  
// Return the reduced ChiSquare (and the I ratio):
  Err = Py_BuildValue("[d,d]",ChiSq,TotFlux);

  NITER += 1;
  printf("\r ITER %i. ChiSq %.3e ; Model Factor: %.3e",NITER,ChiSq,IRatio);
  fflush(stdout);
  return Err;

};
