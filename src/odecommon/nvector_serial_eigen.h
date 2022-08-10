#pragma once

#include <sundials/sundials_nvector.h>

#define NV_CONTENT_S(v)  ( (VectorReal*)(v->content) )
#define NV_LENGTH_S(v)   ( NV_CONTENT_S(v)->size() )
#define NV_DATA_S(v)     ( NV_CONTENT_S(v)->data() )
#define NV_Ith_S(v,i)	 ( (*NV_CONTENT_S(v))(i) )
#define EIGV(v)			 (*(VectorReal*)(v->content))

N_Vector N_VNew_Eigen(sunindextype vec_length);
N_Vector N_VNewEmpty_Eigen(sunindextype vec_length);
N_Vector* N_VCloneVectorArray_Eigen(int count, N_Vector w);
N_Vector* N_VCloneVectorArrayEmpty_Eigen(int count, N_Vector w);
void N_VDestroyVectorArray_Eigen(N_Vector* vs, int count);
sunindextype N_VGetLength_Eigen(N_Vector v);
void N_VPrint_Eigen(N_Vector v);
void N_VPrintFile_Eigen(N_Vector v, FILE *outfile);

N_Vector_ID N_VGetVectorID_Eigen(N_Vector v);
N_Vector N_VCloneEmpty_Eigen(N_Vector w);
N_Vector N_VClone_Eigen(N_Vector w);
void N_VDestroy_Eigen(N_Vector v);
void N_VSpace_Eigen(N_Vector v, sunindextype *lrw, sunindextype *liw);
realtype *N_VGetArrayPointer_Eigen(N_Vector v);

/* standard vector operations */
void N_VAdd_Eigen(N_Vector x, N_Vector y);
void N_VLinearSum_Eigen(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z);
void N_VConst_Eigen(realtype c, N_Vector z);
void N_VProd_Eigen(N_Vector x, N_Vector y, N_Vector z);
void N_VDiv_Eigen(N_Vector x, N_Vector y, N_Vector z);
void N_VScale_Eigen(realtype c, N_Vector x, N_Vector z);
void N_VAbs_Eigen(N_Vector x, N_Vector z);
void N_VInv_Eigen(N_Vector x, N_Vector z);
void N_VAddConst_Eigen(N_Vector x, realtype b, N_Vector z);
realtype N_VDotProd_Eigen(N_Vector x, N_Vector y);
realtype N_VMaxNorm_Eigen(N_Vector x);
realtype N_VWrmsNorm_Eigen(N_Vector x, N_Vector w);
realtype N_VWrmsNormMask_Eigen(N_Vector x, N_Vector w, N_Vector id);
realtype N_VMin_Eigen(N_Vector x);
realtype N_VWL2Norm_Eigen(N_Vector x, N_Vector w);
realtype N_VL1Norm_Eigen(N_Vector x);
void N_VCompare_Eigen(realtype c, N_Vector x, N_Vector z);
booleantype N_VInvTest_Eigen(N_Vector x, N_Vector z);
booleantype N_VConstrMask_Eigen(N_Vector c, N_Vector x, N_Vector m);
realtype N_VMinQuotient_Eigen(N_Vector num, N_Vector denom);

/* fused vector operations */
int N_VLinearCombination_Eigen(int nvec, realtype* c, N_Vector* V,
                                                N_Vector z);
int N_VScaleAddMulti_Eigen(int nvec, realtype* a, N_Vector x,
                                            N_Vector* Y, N_Vector* Z);
int N_VDotProdMulti_Eigen(int nvec, N_Vector x,
                                           N_Vector* Y, realtype* dotprods);

/* vector array operations */
int N_VLinearSumVectorArray_Eigen(int nvec, 
                                                   realtype a, N_Vector* X,
                                                   realtype b, N_Vector* Y,
                                                   N_Vector* Z);
int N_VScaleVectorArray_Eigen(int nvec, realtype* c,
                                               N_Vector* X, N_Vector* Z);
int N_VConstVectorArray_Eigen(int nvecs, realtype c,
                                               N_Vector* Z);
int N_VWrmsNormVectorArray_Eigen(int nvecs, N_Vector* X,
                                                  N_Vector* W, realtype* nrm);
int N_VWrmsNormMaskVectorArray_Eigen(int nvecs, N_Vector* X,
                                                      N_Vector* W, N_Vector id,
                                                      realtype* nrm);
int N_VScaleAddMultiVectorArray_Eigen(int nvec, int nsum,
                                                       realtype* a,
                                                       N_Vector* X,
                                                       N_Vector** Y,
                                                       N_Vector** Z);
int N_VLinearCombinationVectorArray_Eigen(int nvec, int nsum,
                                                           realtype* c,
                                                           N_Vector** X,
                                                           N_Vector* Z);

/* OPTIONAL local reduction kernels (no parallel communication) */
realtype N_VWSqrSumLocal_Eigen(N_Vector x, N_Vector w);
realtype N_VWSqrSumMaskLocal_Eigen(N_Vector x, N_Vector w, N_Vector id);
  
/*
 * -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * -----------------------------------------------------------------
 */

int N_VEnableFusedOps_Eigen(N_Vector v, booleantype tf);

int N_VEnableLinearCombination_Eigen(N_Vector v, booleantype tf);
int N_VEnableScaleAddMulti_Eigen(N_Vector v, booleantype tf);
int N_VEnableDotProdMulti_Eigen(N_Vector v, booleantype tf);

int N_VEnableLinearSumVectorArray_Eigen(N_Vector v, booleantype tf);
int N_VEnableScaleVectorArray_Eigen(N_Vector v, booleantype tf);
int N_VEnableConstVectorArray_Eigen(N_Vector v, booleantype tf);
int N_VEnableWrmsNormVectorArray_Eigen(N_Vector v, booleantype tf);
int N_VEnableWrmsNormMaskVectorArray_Eigen(N_Vector v, booleantype tf);
int N_VEnableScaleAddMultiVectorArray_Eigen(N_Vector v, booleantype tf);
int N_VEnableLinearCombinationVectorArray_Eigen(N_Vector v, booleantype tf);
