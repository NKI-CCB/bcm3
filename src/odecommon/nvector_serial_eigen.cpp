#include "Utils.h"

#include "nvector_serial_eigen.h"
#include <sundials/sundials_math.h>

#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)


/* Private functions for special cases of vector operations */
static void VScaleSum_Eigen(realtype c, N_Vector x, N_Vector y, N_Vector z);  /* z=c(x+y)  */
static void VScaleDiff_Eigen(realtype c, N_Vector x, N_Vector y, N_Vector z); /* z=c(x-y)  */
static void VLin1_Eigen(realtype a, N_Vector x, N_Vector y, N_Vector z);      /* z=ax+y    */
static void VLin2_Eigen(realtype a, N_Vector x, N_Vector y, N_Vector z);      /* z=ax-y    */

/* Private functions for special cases of vector array operations */
static int VSumVectorArray_Eigen(int nvec, N_Vector* X, N_Vector* Y, N_Vector* Z);                   /* Z=X+Y     */
static int VDiffVectorArray_Eigen(int nvec, N_Vector* X, N_Vector* Y, N_Vector* Z);                  /* Z=X-Y     */
static int VScaleSumVectorArray_Eigen(int nvec, realtype c, N_Vector* X, N_Vector* Y, N_Vector* Z);  /* Z=c(X+Y)  */
static int VScaleDiffVectorArray_Eigen(int nvec, realtype c, N_Vector* X, N_Vector* Y, N_Vector* Z); /* Z=c(X-Y)  */
static int VLin1VectorArray_Eigen(int nvec, realtype a, N_Vector* X, N_Vector* Y, N_Vector* Z);      /* Z=aX+Y    */
static int VLin2VectorArray_Eigen(int nvec, realtype a, N_Vector* X, N_Vector* Y, N_Vector* Z);      /* Z=aX-Y    */
static int VaxpyVectorArray_Eigen(int nvec, realtype a, N_Vector* X, N_Vector* Y);                    /* Y <- aX+Y */

N_Vector_ID N_VGetVectorID_Eigen(N_Vector v)
{
	return SUNDIALS_NVEC_EIGEN;
}

N_Vector N_VNewEmpty_Eigen(sunindextype length)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty();
  if (v == NULL) return(NULL);

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid     = N_VGetVectorID_Eigen;
  v->ops->nvclone           = N_VClone_Eigen;
  v->ops->nvcloneempty      = N_VCloneEmpty_Eigen;
  v->ops->nvdestroy         = N_VDestroy_Eigen;
  v->ops->nvspace           = N_VSpace_Eigen;
  v->ops->nvgetarraypointer = N_VGetArrayPointer_Eigen;
  v->ops->nvsetarraypointer = NULL;
  v->ops->nvgetlength       = N_VGetLength_Eigen;

  /* standard vector operations */
  v->ops->nvadd			 = N_VAdd_Eigen;
  v->ops->nvlinearsum    = N_VLinearSum_Eigen;
  v->ops->nvconst        = N_VConst_Eigen;
  v->ops->nvprod         = N_VProd_Eigen;
  v->ops->nvdiv          = N_VDiv_Eigen;
  v->ops->nvscale        = N_VScale_Eigen;
  v->ops->nvabs          = N_VAbs_Eigen;
  v->ops->nvinv          = N_VInv_Eigen;
  v->ops->nvaddconst     = N_VAddConst_Eigen;
  v->ops->nvdotprod      = N_VDotProd_Eigen;
  v->ops->nvmaxnorm      = N_VMaxNorm_Eigen;
  v->ops->nvwrmsnormmask = N_VWrmsNormMask_Eigen;
  v->ops->nvwrmsnorm     = N_VWrmsNorm_Eigen;
  v->ops->nvmin          = N_VMin_Eigen;
  v->ops->nvwl2norm      = N_VWL2Norm_Eigen;
  v->ops->nvl1norm       = N_VL1Norm_Eigen;
  v->ops->nvcompare      = N_VCompare_Eigen;
  v->ops->nvinvtest      = N_VInvTest_Eigen;
  v->ops->nvconstrmask   = N_VConstrMask_Eigen;
  v->ops->nvminquotient  = N_VMinQuotient_Eigen;

  /* fused and vector array operations are disabled (NULL) by default */

  /* local reduction operations */
  v->ops->nvdotprodlocal     = N_VDotProd_Eigen;
  v->ops->nvmaxnormlocal     = N_VMaxNorm_Eigen;
  v->ops->nvminlocal         = N_VMin_Eigen;
  v->ops->nvl1normlocal      = N_VL1Norm_Eigen;
  v->ops->nvinvtestlocal     = N_VInvTest_Eigen;
  v->ops->nvconstrmasklocal  = N_VConstrMask_Eigen;
  v->ops->nvminquotientlocal = N_VMinQuotient_Eigen;
  v->ops->nvwsqrsumlocal     = N_VWSqrSumLocal_Eigen;
  v->ops->nvwsqrsummasklocal = N_VWSqrSumMaskLocal_Eigen;

  /* debugging functions */
  v->ops->nvprint     = N_VPrint_Eigen;
  v->ops->nvprintfile = N_VPrintFile_Eigen;

  /* Create content */
  VectorReal* content = new VectorReal(length);
  if (content == NULL) {
	  N_VDestroy(v);
	  return(NULL);
  }

  /* Attach content */
  v->content = content;

  return(v);
}

N_Vector N_VNew_Eigen(sunindextype length)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Eigen(length);
  if (v == NULL) return(NULL);

  if (length > 0) {
	  EIGV(v).resize(length);
  }

  return(v);
}

N_Vector* N_VCloneVectorArray_Eigen(int count, N_Vector w)
{
	if (count <= 0) return(NULL);

	N_Vector* vs = (N_Vector*) malloc(count * sizeof(N_Vector));
	if(vs == NULL) return(NULL);

	for (int j = 0; j < count; j++) {
		vs[j] = N_VClone_Eigen(w);
		if (vs[j] == NULL) {
			N_VDestroyVectorArray_Eigen(vs, j-1);
			return(NULL);
		}
	}

	return(vs);
}

/* ----------------------------------------------------------------------------
 * Function to create an array of new serial vectors with NULL data array.
 */

N_Vector* N_VCloneVectorArrayEmpty_Eigen(int count, N_Vector w)
{
	if (count <= 0) return(NULL);

	N_Vector* vs = (N_Vector*) malloc(count * sizeof(N_Vector));
	if(vs == NULL) return(NULL);

	for (int j = 0; j < count; j++) {
		vs[j] = NULL;
		vs[j] = N_VCloneEmpty_Eigen(w);
		if (vs[j] == NULL) {
			N_VDestroyVectorArray_Eigen(vs, j-1);
			return(NULL);
		}
	}

	return(vs);
}

/* ----------------------------------------------------------------------------
 * Function to free an array created with N_VCloneVectorArray_Eigen
 */

void N_VDestroyVectorArray_Eigen(N_Vector* vs, int count)
{
  for (int j = 0; j < count; j++) {
	  N_VDestroy_Eigen(vs[j]);
  }

  free(vs); vs = NULL;

  return;
}

/* ----------------------------------------------------------------------------
 * Function to return number of vector elements
 */
sunindextype N_VGetLength_Eigen(N_Vector v)
{
	return EIGV(v).size();
}

/* ----------------------------------------------------------------------------
 * Function to print the a serial vector to stdout
 */

void N_VPrint_Eigen(N_Vector x)
{
	N_VPrintFile_Eigen(x, stdout);
}

/* ----------------------------------------------------------------------------
 * Function to print the a serial vector to outfile
 */

void N_VPrintFile_Eigen(N_Vector x, FILE* outfile)
{
  for (int i = 0; i < EIGV(x).size(); i++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(outfile, "%35.32Lg\n", EIGV(x)(i));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(outfile, "%19.16g\n", EIGV(x)(i));
#else
    fprintf(outfile, "%11.8g\n", EIGV(x)(i));
#endif
  }
  fprintf(outfile, "\n");

  return;
}

/*
 * -----------------------------------------------------------------
 * implementation of vector operations
 * -----------------------------------------------------------------
 */

N_Vector N_VCloneEmpty_Eigen(N_Vector w)
{
	N_Vector v;

	if (w == NULL) return(NULL);

	/* Create vector */
	v = NULL;
	v = N_VNewEmpty();
	if (v == NULL) return(NULL);

	/* Attach operations */
	if (N_VCopyOps(w, v)) { N_VDestroy(v); return(NULL); }

	/* Create content */
	VectorReal* content = new VectorReal(NV_LENGTH_S(w));
	if (content == NULL) { N_VDestroy(v); return(NULL); }

	/* Attach content */
	v->content = content;

	return(v);
}

N_Vector N_VClone_Eigen(N_Vector w)
{
	N_Vector v;
	realtype *data;
	sunindextype length;

	v = NULL;
	v = N_VCloneEmpty_Eigen(w);
	if (v == NULL) return(NULL);
	(*NV_CONTENT_S(v)) = (*NV_CONTENT_S(w));
	return(v);
}

void N_VDestroy_Eigen(N_Vector v)
{
	if (v == NULL) return;

	/* free content */
	if (v->content != NULL) {
		delete (VectorReal*)v->content;
		v->content = NULL;
	}

	/* free ops and vector */
	if (v->ops != NULL) { free(v->ops); v->ops = NULL; }
	free(v); v = NULL;

	return;
}

void N_VSpace_Eigen(N_Vector v, sunindextype *lrw, sunindextype *liw)
{
	*lrw = NV_LENGTH_S(v);
	*liw = 1;
}

realtype *N_VGetArrayPointer_Eigen(N_Vector v)
{
	return (realtype *)NV_DATA_S(v);
}

void N_VAdd_Eigen(N_Vector x, N_Vector y)
{
	EIGV(x).noalias() += EIGV(y);
}

void N_VLinearSum_Eigen(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
	booleantype test;
	if ((b == ONE) && (z == y)) {    /* BLAS usage: axpy y <- ax+y */
		EIGV(y) += a * EIGV(x);
	} else if ((a == ONE) && (z == x)) {    /* BLAS usage: axpy x <- by+x */
		EIGV(x) += b * EIGV(y);
	} else if ((a == ONE) && (b == ONE)) {
		EIGV(z) = EIGV(x) + EIGV(y);
	} else if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE))) {
		N_Vector v1 = test ? y : x;
		N_Vector v2 = test ? x : y;
		EIGV(z) = EIGV(v2) - EIGV(v1);
	} else if ((test = (a == ONE)) || (b == ONE)) {
		realtype c  = test ? b : a;
		N_Vector v1 = test ? y : x;
		N_Vector v2 = test ? x : y;
		VLin1_Eigen(c, v1, v2, z);
	} else if ((test = (a == -ONE)) || (b == -ONE)) {
		realtype c = test ? b : a;
		N_Vector v1 = test ? y : x;
		N_Vector v2 = test ? x : y;
		VLin2_Eigen(c, v1, v2, z);
	} else if (a == b) {
		VScaleSum_Eigen(a, x, y, z);
	} else if (a == -b) {
		VScaleDiff_Eigen(a, x, y, z);
	} else {
		EIGV(z) = a * EIGV(x) + b * EIGV(y);
	}
}

void N_VConst_Eigen(realtype c, N_Vector z)
{
	EIGV(z).setConstant(c);
}

void N_VProd_Eigen(N_Vector x, N_Vector y, N_Vector z)
{
	EIGV(z) = EIGV(x).cwiseProduct(EIGV(y));
}

void N_VDiv_Eigen(N_Vector x, N_Vector y, N_Vector z)
{
	EIGV(z) = EIGV(x).cwiseQuotient(EIGV(y));
}

void N_VScale_Eigen(realtype c, N_Vector x, N_Vector z)
{
	if (z == x) {  /* BLAS usage: scale x <- cx */
		EIGV(x) *= c;
	} else if (c == ONE) {
		EIGV(z) = EIGV(x);
	} else if (c == -ONE) {
		EIGV(z) = -EIGV(x);
	} else {
		EIGV(z).noalias() = c * EIGV(x);
	}
}

void N_VAbs_Eigen(N_Vector x, N_Vector z)
{
	EIGV(z) = EIGV(x).cwiseAbs();
}

void N_VInv_Eigen(N_Vector x, N_Vector z)
{
	EIGV(z) = EIGV(x).cwiseInverse();
}

void N_VAddConst_Eigen(N_Vector x, realtype b, N_Vector z)
{
	EIGV(z).array() = EIGV(x).array() + b;
}

realtype N_VDotProd_Eigen(N_Vector x, N_Vector y)
{
	return EIGV(x).dot(EIGV(y));
}

realtype N_VMaxNorm_Eigen(N_Vector x)
{
	return EIGV(x).maxCoeff();
}

realtype N_VWrmsNorm_Eigen(N_Vector x, N_Vector w)
{
	return(SUNRsqrt(N_VWSqrSumLocal_Eigen(x, w)/(NV_LENGTH_S(x))));
}

realtype N_VWSqrSumLocal_Eigen(N_Vector x, N_Vector w)
{
	realtype sum = ZERO;
	for (int i = 0; i < EIGV(x).size(); i++) {
		realtype p = EIGV(x)[i] * EIGV(w)[i];
		sum += p * p;
	}
	return(sum);
}

realtype N_VWrmsNormMask_Eigen(N_Vector x, N_Vector w, N_Vector id)
{
	return(SUNRsqrt(N_VWSqrSumMaskLocal_Eigen(x, w, id) / (NV_LENGTH_S(x))));
}

realtype N_VWSqrSumMaskLocal_Eigen(N_Vector x, N_Vector w, N_Vector id)
{
	realtype sum = ZERO;
	for (int i = 0; i < EIGV(x).size(); i++) {
		if (EIGV(id)[i] > ZERO) {
			realtype p = EIGV(x)[i] * EIGV(w)[i];
			sum += p * p;
		}
	}
	return(sum);
}

realtype N_VMin_Eigen(N_Vector x)
{
	return EIGV(x).minCoeff();
}

realtype N_VWL2Norm_Eigen(N_Vector x, N_Vector w)
{
	return EIGV(x).cwiseProduct(EIGV(w)).norm();
}

realtype N_VL1Norm_Eigen(N_Vector x)
{
	return EIGV(x).cwiseAbs().sum();
}

void N_VCompare_Eigen(realtype c, N_Vector x, N_Vector z)
{
	for (int i = 0; i < EIGV(x).size(); i++) {
		EIGV(z)(i) = (SUNRabs(EIGV(x)(i)) >= c) ? ONE : ZERO;
	}
}

booleantype N_VInvTest_Eigen(N_Vector x, N_Vector z)
{
	booleantype no_zero_found = SUNTRUE;
	for (int i = 0; i < EIGV(x).size(); i++) {
		if (EIGV(x)(i) == ZERO)
			no_zero_found = SUNFALSE;
		else
			EIGV(z)(i) = ONE / EIGV(x)(i);
	}
	return no_zero_found;
}

booleantype N_VConstrMask_Eigen(N_Vector c, N_Vector x, N_Vector m)
{
	realtype temp = ZERO;
	for (int i = 0; i < EIGV(x).size(); i++) {
		EIGV(m)(i) = ZERO;

		/* Continue if no constraints were set for the variable */
		if (EIGV(c)(i) == ZERO)
			continue;

		/* Check if a set constraint has been violated */
		bool test = (SUNRabs(EIGV(c)(i)) > ONEPT5 && EIGV(x)(i) * EIGV(c)(i) <= ZERO) ||
			        (SUNRabs(EIGV(c)(i)) > HALF   && EIGV(x)(i) * EIGV(c)(i) < ZERO);
		if (test) {
			temp = EIGV(m)(i) = ONE;
		}
	}
	return (temp == ONE) ? SUNFALSE : SUNTRUE;
}

realtype N_VMinQuotient_Eigen(N_Vector num, N_Vector denom)
{
  booleantype notEvenOnce = SUNTRUE;
  realtype min = BIG_REAL;

  for (int i = 0; i < EIGV(num).size(); i++) {
    if (EIGV(denom)(i) == ZERO) continue;
    else {
      if (!notEvenOnce) {
		  min = SUNMIN(min, EIGV(num)(i)  / EIGV(denom)(i));
      } else {
	    min = EIGV(num)(i) / EIGV(denom)(i);
        notEvenOnce = SUNFALSE;
      }
    }
  }
  
  return(min);
}


/*
 * -----------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------
 */

int N_VLinearCombination_Eigen(int nvec, realtype* c, N_Vector* X, N_Vector z)
{
  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VScale */
  if (nvec == 1) {
    N_VScale_Eigen(c[0], X[0], z);
    return(0);
  }

  /* should have called N_VLinearSum */
  if (nvec == 2) {
    N_VLinearSum_Eigen(c[0], X[0], c[1], X[1], z);
    return(0);
  }

  /*
   * X[0] += c[i]*X[i], i = 1,...,nvec-1
   */
  if ((X[0] == z) && (c[0] == ONE)) {
    for (int i = 1; i < nvec; i++) {
		EIGV(z) += c[i] * EIGV(X[i]);
    }
    return(0);
  }

  /*
   * X[0] = c[0] * X[0] + sum{ c[i] * X[i] }, i = 1,...,nvec-1
   */
	if (X[0] == z) {
		EIGV(z) *= c[0];
		for (int i = 1; i < nvec; i++) {
			EIGV(z) += c[i] * EIGV(X[i]);
		}
		return(0);
	}

  /*
   * z = sum{ c[i] * X[i] }, i = 0,...,nvec-1
   */
	EIGV(z) = c[0] * EIGV(X[0]);
	for (int i = 1; i < nvec; i++) {
		EIGV(z) += c[i] * EIGV(X[i]);
	}

	return(0);
}

int N_VScaleAddMulti_Eigen(int nvec, realtype* a, N_Vector x, N_Vector* Y, N_Vector* Z)
{
	/* invalid number of vectors */
	if (nvec < 1) return(-1);

	/* should have called N_VLinearSum */
	if (nvec == 1) {
		N_VLinearSum_Eigen(a[0], x, ONE, Y[0], Z[0]);
		return(0);
	}

	/*
	* Y[i][j] += a[i] * x[j]
	*/
	if (Y == Z) {
		for (int i = 0; i < nvec; i++) {
			EIGV(Y[i]) += a[i] * EIGV(x);
		}
		return(0);
	}

	/*
	* Z[i][j] = Y[i][j] + a[i] * x[j]
	*/
	for (int i = 0; i < nvec; i++) {
		EIGV(Z[i]) = a[i] * EIGV(x) + EIGV(Y[i]);
	}
	return(0);
}


int N_VDotProdMulti_Eigen(int nvec, N_Vector x, N_Vector* Y, realtype* dotprods)
{
	/* invalid number of vectors */
	if (nvec < 1) return(-1);

	/* should have called N_VDotProd */
	if (nvec == 1) {
		dotprods[0] = N_VDotProd_Eigen(x, Y[0]);
		return(0);
	}

	/* compute multiple dot products */
	for (int i = 0; i < nvec; i++) {
		dotprods[i] = EIGV(x).dot(EIGV(Y[i]));
	}

	return(0);
}


/*
 * -----------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------
 */

int N_VLinearSumVectorArray_Eigen(int nvec,
                                   realtype a, N_Vector* X,
                                   realtype b, N_Vector* Y,
                                   N_Vector* Z)
{
	/* invalid number of vectors */
	if (nvec < 1) return(-1);

	/* should have called N_VLinearSum */
	if (nvec == 1) {
		N_VLinearSum_Eigen(a, X[0], b, Y[0], Z[0]);
		return(0);
	}

	/* BLAS usage: axpy y <- ax+y */
	if ((b == ONE) && (Z == Y))
		return(VaxpyVectorArray_Eigen(nvec, a, X, Y));

	/* BLAS usage: axpy x <- by+x */
	if ((a == ONE) && (Z == X))
		return(VaxpyVectorArray_Eigen(nvec, b, Y, X));

	/* Case: a == b == 1.0 */
	if ((a == ONE) && (b == ONE))
		return(VSumVectorArray_Eigen(nvec, X, Y, Z));

	/* Cases:                    */
	/*   (1) a == 1.0, b = -1.0, */
	/*   (2) a == -1.0, b == 1.0 */
	booleantype  test;
	if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE))) {
		N_Vector* V1 = test ? Y : X;
		N_Vector* V2 = test ? X : Y;
		return(VDiffVectorArray_Eigen(nvec, V2, V1, Z));
	}

	/* Cases:                                                  */
	/*   (1) a == 1.0, b == other or 0.0,                      */
	/*   (2) a == other or 0.0, b == 1.0                       */
	/* if a or b is 0.0, then user should have called N_VScale */
	if ((test = (a == ONE)) || (b == ONE)) {
		realtype c = test ? b : a;
		N_Vector* V1 = test ? Y : X;
		N_Vector* V2 = test ? X : Y;
		return(VLin1VectorArray_Eigen(nvec, c, V1, V2, Z));
	}

	/* Cases:                     */
	/*   (1) a == -1.0, b != 1.0, */
	/*   (2) a != 1.0, b == -1.0  */
	if ((test = (a == -ONE)) || (b == -ONE)) {
		realtype c = test ? b : a;
		N_Vector* V1 = test ? Y : X;
		N_Vector* V2 = test ? X : Y;
		return(VLin2VectorArray_Eigen(nvec, c, V1, V2, Z));
	}

	/* Case: a == b                                                         */
	/* catches case both a and b are 0.0 - user should have called N_VConst */
	if (a == b)
		return(VScaleSumVectorArray_Eigen(nvec, a, X, Y, Z));

	/* Case: a == -b */
	if (a == -b)
		return(VScaleDiffVectorArray_Eigen(nvec, a, X, Y, Z));

	/* Do all cases not handled above:                               */
	/*   (1) a == other, b == 0.0 - user should have called N_VScale */
	/*   (2) a == 0.0, b == other - user should have called N_VScale */
	/*   (3) a,b == other, a !=b, a != -b                            */

	/* compute linear sum for each vector pair in vector arrays */
	for (int i = 0; i < nvec; i++) {
		EIGV(Z[i]) = a * EIGV(X[i]) + b * EIGV(Y[i]);
	}

	return(0);
}

int N_VScaleVectorArray_Eigen(int nvec, realtype* c, N_Vector* X, N_Vector* Z)
{
	/* invalid number of vectors */
	if (nvec < 1) return(-1);

	/* should have called N_VScale */
	if (nvec == 1) {
		N_VScale_Eigen(c[0], X[0], Z[0]);
		return(0);
	}

	/*
	* X[i] *= c[i]
	*/
	if (X == Z) {
		for (int i = 0; i < nvec; i++) {
			EIGV(X[i]) *= c[i];
		}
		return(0);
	}

	/*
	* Z[i] = c[i] * X[i]
	*/
	for (int i = 0; i < nvec; i++) {
		EIGV(Z[i]).noalias() = c[i] * EIGV(X[i]);
	}
	return(0);
}

int N_VConstVectorArray_Eigen(int nvec, realtype c, N_Vector* Z)
{
	/* invalid number of vectors */
	if (nvec < 1) return(-1);

	/* should have called N_VConst */
	if (nvec == 1) {
		N_VConst_Eigen(c, Z[0]);
		return(0);
	}

	/* set each vector in the vector array to a constant */
	for (int i = 0; i < nvec; i++) {
		EIGV(Z[i]).setConstant(c);
	}

	return(0);
}

int N_VWrmsNormVectorArray_Eigen(int nvec, N_Vector* X, N_Vector* W,
                                  realtype* nrm)
{
	/* invalid number of vectors */
	if (nvec < 1) return(-1);

	/* should have called N_VWrmsNorm */
	if (nvec == 1) {
		nrm[0] = N_VWrmsNorm_Eigen(X[0], W[0]);
		return(0);
	}

	/* compute the WRMS norm for each vector in the vector array */
	for (int i = 0; i < nvec; i++) {
		nrm[i] = (EIGV(X[i]) * EIGV(W[i])).norm();
	}

	return(0);
}


int N_VWrmsNormMaskVectorArray_Eigen(int nvec, N_Vector* X, N_Vector* W,
                                      N_Vector id, realtype* nrm)
{
	/* invalid number of vectors */
	if (nvec < 1) return(-1);

	/* should have called N_VWrmsNorm */
	if (nvec == 1) {
		nrm[0] = N_VWrmsNormMask_Eigen(X[0], W[0], id);
		return(0);
	}

	/* compute the WRMS norm for each vector in the vector array */
	for (int i = 0; i < nvec; i++) {
		nrm[i] = ZERO;
		for (int j = 0; j < EIGV(X[i]).size(); j++) {
			if (EIGV(id)[j] > ZERO) {
				realtype p = EIGV(X[i])[j] * EIGV(W[i])[j];
				nrm[i] += p * p;
			}
		}
		nrm[i] = SUNRsqrt(nrm[i] / EIGV(X[i]).size());
	}

	return(0);
}

int N_VScaleAddMultiVectorArray_Eigen(int nvec, int nsum, realtype* a,
                                        N_Vector* X, N_Vector** Y, N_Vector** Z)
{
	/* invalid number of vectors */
	if (nvec < 1) return(-1);
	if (nsum < 1) return(-1);

	/* ---------------------------
	* Special cases for nvec == 1
	* --------------------------- */

	if (nvec == 1) {
		/* should have called N_VLinearSum */
		if (nsum == 1) {
			N_VLinearSum_Eigen(a[0], X[0], ONE, Y[0][0], Z[0][0]);
			return(0);
		}

		/* should have called N_VScaleAddMulti */
		N_Vector* YY = (N_Vector*) malloc(nsum * sizeof(N_Vector));
		N_Vector* ZZ = (N_Vector*) malloc(nsum * sizeof(N_Vector));

		for (int j = 0; j < nsum; j++) {
			YY[j] = Y[j][0];
			ZZ[j] = Z[j][0];
		}

		int retval = N_VScaleAddMulti_Eigen(nsum, a, X[0], YY, ZZ);
		free(YY);
		free(ZZ);
		return(retval);
	}

	/* --------------------------
	* Special cases for nvec > 1
	* -------------------------- */

	/* should have called N_VLinearSumVectorArray */
	if (nsum == 1) {
		return N_VLinearSumVectorArray_Eigen(nvec, a[0], X, ONE, Y[0], Z[0]);
	}

	/* ----------------------------
	* Compute multiple linear sums
	* ---------------------------- */

	/*
	* Y[i][j] += a[i] * x[j]
	*/
	if (Y == Z) {
		for (int i = 0; i < nvec; i++) {
			for (int j = 0; j < nsum; j++) {
				EIGV(Y[j][i]) += a[j] * EIGV(X[i]);
			}
		}
		return(0);
	}

	/*
	* Z[i][j] = Y[i][j] + a[i] * x[j]
	*/
	for (int i = 0; i < nvec; i++) {
		for (int j = 0; j < nsum; j++) {
			EIGV(Z[j][i]) = a[j] * EIGV(X[i]) + EIGV(Y[j][i]);
		}
	}
	return(0);
}

int N_VLinearCombinationVectorArray_Eigen(int nvec, int nsum, realtype* c,
                                           N_Vector** X, N_Vector* Z)
{
	/* invalid number of vectors */
	if (nvec < 1) return(-1);
	if (nsum < 1) return(-1);

	/* ---------------------------
	* Special cases for nvec == 1
	* --------------------------- */

	if (nvec == 1) {
		/* should have called N_VScale */
		if (nsum == 1) {
			N_VScale_Eigen(c[0], X[0][0], Z[0]);
			return(0);
		}

		/* should have called N_VLinearSum */
		if (nsum == 2) {
			N_VLinearSum_Eigen(c[0], X[0][0], c[1], X[1][0], Z[0]);
			return(0);
		}

		/* should have called N_VLinearCombination */
		N_Vector* Y = (N_Vector*) malloc(nsum * sizeof(N_Vector));
		for (int i = 0; i < nsum; i++) {
			Y[i] = X[i][0];
		}
		int retval = N_VLinearCombination_Eigen(nsum, c, Y, Z[0]);
		free(Y);
		return(retval);
	}

	/* --------------------------
	* Special cases for nvec > 1
	* -------------------------- */

	/* should have called N_VScaleVectorArray */
	if (nsum == 1) {
		realtype* ctmp = (realtype*)malloc(nvec * sizeof(realtype));
		for (int j = 0; j < nvec; j++) {
			ctmp[j] = c[0];
		}
		int retval = N_VScaleVectorArray_Eigen(nvec, ctmp, X[0], Z);
		free(ctmp);
		return(retval);
	}

	/* should have called N_VLinearSumVectorArray */
	if (nsum == 2) {
		return N_VLinearSumVectorArray_Eigen(nvec, c[0], X[0], c[1], X[1], Z);
	}

	/* --------------------------
	* Compute linear combination
	* -------------------------- */

	/*
	* X[0][j] += c[i]*X[i][j], i = 1,...,nvec-1
	*/
	if ((X[0] == Z) && (c[0] == ONE)) {
		for (int j = 0; j < nvec; j++) {
			for (int i = 1; i < nsum; i++) {
				EIGV(Z[j]) += c[i] * EIGV(X[i][j]);
			}
		}
		return(0);
	}

	/*
	* X[0][j] = c[0] * X[0][j] + sum{ c[i] * X[i][j] }, i = 1,...,nvec-1
	*/
	if (X[0] == Z) {
		for (int j = 0; j < nvec; j++) {
			EIGV(Z[j]) *= c[0];
			for (int i = 1; i < nsum; i++) {
				EIGV(Z[j]) += c[i] * EIGV(X[i][j]);
			}
		}
		return(0);
	}

	/*
	* Z[j] = sum{ c[i] * X[i][j] }, i = 0,...,nvec-1
	*/
	for (int j = 0; j < nvec; j++) {
		EIGV(Z[j]) = c[0] * EIGV(X[0][j]);
		for (int i = 1; i < nsum; i++) {
			EIGV(Z[j]) += c[i] * EIGV(X[i][j]);
		}
	}
	return(0);
}

/*
 * -----------------------------------------------------------------
 * private functions for special cases of vector operations
 * -----------------------------------------------------------------
 */

static void VScaleSum_Eigen(realtype c, N_Vector x, N_Vector y, N_Vector z)
{
	EIGV(z) = c * (EIGV(x) + EIGV(y));
}

static void VScaleDiff_Eigen(realtype c, N_Vector x, N_Vector y, N_Vector z)
{
	EIGV(z) = c * (EIGV(x) - EIGV(y));
}

static void VLin1_Eigen(realtype a, N_Vector x, N_Vector y, N_Vector z)
{
	EIGV(z) = a * EIGV(x) + EIGV(y);
}

static void VLin2_Eigen(realtype a, N_Vector x, N_Vector y, N_Vector z)
{
	EIGV(z) = a * EIGV(x) - EIGV(y);
}

/*
 * -----------------------------------------------------------------
 * private functions for special cases of vector array operations
 * -----------------------------------------------------------------
 */

static int VSumVectorArray_Eigen(int nvec, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
	for (int i = 0; i < nvec; i++) {
		EIGV(Z[i]) = EIGV(X[i]) + EIGV(Y[i]);
	}
	return 0;
}

static int VDiffVectorArray_Eigen(int nvec, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
	for (int i = 0; i < nvec; i++) {
		EIGV(Z[i]) = EIGV(X[i]) - EIGV(Y[i]);
	}
	return 0;
}

static int VScaleSumVectorArray_Eigen(int nvec, realtype c, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
	for (int i = 0; i < nvec; i++) {
		EIGV(Z[i]) = c * (EIGV(X[i]) + EIGV(Y[i]));
	}
	return 0;
}

static int VScaleDiffVectorArray_Eigen(int nvec, realtype c, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
	for (int i = 0; i < nvec; i++) {
		EIGV(Z[i]) = c * (EIGV(X[i]) - EIGV(Y[i]));
	}
	return 0;
}

static int VLin1VectorArray_Eigen(int nvec, realtype a, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
	for (int i = 0; i < nvec; i++) {
		EIGV(Z[i]) = a * EIGV(X[i]) + EIGV(Y[i]);
	}
	return 0;
}

static int VLin2VectorArray_Eigen(int nvec, realtype a, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
	for (int i = 0; i < nvec; i++) {
		EIGV(Z[i]) = a * EIGV(X[i]) - EIGV(Y[i]);
	}
	return 0;
}

static int VaxpyVectorArray_Eigen(int nvec, realtype a, N_Vector* X, N_Vector* Y)
{
	if (a == ONE) {
		for (int i = 0; i < nvec; i++) {
			EIGV(Y[i]) += EIGV(X[i]);
		}
	} else if (a == -ONE) {
		for (int i = 0; i < nvec; i++) {
			EIGV(Y[i]) -= EIGV(X[i]);
		}
	} else {
		for (int i = 0; i < nvec; i++) {
			EIGV(Y[i]) += a * EIGV(X[i]);
		}
	}
	return 0;
}

int N_VEnableFusedOps_Eigen(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  if (tf) {
    /* enable all fused vector operations */
    v->ops->nvlinearcombination = N_VLinearCombination_Eigen;
    v->ops->nvscaleaddmulti     = N_VScaleAddMulti_Eigen;
    v->ops->nvdotprodmulti      = N_VDotProdMulti_Eigen;
    /* enable all vector array operations */
    v->ops->nvlinearsumvectorarray         = N_VLinearSumVectorArray_Eigen;
    v->ops->nvscalevectorarray             = N_VScaleVectorArray_Eigen;
    v->ops->nvconstvectorarray             = N_VConstVectorArray_Eigen;
    v->ops->nvwrmsnormvectorarray          = N_VWrmsNormVectorArray_Eigen;
    v->ops->nvwrmsnormmaskvectorarray      = N_VWrmsNormMaskVectorArray_Eigen;
    v->ops->nvscaleaddmultivectorarray     = N_VScaleAddMultiVectorArray_Eigen;
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Eigen;
  } else {
    /* disable all fused vector operations */
    v->ops->nvlinearcombination = NULL;
    v->ops->nvscaleaddmulti     = NULL;
    v->ops->nvdotprodmulti      = NULL;
    /* disable all vector array operations */
    v->ops->nvlinearsumvectorarray         = NULL;
    v->ops->nvscalevectorarray             = NULL;
    v->ops->nvconstvectorarray             = NULL;
    v->ops->nvwrmsnormvectorarray          = NULL;
    v->ops->nvwrmsnormmaskvectorarray      = NULL;
    v->ops->nvscaleaddmultivectorarray     = NULL;
    v->ops->nvlinearcombinationvectorarray = NULL;
  }

  /* return success */
  return(0);
}

int N_VEnableLinearCombination_Eigen(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombination = N_VLinearCombination_Eigen;
  else
    v->ops->nvlinearcombination = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMulti_Eigen(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmulti = N_VScaleAddMulti_Eigen;
  else
    v->ops->nvscaleaddmulti = NULL;

  /* return success */
  return(0);
}

int N_VEnableDotProdMulti_Eigen(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvdotprodmulti = N_VDotProdMulti_Eigen;
  else
    v->ops->nvdotprodmulti = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearSumVectorArray_Eigen(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearsumvectorarray = N_VLinearSumVectorArray_Eigen;
  else
    v->ops->nvlinearsumvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleVectorArray_Eigen(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscalevectorarray = N_VScaleVectorArray_Eigen;
  else
    v->ops->nvscalevectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableConstVectorArray_Eigen(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvconstvectorarray = N_VConstVectorArray_Eigen;
  else
    v->ops->nvconstvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableWrmsNormVectorArray_Eigen(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvwrmsnormvectorarray = N_VWrmsNormVectorArray_Eigen;
  else
    v->ops->nvwrmsnormvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableWrmsNormMaskVectorArray_Eigen(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvwrmsnormmaskvectorarray = N_VWrmsNormMaskVectorArray_Eigen;
  else
    v->ops->nvwrmsnormmaskvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMultiVectorArray_Eigen(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmultivectorarray = N_VScaleAddMultiVectorArray_Eigen;
  else
    v->ops->nvscaleaddmultivectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearCombinationVectorArray_Eigen(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Eigen;
  else
    v->ops->nvlinearcombinationvectorarray = NULL;

  /* return success */
  return(0);
}
