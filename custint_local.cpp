/*

	simple interval struct with end-points double
	based on constant upward rounding, no
	rounding mode switch necessary during operations

	Marc Meidlinger
	May-July 2020

	NO WARRANTY OF CORRECT RESULTS (although I did my best to
	ensure correctness of the implementation)

	CAUTION: Every arithmetics routine cust^NNN_ZAB()
		 ASSUMES rounding mode is set to UPWARDS

		 It is at the discretion of the user to
		 ensure this by defining a struct of
		 type CustRoundingUpwards

		(except for sqrt-functions uses rounding mode switch)

	based on the article:

	Interval arithmetic with fixed rounding modw
	by S. Rump et al., 2016

*/

#ifndef _CUSTOMIZED_INTERVALC
#define _CUSTOMIZED_INTERVALC

#include "custint.h"
#include "stdio.h"
#include "stdlib.h"
#include "stdint.h"
#include "math.h"
#include "fenv.h"
#include "float.h"

using namespace std;


// forward

// general functions
void terminate_program(void);
inline double minimumDouble(const double,const double,const double,const double);
inline double maximumDouble(const double,const double,const double,const double);

// CustInterval operations
int32_t custAdd_ZAB(CustInterval&,CustInterval&,CustInterval&);
int32_t custSub_ZAB(CustInterval&,CustInterval&,CustInterval&);
int32_t custMul_ZAB(CustInterval&,CustInterval&,CustInterval&);
int32_t custDiv_ZAB(CustInterval&,CustInterval&,CustInterval&);
int32_t custSqrt_ZA_switch(CustInterval&,CustInterval&);
int32_t custSqrt_ZA(CustInterval&,CustInterval&);
int32_t custPow2_ZA(CustInterval&,CustInterval&);
void custPrint(FILE*,CustInterval&);
int32_t custMidpoint_ZA(CustInterval&,CustInterval&);

// CustComplex operations
int32_t custCplxAdd_ZAB(CustComplex&,CustComplex&,CustComplex&);
int32_t custCplxSub_ZAB(CustComplex&,CustComplex&,CustComplex&);
int32_t custCplxMul_ZAB(CustComplex&,CustComplex&,CustComplex&);
int32_t custCplxMul_ZArealint64(CustComplex&,CustComplex&,const int64_t);
int32_t custCplxDiv_ZAB(CustComplex&,CustComplex&,CustComplex&);
int32_t custCplxNormQ_ZA(CustInterval&,CustComplex&);
int32_t custCplxSqrt_ZA(CustComplex&,CustComplex&);
int32_t custCplxMidpoint_ZA(CustComplex&,CustComplex&);

// custMatrix routines
int32_t custMatrixAdd_ZAB(CustIntervalMatrix2x2&,CustIntervalMatrix2x2&,CustIntervalMatrix2x2&);
int32_t custMatrixSub_ZAB(CustIntervalMatrix2x2&,CustIntervalMatrix2x2&,CustIntervalMatrix2x2&);
int32_t custMatrixMul_ZAB(CustIntervalMatrix2x2&,CustIntervalMatrix2x2&,CustIntervalMatrix2x2&);
int32_t custMatrixCplxMul_ZMA(CustComplex&,CustIntervalMatrix2x2&,CustComplex&);
int32_t custMatrixInverse_ZA(CustIntervalMatrix2x2&,CustIntervalMatrix2x2&);


// general function

void terminate_program(void) {
	fprintf(stderr,"Program terminates.\n");
	exit(99);
}

double minimumDouble(
	const double a,const double b,
	const double c,const double d
) {
	double m=a;
	if (b < m) m=b;
	if (c < m) m=c;
	if (d < m) m=d;

	return m;
}

double maximumDouble(
	const double a,const double b,
	const double c,const double d
) {
	double m=a;
	if (b > m) m=b;
	if (c > m) m=c;
	if (d > m) m=d;

	return m;
}


// struct CustRoundingUpwards
// the constructor retrieves the current rounding mode and
// saves it, then sets to upward infinity
// the destructor resets to the old rounding mode

CustRoundingUpwards::CustRoundingUpwards() {
	mode=fegetround();
	if (mode<0) {
		fprintf(stderr,"CustRoundingUpwards::construct. Cannot get current rounding mode.\n");
		terminate_program();
	}
	if (fesetround(FE_UPWARD) != 0) {
		fprintf(stderr,"CustRoundingUpwards::constructor. Cannot set rounding mode to upwards.\n");
		terminate_program();
	}

	//printf("\nCustRoundingUpwards::set %i\n",fegetround());
}

CustRoundingUpwards::~CustRoundingUpwards() {
	if (fesetround(mode) != 0) {
		fprintf(stderr,"CustRoundingUpwards::destructor. Cannot reset rounding mode to %i.\n",mode);
		terminate_program();
	}
}

// CustRoundingDownards
CustRoundingDownwards::CustRoundingDownwards() {
	mode=fegetround();
	if (mode<0) {
		fprintf(stderr,"CustRoundingDownwards::construct. Cannot get current rounding mode.\n");
		terminate_program();
	}
	if (fesetround(FE_DOWNWARD) != 0) {
		fprintf(stderr,"CustRoundingDownwards::constructor. Cannot set rounding mode to upwards.\n");
		terminate_program();
	}
}

CustRoundingDownwards::~CustRoundingDownwards() {
	if (fesetround(mode) != 0) {
		fprintf(stderr,"CustRoundingDownwards::destructor. Cannot reset rounding mode to %i.\n",mode);
		terminate_program();
	}
}

// CustInterval
void CustInterval::copyFrom( CustInterval& A ) {
    left=A.left;
    right=A.right;
}

void CustInterval::ausgabe(FILE* f) {
    fprintf(f,"[%.20lg..%.20lg]",left,right);
}

CustInterval::CustInterval() {
	// do nothing for speed reasons
}

CustInterval::CustInterval(const double a) {
	left=right=a;
	int32_t fpa=fpclassify(a);
	if (
		(fpa != FP_NORMAL) &&
		(fpa != FP_ZERO)
	) {
		fprintf(stderr,"CustInterval::constructor. Error, only normal floating-points (and 0) are allowed.\n");
		terminate_program();
	}
}

CustInterval::CustInterval(const double a,const double b) {
	left=a;
	right=b;
	int32_t fpa=fpclassify(a);
	int32_t fpb=fpclassify(b);
	if (
		(
			(fpa != FP_NORMAL) &&
			(fpa != FP_ZERO)
		) ||
		(
			(fpb != FP_NORMAL) &&
			(fpb != FP_ZERO)
		)
	) {
		fprintf(stderr,"CustInterval::constructor. Error, only normal floating-points (and 0) are allowed.\n");
		fprintf(stderr,": %.20lg..%.20lg\n",a,b);
		terminate_program();
	}

	if (a > b) {
		fprintf(stderr,"CustInterval::constructor. Not an interval. End points in wrong order.\n");
		terminate_program();
	}
}


// arithmetics

// all routines assume that the given parameters
// are normal or zero (i.e. no Inf, NaN or subnormals)
// (except division, that checks for the dividend
// not containing zero)
// resulting intervals are always normals or zero at
// end points (returning value 0).
// If not possible, the return value is -1
// and should be handled by the calling function appropriately

// a subnormal resulting interval end point is
// moved towards a circumferencing next normal (or zero)
#define IAVALIDATE_LEFT(VAR) \
{\
	switch (fpclassify(VAR)) {\
		case FP_INFINITE:\
		case FP_NAN: return -1; /* error */\
		case FP_SUBNORMAL: {\
			if (VAR > 0.0) {\
				/* left end, positive subnormal */\
				/* move to 0 */\
				VAR=0.0;\
			} else if (VAR < 0.0) {\
				/* left end, negative subnormal */\
				/* move to -FMIN */\
				VAR=-DBL_MIN;\
			}\
		}\
		default: {\
			break;\
		}\
	}\
}

#define IAVALIDATE_RIGHT(VAR) \
{\
	switch (fpclassify(VAR)) {\
		case FP_INFINITE:\
		case FP_NAN: return -1; /* error */\
		case FP_SUBNORMAL: {\
			if (VAR > 0.0) {\
				/* right end, positive subnormal */\
				/* move to smallest double */\
				VAR=DBL_MIN;\
			} else if (VAR < 0.0) {\
				/* right end, negative subnormal */\
				/* move to 0 */\
				VAR=0.0;\
			}\
		}\
		default: {\
			break;\
		}\
	}\
}

int32_t custAdd_ZAB(
	CustInterval& res,
	CustInterval& A,
	CustInterval& B
) {
	// /////////////////////////////////////////
	// rounding mode is assumed UPWARD
	// /////////////////////////////////////////

	// it is assumed that res is not the same object as A or B

	// left end: simulates downwards rounding by
	// consecutive negation (see articel by S.Rump)
	double a1=A.left; a1=-a1; // no rounding occurs
	double b1=B.left; b1=-b1; // no rounding occurs

	// ROUNDING upwards
	res.left=a1+b1;

	res.left=-res.left; // no rounding
	// now the result is DOWNWARD-rounded A.left+B.left

	// right end-point: uses upward rounding
	res.right=A.right+B.right;

	// if result is subnormal => enlarge interval or
	// return error value -1 if infinity or NaN
	IAVALIDATE_LEFT(res.left);
	IAVALIDATE_RIGHT(res.right);

	// everything worked, valid interval in variable res
	return 0;
}

int32_t custSub_ZAB(
	CustInterval& res,
	CustInterval& A,
	CustInterval& B
) {
	// /////////////////////////////////////////
	// rounding mode is assumed UPWARD
	// /////////////////////////////////////////

	// it is assumed that res is not the same object as A or B

	// [aL..aR] - [bL..bR] = [aL-bR..aR-bL]

	// left end: simulates downwards rounding by
	// consecutive negation (see articel by S.Rump)
	// aL -DOWN bR = -( (-aL) UP- (-bR) )
	double a1=A.left;  a1=-a1; // no rounding occurs
	double b1=B.right; b1=-b1; // no rounding occurs

	// ROUNDING upwards
	res.left=a1-b1;

	res.left=-res.left; // no rounding
	// now the result is DOWNWARD-rounded A.left-B.right

	// right end-point: uses upward rounding
	res.right=A.right-B.left;

	// if result is subnormal => enlarge interval or
	// return error value -1 if infinity or NaN
	IAVALIDATE_LEFT(res.left);
	IAVALIDATE_RIGHT(res.right);

	// everything worked, valid interval in variable res
	return 0;

}

int32_t custMul_ZAB(
	CustInterval& res,
	CustInterval& A,
	CustInterval& B
) {
	// /////////////////////////////////////////
	// rounding mode is assumed UPWARD
	// /////////////////////////////////////////

	// it is assumed that res is not the same object as A or B

	// could be implemented using sign checks

	// downward-mul: -( upward(a*(-b)) )

	/*
	// the following method is about 10% SLOWER
	// than the general one, probably because of disrupting
	// the flow of instructions in the CPU pipeline
	// the double number type is fast enough as it is
	// saving some multiplications doesn't affect
	// the running time positively
	// TODO: for bigfpiacplx e.g. it could
	// be faster
	// check A, B fully positive or negative
	// then the result can be obtained by a total of
	// two multiplications instead of 8

	if (A.left > 0.0) {
		// A is positive
		if (B.left > 0.0) {
			// b positive
			// 0 : A,B
			// result: Aleft *DOWN Bleft .. Aright *UP Bright

			// forced execution order to
			// obtain the desired ROUNDING effect
			double w1=B.left; w1=-w1; // no rounding

			// ROUNDING upward
			w1 *= A.left;

			// no rounding
			w1=-w1;
			// result is DOWNWARD-rounded A.left*B.left;
			res.left = w1;

			res.right = A.right * B.right; // UPWARD rpunding

			IAVALIDATE_LEFT(res.left);
			IAVALIDATE_RIGHT(res.right);

			return 0;
		} else if (B.right < 0.0) {
			// b negative
			// B : 0 : A
			// result: Aright *DOWN Bleft .. Aleft *UP Bright

			double w1=B.left; w1=-w1; w1 *= A.right;  w1=-w1;
			res.left = w1;

			res.right = A.left * B.right; // UPWARD

			IAVALIDATE_LEFT(res.left);
			IAVALIDATE_RIGHT(res.right);

			return 0;
		}
	} else if (A.right < 0.0) {
		// A fully negative
		if (B.left > 0.0) {
			// B positive
			// A : 0 : B
			// result: Aleft *DOWN Bright . Aright *UP Bleft

			double w1=B.right; w1=-w1; w1 *= A.left;  w1=-w1;
			res.left = w1;

			res.right = A.right * B.left; // UPWARD

			IAVALIDATE_LEFT(res.left);
			IAVALIDATE_RIGHT(res.right);

			return 0;
		} else if (B.right < 0.0) {
			// B negative
			// A,B : 0
			// result: Aright *DOWN Bright .. Aleft *UP Bleft

			double w1=B.right; w1=-w1; w1 *= A.right;  w1=-w1;
			res.left = w1;

			res.right = A.left * B.left; // UPWARD

			IAVALIDATE_LEFT(res.left);
			IAVALIDATE_RIGHT(res.right);

			return 0;
		} // B
	} // A
	*/


	// general method
	// (with double type the fastest method, better
	// than reducing number of multiplications
	// by checking where the intervals lie

	// left end

	// forced execution order to
	// obtain the desired ROUNDING effect
	double w1=B.left; w1=-w1; // no rounding

	// ROUNDING upward
	w1 *= A.left;

	// no rounding
	w1=-w1;
	// result is DOWNWARD-rounded A.left*B.left;

	// and the others
	double w2=B.right; w2=-w2; w2 *= A.left;  w2=-w2;
	double w3=B.left;  w3=-w3; w3 *= A.right; w3=-w3;
	double w4=B.right; w4=-w4; w4 *= A.right; w4=-w4;
	res.left=minimumDouble(w1,w2,w3,w4);

	// right end: only binary operation, so execution
	// order is clear and temporary results can be passed directly
	// to maximumDouble
	res.right=maximumDouble(
		A.left*B.left,
		A.left*B.right,
		A.right*B.left,
		A.right*B.right
	);

	// if resulting end point(s) are subnormal => adjust
	// if inf,nan => return error value -1
	IAVALIDATE_LEFT(res.left);
	IAVALIDATE_RIGHT(res.right);

	// everything worked, valid interval in variable res
	return 0;
}

int32_t custDiv_ZAB(
	CustInterval& res,
	CustInterval& A,
	CustInterval& B
) {
	// /////////////////////////////////////////
	// rounding mode is assumed UPWARD
	// /////////////////////////////////////////

	// it is assumed that res is not the same object as A or B

	// [a0..a1] / [b0..b1] = [a0..b0] * [1/b1..1/b0]
	// b must not contain zero

	// end-point test:
	// here also a check for inf/nan is possible without
	// speed-loss as the check for zero at the end-points
	// has to be done anyways
	if (fpclassify(B.left) != FP_NORMAL) return -1; // also checks for zero
	if (fpclassify(B.right)!= FP_NORMAL) return -1; // also checks for zero

	// zero within?
	if ( (B.left < 0.0) && (B.right > 0.0) ) return -1;

	CustInterval binv;
	// left end: downward 1/b.right = -( UPWARD( (-1)/B.right ) )
	binv.left=-1.0; // no rounding

	// ROUNDING upward
	binv.left /= B.right;

	binv.left = -binv.left; // no rounding
	// now binv.left is downward rounded 1/B.right

	// right end point: upward rounding
	binv.right=1.0/B.left;

	// subnormals can be moved to zero or DBL_min
	// i.e. binv can be a superset of the true 1/B
	// as basic interval operations are inclusion monotone,
	// multiplication below is a valid result
	IAVALIDATE_LEFT(binv.left)
	IAVALIDATE_RIGHT(binv.right)

	return custMul_ZAB(res,A,binv);
}

int32_t custSqrt_ZA_switch(
	CustInterval& res,
	CustInterval& A
) {
	// currently uses switching rounding mode
	//printf("CustIA::sqrt: experimental\n");

	{
		CustRoundingUpwards rdup;
		if (A.right < 0.0) return -1; // error
		res.right=sqrt(A.right);

		// implicit destruction
	}

	{
		CustRoundingDownwards rddown;
		if (A.left < 0.0) return -1; // error
		res.left=sqrt(A.left);

		// implicit destruction
	}

	// now rounding mode is set back to up
	if (fegetround() != FE_UPWARD) {
		fprintf(stderr,"custSqrt: error, not upward rounding\n");
		return -1;
	}

	IAVALIDATE_LEFT(res.left);
	IAVALIDATE_RIGHT(res.right);

	return 0;
}

int32_t custSqrt_ZA(
	CustInterval& res,
	CustInterval& A
) {
	fprintf(stderr,"CustIA::sqrt: noswitch experimental\n");
	return -1;

	// uses upward computed sqrt-function
	// and substracts (downard) 1 ulp
	// checks for underflow


	return 0;
}

int32_t custPow_ZAE(
	CustInterval& erg,
	CustInterval& A,
	const int64_t exponent
) {
	if (exponent < 0) return -1; // not possible

	if (exponent == 0) {
		erg.left=erg.right=1.0;
		return 0;
	}

	erg.left=A.left;
	erg.right=A.right;

	if (exponent == 1) {
		return 0;
	}

	int32_t error=0;

	CustInterval tmp;

	for(int32_t i=2;i<=exponent;i++) {
		error += custMul_ZAB(tmp,erg,A);

		erg.left=tmp.left;
		erg.right=tmp.right;
	} // i

	return error;
}

int32_t custPow2_ZA(CustInterval& erg,CustInterval& A) {
	// if straddling zero => result: 0..max(A.left^2,A.right^2)
	// else use multiplication routine

	int32_t error=0;

	if ( (A.left < 0.0) && (A.right > 0.0 ) ) {
		erg.left=0.0;

		// UPWARD rounding
		double l2=A.left*A.left;
		double r2=A.right*A.right;

		if (l2 > r2) erg.right=l2;
		else erg.right=r2;
	} else {
		error += custMul_ZAB(erg,A,A);
	}

	return error;
}

// CustComplex
void CustComplex::copyFrom( CustComplex& A ) {
    re.copyFrom( A.re );
    im.copyFrom( A.im );
};

void CustComplex::ausgabe(FILE* f) {
    re.ausgabe(f);
    fprintf(f,"+i*");
    im.ausgabe(f);
}

CustComplex::CustComplex(const CustInterval& v1,const CustInterval& v2) {
	re.left=v1.left;
	re.right=v1.right;
	im.left=v2.left;
	im.right=v2.right;
}

CustComplex::CustComplex(const double a,const double b) {
	// complex point
	re.left=re.right=a;
	im.left=im.right=b;
}

CustComplex::CustComplex() {
	// no initialization for speed up
}

int32_t custCplxAdd_ZAB(CustComplex& erg,CustComplex& A,CustComplex& B) {
	int32_t error=0;
	// real and imag interval are added separately
	error += custAdd_ZAB(erg.re,A.re,B.re);
	error += custAdd_ZAB(erg.im,A.im,B.im);

	return error;
}

int32_t custCplxSub_ZAB(CustComplex& erg,CustComplex& A,CustComplex& B) {
	int32_t error=0;
	// real and imag interval are subtracted separately
	error += custSub_ZAB(erg.re,A.re,B.re);
	error += custSub_ZAB(erg.im,A.im,B.im);

	return error;
}

int32_t custCplxMul_ZAB(CustComplex& erg,CustComplex& A,CustComplex& B) {
	int32_t error=0;
	// (Are+Aim*i) * (Bre+Bim*i) =
	// (Are*Bre-Aim*Bim + (Are*Bim+Aim*Bre)*i)
	CustInterval arebre,arebim,aimbre,aimbim;

	error += custMul_ZAB(arebre,A.re,B.re);
	error += custMul_ZAB(arebim,A.re,B.im);
	error += custMul_ZAB(aimbre,A.im,B.re);
	error += custMul_ZAB(aimbim,A.im,B.im);

	error += custSub_ZAB(erg.re,arebre,aimbim);
	error += custAdd_ZAB(erg.im,arebim,aimbre);

	return error;
}

int32_t custCplxDiv_ZAB(CustComplex& erg,CustComplex& A,CustComplex& B) {
	int32_t error=0;
	// (a+b*i)/(c+d*i) =: (e+f*i)
	// by maxima:
	// realpart: (b*d+a*c)/(d^2+c^2)
	// imagpart: (b*c-a*d)/(d^2+c^2)
	CustInterval d2c2,d2,c2;
	error += custPow2_ZA(c2,B.re); // interval does not straddle 0, but may end-point contain it
	error += custPow2_ZA(d2,B.im);
	error += custAdd_ZAB(d2c2,d2,c2);

	// realpart
	CustInterval bd,ac,bdac;
	error += custMul_ZAB(bd,A.im,B.im);
	error += custMul_ZAB(ac,A.re,B.re);
	error += custAdd_ZAB(bdac,bd,ac);
	error += custDiv_ZAB(erg.re,bdac,d2c2);

	// imagpart: (b*c-a*d)/d2c2
	CustInterval bc,ad,bcad;
	error += custMul_ZAB(bc,A.im,B.re);
	error += custMul_ZAB(ad,A.re,B.im);
	error += custSub_ZAB(bcad,bc,ad);
	error += custDiv_ZAB(erg.im,bcad,d2c2);

	return error;
}

int32_t custCplxNormQ_ZA(CustInterval& normq,CustComplex& A) {
	// coputes the SQUARE of the norm [a..b]^2+[c..d]^2
	// from the input complex numer [a..b]x[c..d]
	int32_t error=0;

	CustInterval realq,imagq;
	// ATTN: re,im CAN contain zero, so
	// re*re is NOT applicable as it might result
	// in an norm-interval with negative endpoint
	// one has to use a power-definition
	error += custPow2_ZA(realq,A.re);
	error += custPow2_ZA(imagq,A.im);
	error += custAdd_ZAB(normq,realq,imagq);

	return error;
}

#define OUTIA(TEXT,WW) \
{\
	printf("%s; %.20lg..%.20lg\n",TEXT,WW.left,WW.right);\
}

// using princiapl root formulka
int32_t custCplxSqrt_ZA_principal(
	CustComplex& erg,
	CustComplex& val
) {
	// jump back with error as soon as possible

	if (
		(fpclassify(val.im.left) == FP_ZERO) &&
		(fpclassify(val.im.right) == FP_ZERO)
	) {
		// pure real root
		erg.im.left=erg.im.right=0.0;

		return custSqrt_ZA_switch(erg.re,val.re);
	}
	// using principal root formula from
	// wikipedia:

	// sqrt(a+b*i) =: e+f*i
	// e=sqrt( 0.5*( sqrt(a^2+b^2) + a) )
	// f=sign(b)*sqrt( 0.5*(sqrt(a^2+b^2) - a) )

	// so if val.im contains 0 at an end or within
	// return with error. Those have to be computed
	// using a different IA library
	// currently not applicable
	#define CUSTCONTAINSZERO(II) \
	(\
		(II.left <= 0.0) && (II.right >= 0.0)\
	)

	if (CUSTCONTAINSZERO(val.im) > 0) {
		return -1;
	}



	CustInterval a2,b2,a2b2,sqrta2b2;

	// direct return to reduce number of
	// operations, as sqrt is then compüted by
	// a different IA library, e.g. kv
	if (custMul_ZAB(a2,val.re,val.re) != 0) return -1;
	if (custMul_ZAB(b2,val.im,val.im) != 0) return -1;
	if (custAdd_ZAB(a2b2,a2,b2) != 0) return -1;
	if (custSqrt_ZA_switch(sqrta2b2,a2b2) != 0) return -1;
	CustInterval t1,t2;
	CustInterval cust05(0.5,0.5);
	CustInterval custminus1(-1.0,-1.0);
	if (custSub_ZAB(t1,sqrta2b2,val.re) != 0) return -1;
	if (custMul_ZAB(t2,t1,cust05) != 0) return -1;

	if (val.im.right < 0.0) {
		CustInterval t3;
		if (custSqrt_ZA_switch(t3,t2) != 0) return -1;
		if (custMul_ZAB(erg.im,t3,custminus1) != 0) return -1;;
	} else {
		if (custSqrt_ZA_switch(erg.im,t2) != 0) return -1;
	}

	// real part
	// e=sqrt( t5 )
	CustInterval t4,t5;
	if (custAdd_ZAB(t4,sqrta2b2,val.re) != 0) return -1;
	if (custMul_ZAB(t5,t4,cust05) != 0) return -1;
	if (custSqrt_ZA_switch(erg.re,t5) != 0) return -1;

	return 0;
}

void custPrint(FILE* f,CustInterval& a) {
	fprintf(f,"[%.20lg..%.20lg]",a.left,a.right);
}

/* custIntervalMatrix2x2 */

int32_t custMatrixAdd_ZAB(
	CustIntervalMatrix2x2& erg,
	CustIntervalMatrix2x2& A,
	CustIntervalMatrix2x2& B
) {
	int32_t error=0;

	error += custAdd_ZAB(erg.x1y1,A.x1y1,B.x1y1);
	error += custAdd_ZAB(erg.x2y1,A.x2y1,B.x2y1);
	error += custAdd_ZAB(erg.x1y2,A.x1y2,B.x1y2);
	error += custAdd_ZAB(erg.x2y2,A.x2y2,B.x2y2);

	return error;
}

int32_t custMatrixSub_ZAB(
	CustIntervalMatrix2x2& erg,
	CustIntervalMatrix2x2& A,
	CustIntervalMatrix2x2& B
) {
	int32_t error=0;

	error += custSub_ZAB(erg.x1y1,A.x1y1,B.x1y1);
	error += custSub_ZAB(erg.x2y1,A.x2y1,B.x2y1);
	error += custSub_ZAB(erg.x1y2,A.x1y2,B.x1y2);
	error += custSub_ZAB(erg.x2y2,A.x2y2,B.x2y2);

	return error;
}

int32_t custMatrixMul_ZAB(
	CustIntervalMatrix2x2& erg,
	CustIntervalMatrix2x2& A,
	CustIntervalMatrix2x2& B
) {
	int32_t error=0;

	/*
		(
			A11[a1..a2]	A21[b1..b2]
			A12[c1..c2]	A22[d1..d2]
		)

		*

		(
			B11[e1..e2]	B21[f1..f2]
			B12[g1..g2]	B22[h1..h2]
		)

		=

		(
			ae+bg		af+bh
			ce+dg		cf+dh
		)

		=

		(
			A11B11+A21B12		A11B21+A21B22
			A12B11+A22B12		A12B21+A22B22
		)

	*/

	CustInterval t1,t2;
	error += custMul_ZAB(t1,A.x1y1,B.x1y1);
	error += custMul_ZAB(t2,A.x2y1,B.x1y2);
	error += custAdd_ZAB(erg.x1y1,t1,t2);

	error += custMul_ZAB(t1,A.x1y1,B.x2y1);
	error += custMul_ZAB(t2,A.x2y1,B.x2y2);
	error += custAdd_ZAB(erg.x2y1,t1,t2);

	error += custMul_ZAB(t1,A.x1y2,B.x1y1);
	error += custMul_ZAB(t2,A.x2y2,B.x1y2);
	error += custAdd_ZAB(erg.x1y2,t1,t2);

	error += custMul_ZAB(t1,A.x1y2,B.x2y1);
	error += custMul_ZAB(t2,A.x2y2,B.x2y2);
	error += custAdd_ZAB(erg.x2y2,t1,t2);

	return error;
}

int32_t custMatrixCplxMul_ZMA(
	CustComplex& erg,
	CustIntervalMatrix2x2& M,
	CustComplex& A
) {
	// erg=M*A; matrix of CustInterval to the left

	/*

			M11	M21		Ar		M11*Ar+M21*Ai
					*		=
			M12	M22		Ai		M12*Ar+M22*Ai

	*/

	int32_t error=0;

	CustInterval t1,t2;
	error += custMul_ZAB(t1,M.x1y1,A.re);
	error += custMul_ZAB(t2,M.x2y1,A.im);
	error += custAdd_ZAB(erg.re,t1,t2);

	CustInterval t3,t4;
	error += custMul_ZAB(t3,M.x1y2,A.re);
	error += custMul_ZAB(t4,M.x2y2,A.im);
	error += custAdd_ZAB(erg.im,t3,t4);


	return error;

}

int32_t custMidpoint_ZA(CustInterval& erg,CustInterval& A) {
	// compue via IAmul
	CustInterval tx0,tx1,halb,t1;
	halb.left=halb.right=0.5;

	tx0.left=tx0.right=A.left;
	tx1.left=tx1.right=A.right;

	int32_t error=0;

	error += custAdd_ZAB(t1,tx0,tx1);
	error += custMul_ZAB(erg,halb,t1);

	// erg must be a subset (or equal) to A
	// so if rounding results in the "midpoint"
	// been shifted outside => error

	if (
		(erg.left < A.left) ||
		(erg.right > A.right)
	) error+=-1;

	//printf("\nmidp of "); A.ausgabe(stdout);
	//printf("\nerr%i => ",error);
	//erg.ausgabe(stdout);

	return error;
}

int32_t custCplxMidpoint_ZA(
	CustComplex& erg,
	CustComplex& A
) {
	int32_t error=0;

	error += custMidpoint_ZA(erg.re,A.re);
	error += custMidpoint_ZA(erg.im,A.im);

	return error;
}

int32_t custMatrixInverse_ZA(
	CustIntervalMatrix2x2& erg,
	CustIntervalMatrix2x2& A
) {
	/*
					)-1
		x11		x21	)		1	x22		-x21
					)	=  ---
		x12		x22	)	   det	-x12	x11
					)

		det=x11*x22 - x21*x12

	*/

	int32_t error=0;

	CustInterval det,t1,t2;
	error += custMul_ZAB(t1,A.x1y1,A.x2y2);
	error += custMul_ZAB(t2,A.x2y1,A.x1y2);
	error += custSub_ZAB(det,t1,t2);

	if (error != 0) return -1;

	// 0 in detrminant ?
	if (
		(fpclassify(det.left) == FP_ZERO) ||
		(fpclassify(det.right) == FP_ZERO) ||
		(
			(det.left <= 0.0) &&
			(det.right >= 0.0)
		)
	) {
		return -2; // special return value: not invertible
	}

	CustInterval t3,t4,minus1;
	minus1.left=minus1.right=-1.0;

	error += custDiv_ZAB(erg.x1y1,A.x2y2,det);
	error += custDiv_ZAB(t3,A.x2y1,det);
	error += custMul_ZAB(erg.x2y1,minus1,t3);
	error += custDiv_ZAB(t4,A.x1y2,det);
	error += custMul_ZAB(erg.x1y2,minus1,t4);
	error += custDiv_ZAB(erg.x2y2,A.x1y1,det);

	if (error != 0) return -1;

	return error;
}

// matrix
void CustIntervalMatrix2x2::ausgabe(FILE* f) {
    x1y1.ausgabe(f);
    fprintf(f," / ");
    x2y1.ausgabe(f);
    fprintf(f,"\n");
    x1y2.ausgabe(f);
    fprintf(f," / ");
    x2y2.ausgabe(f);
    fprintf(f,"\n");
}

CustIntervalMatrix2x2::CustIntervalMatrix2x2() {
	// no initialization
}

CustIntervalMatrix2x2::CustIntervalMatrix2x2(
	const double a,const double b,
	const double c,const double d
) {
	x1y1=a;
	x2y1=b;
	x1y2=c;
	x2y2=d;
}

// unreliable helper routines
struct psComplex {
    double re,im;

    psComplex(const double,const double);
    psComplex(void) { } // does nothing
    void copyFrom(psComplex&);

    double norm(void);
    void ausgabe(FILE*);
};

void psComplex::ausgabe(FILE* f) {
    fprintf(f,"%.20lg,%.20lg ",
        re,im);
}

void psComplex::copyFrom( psComplex& A ) {
    re=A.re;
    im=A.im;
}

psComplex::psComplex(const double a,const double b) {
    re=a;
    im=b;
}

double psComplex::norm(void) {
    return sqrt(
        re*re+im*im
    );
}

psComplex csqrt3(psComplex c) {
    // convert to polar coordinates
    double R=c.norm();
    double angle=atan2( c.im,c.re );
    if (angle < 0.0) angle += (2*M_PI);

    // third
    R=exp(log(R)/3.0); // cube root of magnitude
    angle /= 3.0; // one third of the angle

    return psComplex(
        R*cos(angle),
        R*sin(angle)
    );
}

psComplex csqrt2(psComplex c) {
    // convert to polar coordinates
    double R=c.norm();
    double angle=atan2( c.im,c.re );
    if (angle < 0.0) angle += (2*M_PI);

    // third
    R=sqrt(R); // square root of magnitude
    angle *= 0.5; // one half of the angle

    return psComplex(
        R*cos(angle),
        R*sin(angle)
    );
}

psComplex cexp(psComplex z) {
	double ex=exp((double)(z.re));
	return psComplex(
		ex*cos(z.im),
		ex*sin(z.im)
	);
}

void CustInterval::setToZero(void) {
    left=right=0.0;
}

int32_t CustInterval::setStr(const char* str) {
    // returns: 0 => successfully set
    // < 0 => error

    // only digits, at most one . nothing else
    int32_t slen=strlen(str);
    int8_t punkt=0;

    for(int32_t k=0;k<slen;k++) {
        if (
            (str[k] == '-') ||
            (str[k] == '+')
        ) {
            if (k > 0) return -1; // error
            continue;
        }

        if (str[k] == '.') {
            if (punkt > 0) return -1; // error
            punkt=1;
            continue;
        }

        if (
            (str[k] < '0') ||
            (str[k] > '9')
        ) {
            return -1; // error
        }

    } // k

    CustInterval erg;
    erg.setToZero();
    int32_t error=0;
    int32_t endidx=strlen(str)-1,startidx=0;
    int8_t sign=1;
    if (str[0] == '-') {
        sign=-1;
        startidx=1;
    } else
    if (str[0] == '+') {
        sign=1;
        startidx=1;
    }

    // integer part
    CustInterval factor,zehn;
    CustInterval fractfactor;
    zehn.left=zehn.right=10.0;
    factor.left=factor.right=1.0;
    punkt=0;
    endidx=slen-1;

    CustInterval t1,t2;
    for(int32_t k=endidx;k>=startidx;k--) {
        if (str[k] == '.') {
            punkt=1;
            fractfactor.copyFrom(factor);
            continue;
        }

        int32_t w=str[k]-'0';

        CustInterval a;
        a.left=a.right=w; // fits
        error += custMul_ZAB(t1,a,factor);

        t2.copyFrom( erg );
        error += custAdd_ZAB(erg,t2,t1);

        t2.copyFrom( factor );
        error += custMul_ZAB(factor,t2,zehn);

    } // k

    if (punkt > 0) {
        // divide by fractfactor
        t1.copyFrom( erg );
        error += custDiv_ZAB(erg,t1,fractfactor);
    }

    if (sign < 0) {
        t1.copyFrom( erg );
        t2.left=t2.right=-1.0;
        error += custMul_ZAB(erg,t1,t2);
    }

    left=erg.left;
    right=erg.right;

    return error;

}

void CustComplex::setToZero(void) {
    re.setToZero();
    im.setToZero();
}

void initializeCustInt(void) {
    // for future use
}

int32_t custCplxHorner_TDCX(
    CustComplex& erg,
    const int32_t adeg,
    CustComplex* acoeff,
    CustComplex& x
) {
    int32_t error=0;
    erg.setToZero();

    for(int32_t k=adeg;k>=0;k--) {
        // erg = erg*x + coeff[k]
        CustComplex t1;
        error += custCplxMul_ZAB(t1,erg,x);
        error += custCplxAdd_ZAB(erg,t1,acoeff[k]);
    } // k

    // erg already set

    return error;
}

int32_t custHorner_TDCX(
    CustInterval& erg,
    const int32_t adeg,
    CustInterval* acoeff,
    CustInterval& x
) {
    int32_t error=0;
    erg.setToZero();

    for(int32_t k=adeg;k>=0;k--) {
        // erg = erg*x + coeff[k]
        CustInterval t1;
        error += custMul_ZAB(t1,erg,x);
        error += custAdd_ZAB(erg,t1,acoeff[k]);
    } // k

    // erg already set

    return error;
}

int32_t CustComplex::setStr(const char* astr) {
    int32_t error=0;
    error += re.setStr(astr); // takes care of values not fitting double

    im.setToZero();

    return error;
}

int32_t CustComplex::circumference_T(CustInterval& erg) {
    // circumference of the complex number's rectangle
    erg.setToZero();

    int32_t error=0;

    CustInterval t1,t2,width,height;

    // width
    t1.left=t1.right=re.left ;
    t2.left=t2.right=re.right;
    error += custSub_ZAB(width,t2,t1);

    // height
    t1.left=t1.right=im.left ;
    t2.left=t2.right=im.right;
    error += custSub_ZAB(height,t2,t1);

    // erg=2*height+2*width
    error += custAdd_ZAB(t1,height,width);
    t2.left=t2.right=2.0;
    error += custMul_ZAB(erg,t1,t2);

    return error;

}

// polynomial
custcplxPolynom::custcplxPolynom() {
    GRAD=-1; // not zero
}

void custcplxPolynom::clearCoeff(void) {
    GRAD=-1;

    for(int32_t k=0;k<MAXCUSTCPLXPOLYNOMGRAD;k++) {
        coeff[k].setToZero();
    } // k

}

int32_t custcplxPolynom::eval_TA(
    CustComplex& erg,
    CustComplex& z
) {
    int32_t error=0;

    // Horner-schme
    erg.setToZero();
    for(int32_t k=GRAD;k>=0;k--) {
        // erg <- erg*z+coeff[k]
        CustComplex t1;
        error += custCplxMul_ZAB(t1,erg,z);
        // erg <- t1+coeff[k]
        error += custCplxAdd_ZAB(erg,t1,coeff[k]);
        if (error != 0) return -1; // error

    } // k

    return error;

}

void custcplxPolynom::setCoeff_NA(
    const int32_t idx,
    CustComplex& A
) {
    if (
        (idx < 0) ||
        (idx >= MAXCUSTCPLXPOLYNOMGRAD)
    ) {
        LOGMSG2("\nerror. coefficient for invalid index %i to set\n",
            idx);
        exit(99);
    }

    coeff[idx].copyFrom( A );
    if (idx > GRAD) GRAD=idx;

}

void custcplxPolynom::ausgabe(FILE* f) {
    for(int32_t k=GRAD;k>=0;k--) {
        if (custcplxpredZero( coeff[k] ) > 0) continue;

        fprintf(f,"+z^%i*(",k);
        coeff[k].ausgabe(f);
        fprintf(f,")");
    } // k

}

int8_t custcplxpredZero(CustComplex& A) {
    // ret: <= 0 => not zero
    // > 1 => zero
    if (
        (fpclassify(A.re.left ) != FP_ZERO) ||
        (fpclassify(A.re.right) != FP_ZERO) ||
        (fpclassify(A.im.left ) != FP_ZERO) ||
        (fpclassify(A.im.right) != FP_ZERO)
    ) return 0;

    return 1;
}


#endif

