/*

	simple interval struct with end-points double
	based on constant upward rounding, no
	rounding mode switch necessary during operations

	Marc Meidlinger
	May 2020

	NO WARRANTY OF CORRECT RESULTS (although I did my best to
	ensure correctness of the implementation)

	CAUTION: Every arithmetics routine cust^NNN_ZAB()
			 ASSUMES rounding mode is set to UPWARDS

			 It is at the discretion of the user to
			 ensure this

	based on the article:

	Interval arithmetic with fixed rounding modw
	by S. Rump et al., 2016

*/

#ifndef _CUSTOMIZED_INTERVALH
#define _CUSTOMIZED_INTERVALH

#include "fenv.h"
#include "stdio.h"
#include "stdint.h"
#include "float.h"


// consts
const int32_t MAXCUSTCPLXPOLYNOMGRAD=64;

// get the current rounding mode, store it and
// set to upward infinity
struct CustRoundingUpwards {
	int32_t mode;

	CustRoundingUpwards();
	virtual ~CustRoundingUpwards();
};

struct CustRoundingDownwards {
	int32_t mode;

	CustRoundingDownwards();
	virtual ~CustRoundingDownwards();
};

// interval arithmetics with double and constant
// rounding upwards (downwards simulated by sign change)

struct CustInterval {
	double left,right;

	CustInterval();
	CustInterval(const double);
	CustInterval(const double,const double);
	void ausgabe(FILE*);
	void copyFrom( CustInterval& );
	int32_t setStr(const char*);
	void setToZero(void);
};

struct CustComplex {
	CustInterval re,im;

	CustComplex(const CustInterval&,const CustInterval&);
	CustComplex(const double,const double);
	CustComplex();
	void setToZero(void);
	void ausgabe(FILE*);
	void copyFrom( CustComplex& );
	int32_t setStr(const char*);
	int32_t circumference_T(CustInterval&);
};

struct CustIntervalMatrix2x2 {
	CustInterval x1y1,x2y1,x1y2,x2y2;

	CustIntervalMatrix2x2();
	CustIntervalMatrix2x2(const double,const double,const double,const double);

	void ausgabe(FILE*);
};

struct custcplxPolynom {
    int32_t GRAD;
    CustComplex coeff[MAXCUSTCPLXPOLYNOMGRAD];

    custcplxPolynom();

    void clearCoeff(void);
    int32_t eval_TA(CustComplex&,CustComplex&);
    void setCoeff_NA(const int32_t,CustComplex&);
    void ausgabe(FILE*);

};


// forward

// general functions
void terminate_program(void);
void initializeCustInt(void);
inline double minimumDouble(const double,const double,const double,const double);
inline double maximumDouble(const double,const double,const double,const double);
int8_t custcplxpredZero(CustComplex&);

// CustInterval operations
int32_t custAdd_ZAB(CustInterval&,CustInterval&,CustInterval&);
int32_t custSub_ZAB(CustInterval&,CustInterval&,CustInterval&);
int32_t custMul_ZAB(CustInterval&,CustInterval&,CustInterval&);
int32_t custDiv_ZAB(CustInterval&,CustInterval&,CustInterval&);
int32_t custSqrt_ZA_switch(CustInterval&,CustInterval&);
int32_t custSqrt_ZA(CustInterval&,CustInterval&);
int32_t custPow2_ZA(CustInterval&,CustInterval&);
int32_t custPow_ZAE(CustInterval&,CustInterval&,const int64_t);
void custPrint(FILE*,CustInterval&);
int32_t custMidpoint_ZA(CustInterval&,CustInterval&);
int32_t custSin_ZA(CustInterval&,CustInterval&);
int32_t custCos_ZA(CustInterval&,CustInterval&);


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


// assorted functions
int32_t custCplxHorner_TDCX(CustComplex&,const int32_t,CustComplex*,CustComplex&);
int32_t custHorner_TDCX(CustInterval&,const int32_t,CustInterval*,CustInterval&);

#endif
