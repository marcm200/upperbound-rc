/*

    affine arithmetics

    return value int32_t: < 0 => error
        == 0 => all valid

    NEEDS custInt and

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    needs constant UPWARD ROUNDING due to use
    of custInt
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    based on the articles:

    J. Stolfi, LH de Figueiredo. Self-validated numerical methods and applications. 1997.
    F. Messine, A. Touhami. A general reliable quadratic form: an extension of Affine Arithmetics. Reliable Computing, 2006.

*/

#ifndef _AFFINECPP
#define _AFFINECPP

#include "custint.h"


// conbst
const int32_t MAXSYMBOLSPERAA=256;
const int32_t MAXAFFINECPLXPOLYNOMGRAD=64;


// struct
struct AffinePath {
    int32_t nextfreesymbol;

    AffinePath();
    void init(void);

    int32_t getFreeSymbol_T(int32_t&);

};

struct AffineForm {
    AffinePath* belongs;
    double xi[MAXSYMBOLSPERAA]; // [0]=central value
    int32_t highestusedsymbol;
    int8_t used[MAXSYMBOLSPERAA]; // used[0] has no meaning
    // as xi[0] is always used as the central value of the affine form

    AffineForm();
    AffineForm(AffinePath*);
    void setBelongs(AffinePath*);
    void ausgabe(FILE*);
    void setToZero(void);
    void copyFrom( AffineForm& );
    int32_t getAbsWidthNonCentral_T(CustInterval&);

};

struct AffineComplex {
    AffineForm real,imag;

    AffineComplex();

    void ausgabe(FILE*);
    int32_t set_RIP(CustInterval&,CustInterval&,AffinePath&);
    int32_t set_by_copy_RI(AffineForm&,AffineForm&);
    AffinePath* getCommonBelongs(void);
    void setToZero(void);
    void copyFrom( AffineComplex& );

};

struct afcplxPolynom {
    int32_t GRAD;
    AffinePath* commonbelongs;
    AffineComplex coeff[MAXAFFINECPLXPOLYNOMGRAD];
    int8_t used[MAXAFFINECPLXPOLYNOMGRAD];

    afcplxPolynom();

    void clearCoeff(void);
    int32_t eval_TA(AffineComplex&,AffineComplex&);
    void setCoeff_NA(const int32_t,AffineComplex&);
    void ausgabe(FILE*);
    int32_t setFromCustPolynom_AP(custcplxPolynom&,AffinePath&);

};


// forward
int32_t convertToCustInt_TA(CustInterval&,AffineForm&);
int32_t setFromCustInt_TAP(AffineForm&,CustInterval&,AffinePath&);
int32_t spanCust_TRA(double&,double&,CustInterval&);


// defines

#define CHECKBELONG2(AA,BB) \
{\
    if (\
        ( (BB).belongs != (AA).belongs) ||\
        ( !((BB).belongs) ) ||\
        ( !((AA).belongs) )\
    ){\
        printf("\nerror. adding affine forms that do not belong to the same path is not valid.\n");\
        exit(99);\
    }\
}

#define CHECKBELONGCPLX2(AA,BB) \
{\
    AffinePath *ap=(AA).getCommonBelongs();\
    AffinePath *bp=(BB).getCommonBelongs();\
    \
    if (\
        (bp != ap) ||\
        (!ap) ||\
        (!bp)\
    ) {\
        LOGMSG("\nerror. belongs/complex not valid (non-identical or null).\n");\
        exit(99);\
    }\
}


// operations
int32_t afAdd_TAB_const(AffineForm&,AffineForm&,AffineForm&);
int32_t afSub_TAB_const(AffineForm&,AffineForm&,AffineForm&);
int32_t afMul_TAB_const(AffineForm&,AffineForm&,AffineForm&);
int32_t afNeg_TA_const(AffineForm&,AffineForm&);

// complex
int32_t afcplxAdd_TAB_const(AffineComplex&,AffineComplex&,AffineComplex&);
int32_t afcplxSub_TAB_const(AffineComplex&,AffineComplex&,AffineComplex&);
int32_t afcplxMul_TAB_const(AffineComplex&,AffineComplex&,AffineComplex&);


// routines

AffinePath::AffinePath() {
    nextfreesymbol=1; // stzarts with 1, noit zero
}

void AffinePath::init(void) {
    nextfreesymbol=1; // stzarts with 1, noit zero
}

int32_t convertToCustInt_TA(
    CustInterval& erg,
    AffineForm& A
) {
    int32_t error=0;

    // erg = central vlue + sum{i=1..highest] [-1..1]*xi if used[i]!= 0
    erg.left=erg.right=A.xi[0]; // central value
    CustInterval unit;
    unit.left =-1.0;
    unit.right= 1.0;

    //printf("\nhui: %i\n",A.highestusedsymbol);

    for(int32_t k=1;k<=A.highestusedsymbol;k++) {
        if (A.used[k] == 0) continue;

        /*
        printf("\ncurrent erg "); erg.ausgabe(stdout);
        printf("\n  plus %.20lg * ", A.xi[k]);
        printf(" * "); unit.ausgabe(stdout);
        */

        CustInterval ergcopy,t1,t2;
        // erg=erg+[xk..xk]*[-1..1]
        ergcopy.copyFrom(erg);
        t1.left=t1.right=A.xi[k];

        // erg=ergcopy+t1*unit
        error += custMul_ZAB(t2,t1,unit);

        //printf("\n  =(one) "); t2.ausgabe(stdout);

        // erg=ergcopy+t2
        error += custAdd_ZAB(erg,ergcopy,t2);

        //printf("\n  (err%i) after adding ",error); erg.ausgabe(stdout);

        if (error != 0) return -1; // error

    } // k

    return error;

}

int32_t setFromCustInt_TAP(
    AffineForm& erg,
    CustInterval& A,
    AffinePath& path
) {
    // A=[a..b] (can straddle 0)
    // compute midpoint m (or some other point inside)
    // compute w s.t. m-w <= a <= b <= m+w
    int32_t error=0;

    double m,R;
    error += spanCust_TRA(m,R,A);
    if (error != 0) return -1;

    erg.setToZero();

    // now it holds: m-R..m+R contains A
    erg.xi[0]=m; // 0 is a fixed index

    if (fpclassify(R) == FP_ZERO) {
        erg.highestusedsymbol=0; // so xi[0] is valud, but
        erg.used[0]=1; // but index 0 is always used
        erg.setBelongs(&path);
    } else {
        // get new symbol
        int32_t newsym;
        if (path.getFreeSymbol_T(newsym) != 0) return -1; // error

        if (newsym <= 0) {
            printf("\nimplementation error. symbol %i <= 0\n",newsym);
            exit(99);
        }

        // now all 1..newsym-1 as UNUSED
        for(int32_t k=1;k<newsym;k++) erg.used[k]=0;
        erg.used[newsym]=1;
        erg.xi[newsym]=R;
        erg.highestusedsymbol=newsym;
        erg.setBelongs(&path);
    }

    return error;

}

int32_t AffinePath::getFreeSymbol_T(int32_t& erg) {
    if (nextfreesymbol >= MAXSYMBOLSPERAA) {
        return -1; // error
    }

    erg=nextfreesymbol;
    nextfreesymbol++;
    //printf("\n\nnext free %i\n",nextfreesymbol);

    return 0;

}

// affineForm
void AffineForm::setBelongs(AffinePath* ap) {
    belongs=ap;
}

AffineForm::AffineForm() {
    highestusedsymbol=-1; // not zero
    belongs=NULL;

}

AffineForm::AffineForm(AffinePath* apath) {
    highestusedsymbol=-1; // not zero
    belongs=apath;
}

void AffineForm::ausgabe(FILE* f) {
    if (
        (!belongs) ||
        (highestusedsymbol < 0)
    ) {
        fprintf(f,"(undefined)");
        return;
    }

    fprintf(f,"(");

    for(int32_t k=0;k<=highestusedsymbol;k++) {
        if (
            (k > 0) &&
            (used[k] <= 0)
        ) continue;

        fprintf(f,"%+.20lg ",xi[k]);
        if (k > 0) fprintf(f,"*e%i",k);
    } // k
    fprintf(f,")");

}

void AffineForm::setToZero(void) {
    highestusedsymbol=0;
    xi[0]=0.0;

}

int32_t afSub_TAB_const(
    AffineForm& erg,
    AffineForm& A,
    AffineForm& B
) {
    int32_t error=0;

    // TODO: faster routine by directly using custSub
    // here, intuitive implementatioin:
    // adding the nagtuion of B

    AffineForm Bneg( B.belongs );
    error += afNeg_TA_const(Bneg,B);
    error += afAdd_TAB_const(erg,A,Bneg);

    return error;

}

int32_t afAdd_TAB_const(
    AffineForm& erg,
    AffineForm& A,
    AffineForm& B
) {
    // adding two affine forms
    // all variables must belong to the same path
    // A,B are NOT changed during the course
    // so A and B could be the same physical object
    // belongs erg will be set to the one of A,B

    CHECKBELONG2(A,B)

    erg.setToZero();
    erg.setBelongs(A.belongs);

    int32_t error=0;
    CustInterval rderroraccum;
    rderroraccum.setToZero();

    #define SFD_TA(TT,AA) \
    {\
        (TT).left=(TT).right=(AA);\
    }

    // add both x0 as an interval
    // take its midpoint (or something else inside)
    // and a covering radius R
    // the midpoint is the result's x0
    // the radius is a new noise symbol
    CustInterval ax0,bx0,t1;
    SFD_TA(ax0,A.xi[0]);
    SFD_TA(bx0,B.xi[0]);
    if (custAdd_ZAB(t1,ax0,bx0) != 0) return -1; // error

    #define ADDRDERROR(WW) \
    {\
        CustInterval tmp1,rdcopy;\
        tmp1.left=tmp1.right=(WW);\
        rdcopy.copyFrom(rderroraccum);\
        if (custAdd_ZAB(rderroraccum,rdcopy,tmp1) != 0) return -1;\
    }

    double x0m,x0R;
    if (spanCust_TRA(x0m,x0R,t1) != 0) return -1;
    erg.xi[0]=x0m;
    ADDRDERROR(x0R)

    // for the parameterx xi i=1..
    // add Axi+Bxi as interval
    // again ake the midpoints and a coverign radius
    // and add that radius the the above symbol (upward rpounding)

    // erg is a valid AffineFOrm,e verything set
    // excapet the rounding accumated error
    // must be assigned to a new noise symbol
    int32_t maxidx=A.highestusedsymbol;
    if (B.highestusedsymbol > maxidx) maxidx=B.highestusedsymbol;

    for(int32_t k=1;k<=maxidx;k++) {
        erg.used[k]=0;
        CustInterval ta,tb;

        if (
            (k > A.highestusedsymbol) ||
            (A.used[k] <= 0)
        ) {
            ta.setToZero();
        } else {
            SFD_TA(ta,A.xi[k]);
        }

        if (
            (k > B.highestusedsymbol) ||
            (B.used[k] <= 0)
        ) {
            tb.setToZero();
        } else {
            SFD_TA(tb,B.xi[k]);
        }

        CustInterval add;
        if (custAdd_ZAB(add,ta,tb) != 0) return -1;

        // midpoint,R
        double m,R;
        if (spanCust_TRA(m,R,add) != 0) return -1;

        if (
            (fpclassify(m) == FP_ZERO) &&
            (fpclassify(R) == FP_ZERO)
        ) {
            // nothing to add
            erg.used[k]=0;
            continue;
        }

        erg.used[k]=1;
        erg.xi[k]=m;

        ADDRDERROR(R)

    } // k

    erg.highestusedsymbol=maxidx;

    // now is there a rounding error ?
    if (
        (fpclassify(rderroraccum.left ) != FP_ZERO) ||
        (fpclassify(rderroraccum.right) != FP_ZERO)
    ) {

        // now new noise symbol for accumalted rounding errors
        int32_t newsym;
        // all belogns are identical
        if (erg.belongs->getFreeSymbol_T(newsym) != 0) return -1; // error

        if (
            (newsym <= 0) ||
            (newsym <= erg.highestusedsymbol)
        ) {
            printf("\nimplementation error. symbol %i <= 0 or already used indices\n",newsym);
            exit(99);
        }

        // now all erg.hgiehst+1..newsym-1 as UNUSED
        for(int32_t k=(erg.highestusedsymbol+1);k<newsym;k++) erg.used[k]=0;
        erg.used[newsym]=1;

        // now take max{|rd.left|,|rd.right|) (could be rounding otuwards into
        // the negative)
        if (rderroraccum.right < 0.0) {
            printf("\naccumulated rounding error negative.\n");
            exit(99);
        }

        double news=rderroraccum.right;
        double might=rderroraccum.left;
        if (might < 0.0) might=-might; // exact operation
        if (might > news) news=might;

        erg.xi[newsym]=news;
        erg.highestusedsymbol=newsym;

    } // rounding error present

    return error;

}

int32_t spanCust_TRA(
    double& m,
    double& R,
    CustInterval& A
) {
    int32_t error=0;

    CustInterval minterval;
    error += custMidpoint_ZA(minterval,A);
    if (error != 0) return -1;

    //printf("\nset int "); A.ausgabe(stdout);
    //printf("\nm= "); minterval.ausgabe(stdout);

    int8_t setm=0;

    #define TESTM(WW) \
    {\
        if (setm <= 0) {\
            if (\
                (A.left <= (WW) ) &&\
                ( (WW) <= A.right )\
            ) {\
                m=WW;\
                setm=1;\
            }\
        }\
    }

    // is m inside the interval (could've been rounded outwards)
    TESTM(minterval.left);
    TESTM(minterval.right);

    if (setm <= 0) {
        // cannot find midpoint inside A
        return -1; // error
    }

    //printf("\nm to %.20lg",m);

    CustInterval ma,mb,ivm,iva,ivb;
    ivm.left=ivm.right=m;
    iva.left=iva.right=A.left;
    ivb.left=ivb.right=A.right;

    error += custSub_ZAB(ma,ivm,iva);
    error += custSub_ZAB(mb,ivb,ivm);
    if (error != 0) return -1;

    //printf("\n|m-a| >= "); ma.ausgabe(stdout);
    //printf("\n|m-b| >= "); mb.ausgabe(stdout);

    if (
        (ma.right < 0.0) ||
        (mb.right < 0.0)
    ) {
        printf("\nerror 1. set\n");
        exit(99);
    }

    // take larger RIGHT end
    R=ma.right;
    if (mb.right > R) R=mb.right;

    // now it holds: m-R..m+R contains A
    // check this
    CustInterval m2,mm,rr;
    rr.left =-R; // exact operation
    rr.right= R;
    mm.left=mm.right=m;
    error += custAdd_ZAB(m2,mm,rr);
    if (error != 0) return -1;

    if (
        !(
            (m2.left <= A.left) &&
            (A.right <= m2.right)
        )
    ) {
        printf("\nerror implementation. not coveriung.\n");
        exit(99);
    }

    // it holds: m-R <= A.left <= A.right <= m+R in exact (and
    // then also in outward rouding) operations

    return error;

}

int32_t afNeg_TA_const(
    AffineForm& erg,
    AffineForm& A
) {
    // negatives an affine form
    // this is an EXACT operation as only signs are switched

    erg.copyFrom( A );
    // start from index 0
    for(int32_t k=0;k<=erg.highestusedsymbol;k++) {
        erg.xi[k] = -erg.xi[k]; // exact operation
    } // k

    return 0; // no error can occur

}

void AffineForm::copyFrom( AffineForm& A ) {
    belongs = A.belongs; // copy blong ptr
    highestusedsymbol = A.highestusedsymbol;

    for(int32_t k=0;k<=A.highestusedsymbol;k++) {
        used[k] = A.used[k];
        xi[k] = A.xi[k];
    } // k

    // no error possible

}

int32_t afMul_TAB_const(
    AffineForm& erg,
    AffineForm& A,
    AffineForm& B
) {
    int32_t error=0;

    /*

        erg = x_0y_0 + sum{ e_i*(x_0y_i + x_iy_0) }
                + e_new*( (sum |xi|) * (sum |yi|) )

        adding all rounding errors into e_new

        from "A general relialbe quadratic form"
        by Fredrik Messine
        p174

    */

    CHECKBELONG2(A,B)

    erg.setToZero();
    erg.setBelongs(A.belongs);
    // erg0=x0y0 => rounding in new
    // erg1..maxidx: (x0yi+xiy0) => rounding in new
    // (sum|xi|) * (sum|<i|) => rounding in new
    #define SFD_TA(TT,AA) \
    {\
        (TT).left=(TT).right=(AA);\
    }

    #define ADDRDERROR(WW) \
    {\
        CustInterval tmp1,rdcopy;\
        tmp1.left=tmp1.right=(WW);\
        rdcopy.copyFrom(rderroraccum);\
        if (custAdd_ZAB(rderroraccum,rdcopy,tmp1) != 0) return -1;\
    }

    CustInterval rderroraccum;
    rderroraccum.setToZero();

    // erg0: x0y0
    erg.highestusedsymbol=0;
    CustInterval x0,y0,x0y0;
    SFD_TA(x0,A.xi[0])
    SFD_TA(y0,B.xi[0])
    if (custMul_ZAB(x0y0,x0,y0) != 0) return -1;
    double m,R;
    if (spanCust_TRA(m,R,x0y0) != 0) return -1; // error
    erg.xi[0]=m;
    ADDRDERROR(R)

    //printf("\nx0y0: %.20lg + %.20lg\n",m,R);

    // erg1..n: (x9yi+xiy0)
    int32_t maxidx=A.highestusedsymbol;
    if (B.highestusedsymbol > maxidx) maxidx=B.highestusedsymbol;

    CustInterval sumabsxi,sumabsyi;
    sumabsxi.setToZero();
    sumabsyi.setToZero();

    for(int32_t k=1;k<=maxidx;k++) {
        erg.used[k]=0;
        // x0yk+xky0
        CustInterval xk,yk,t1,t2,t3;

        if (
            (k > A.highestusedsymbol) ||
            (A.used[k] <= 0)
        ) {
            SFD_TA(xk,0)
        } else {
            SFD_TA(xk,A.xi[k])
        }
        if (xk.left != xk.right) {
            printf("\nimplementation error. interval not a point\n");
            exit(99);
        }

        CustInterval sumcopy;
        sumcopy.copyFrom( sumabsxi );
        CustInterval absxk;
        double absx;
        if (xk.left < 0.0) absx=-xk.left; else absx=xk.left;
        SFD_TA(absxk,absx)
        error += custAdd_ZAB(sumabsxi,sumcopy,absxk);

        if (
            (k > B.highestusedsymbol) ||
            (B.used[k] <= 0)
        ) {
            SFD_TA(yk,0)
        } else {
            SFD_TA(yk,B.xi[k])
        }

        if (yk.left != yk.right) {
            printf("\nimplementation error. interval not a point\n");
            exit(99);
        }

        sumcopy.copyFrom( sumabsyi );
        CustInterval absyk;
        double absy;
        if (yk.left < 0.0) absy=-yk.left; else absy=yk.left;
        SFD_TA(absyk,absy)
        error += custAdd_ZAB(sumabsyi,sumcopy,absyk);

        error += custMul_ZAB(t1,x0,yk);
        // t1+xky0
        error += custMul_ZAB(t2,xk,y0);
        // t1+t2
        error += custAdd_ZAB(t3,t1,t2);
        if (error != 0) return -1;
        // t3
        double m,R;
        if (spanCust_TRA(m,R,t3) != 0) return -1;

        if (
            (fpclassify(m) == FP_ZERO) &&
            (fpclassify(R) == FP_ZERO)
        ) {
            erg.used[k]=0;
            continue;
        }

        erg.xi[k]=m;
        ADDRDERROR(R);
        erg.used[k]=1;
    } // k
    erg.highestusedsymbol=maxidx;

    // affine approx quadratic term
    // new += (sum|xi|) * (sum|yi|)
    CustInterval summul;
    error += custMul_ZAB(summul,sumabsxi,sumabsyi);
    if (error != 0) return -1;
    if (spanCust_TRA(m,R,summul) != 0) return -1;

    if (m < 0.0) return -1; // as |xi|*yi| >= 0, if here < 0 => rounding

    // m is quasi a rounding error
    ADDRDERROR(m)

    // and the width of the interval
    ADDRDERROR(R)

    // now is there a rounding error ?
    if (
        (fpclassify(rderroraccum.left ) != FP_ZERO) ||
        (fpclassify(rderroraccum.right) != FP_ZERO)
    ) {

        // now new noise symbol for accumalted rounding errors
        int32_t newsym;
        // all belogns are identical
        if (erg.belongs->getFreeSymbol_T(newsym) != 0) return -1; // error

        if (
            (newsym <= 0) ||
            (newsym <= erg.highestusedsymbol)
        ) {
            printf("\nimplementation error. symbol %i <= 0 or already used indices\n",newsym);
            exit(99);
        }

        // now all erg.hgiehst+1..newsym-1 as UNUSED
        for(int32_t k=(erg.highestusedsymbol+1);k<newsym;k++) erg.used[k]=0;
        erg.used[newsym]=1;

        // now take max{|rd.left|,|rd.right|) (could be rounding otuwards into
        // the negative)
        if (rderroraccum.right < 0.0) {
            printf("\naccumulated rounding error negative.\n");
            exit(99);
        }

        double news=rderroraccum.right;
        double might=rderroraccum.left;
        if (might < 0.0) might=-might; // exact operation
        if (might > news) news=might;

        erg.xi[newsym]=news;
        erg.highestusedsymbol=newsym;

    } // rounding error present

    return error;

}

void AffineComplex::copyFrom( AffineComplex& A ) {
    real.copyFrom( A.real );
    imag.copyFrom( A.imag );
}

void AffineComplex::setToZero(void) {
    real.setToZero();
    imag.setToZero();
}

AffineComplex::AffineComplex() {
}

void AffineComplex::ausgabe(FILE* f) {
    fprintf(flog,"[");
    real.ausgabe(f);
    fprintf(f,"]+i*[");
    imag.ausgabe(f);
    fprintf(f,"]");
}

int32_t AffineComplex::set_RIP(
    CustInterval& ar,
    CustInterval& ai,
    AffinePath& apath
) {
    int32_t error=0;

    error += setFromCustInt_TAP(real,ar,apath);
    error += setFromCustInt_TAP(imag,ai,apath);

    return error;
}

int32_t AffineComplex::set_by_copy_RI(
    AffineForm& ar,
    AffineForm& ai
) {
    if (
        (!ar.belongs) ||
        (!ai.belongs) ||
        (ar.belongs != ai.belongs)
    ) return -1; // not set

    real.copyFrom( ar );
    imag.copyFrom( ai );

    return 0;

}

int32_t afcplxAdd_TAB_const(
    AffineComplex& erg,
    AffineComplex& A,
    AffineComplex& B
) {
    CHECKBELONGCPLX2(A,B)

    // erg.real=A.real+B.real
    // erg.imag=A.imag+B.imag

    int32_t error=0;

    error += afAdd_TAB_const(erg.real,A.real,B.real);
    error += afAdd_TAB_const(erg.imag,A.imag,B.imag);

    return error;

}

int32_t afcplxSub_TAB_const(
    AffineComplex& erg,
    AffineComplex& A,
    AffineComplex& B
) {
    CHECKBELONGCPLX2(A,B)

    int32_t error=0;

    error += afSub_TAB_const(erg.real,A.real,B.real);
    error += afSub_TAB_const(erg.imag,A.imag,B.imag);

    return error;

}

int32_t afcplxMul_TAB_const(
    AffineComplex& erg,
    AffineComplex& A,
    AffineComplex& B
) {
    CHECKBELONGCPLX2(A,B)

    int32_t error=0;
    // (ar+ai*i)*(br+bi*i)
    // = ar*br + ar*bi*i + ai*br*i + ai*bi*i^2
    // = ar*br-ai*bi + i*(ar*bi+ai*br)

    // erg.real=A.real*B.real - A.imag*B.imag
    AffineForm t1(A.real.belongs),t2(A.real.belongs);
    error += afMul_TAB_const(t1,A.real,B.real);
    // erg.real=t1 - A.imag*B.imag
    error += afMul_TAB_const(t2,A.imag,B.imag);
    // erg.real=t1 - t2
    error += afSub_TAB_const(erg.real,t1,t2);

    // erg.imag = A.real*B.imag + A:imag*B.real
    error += afMul_TAB_const(t1,A.real,B.imag);
    // erg.imag = t1 + A:imag*B.real
    error += afMul_TAB_const(t2,A.imag,B.real);
    // erg.imag = t1 + t2
    error += afAdd_TAB_const(erg.imag,t1,t2);


    return error;

}

AffinePath* AffineComplex::getCommonBelongs(void) {
    if (real.belongs != imag.belongs) return NULL;

    return real.belongs; // might still be zero
}

// polynom
afcplxPolynom::afcplxPolynom() {
    GRAD=-1; // not zero
    commonbelongs=NULL;
}

void afcplxPolynom::clearCoeff(void) {
    GRAD=-1;
    commonbelongs=NULL;
    for(int32_t k=0;k<MAXAFFINECPLXPOLYNOMGRAD;k++) {
        used[k]=0;
    }

}

int32_t afcplxPolynom::eval_TA(
    AffineComplex& erg,
    AffineComplex& z
) {
    if (
        (z.real.belongs != commonbelongs) ||
        (z.imag.belongs != commonbelongs) ||
        (!z.real.belongs) ||
        (!z.imag.belongs) ||
        (!commonbelongs)
    ) {
        LOGMSG("\nerror. cannot evaluate a polynomial with differing AffinePath\n");
        exit(99);
    }

    int32_t error=0;

    erg.setToZero();
    erg.real.setBelongs( z.real.belongs );
    erg.imag.setBelongs( z.real.belongs );
    // erg must be a valid affineForm tzo be used in
    // calculations

    // Horner-schme
    for(int32_t k=GRAD;k>=0;k--) {
        // erg <- erg*z+coeff[k]
        AffineComplex t1;
        error += afcplxMul_TAB_const(t1,erg,z);
        // erg <- t1+coeff[k]
        if (used[k] > 0) {
            error += afcplxAdd_TAB_const(erg,t1,coeff[k]);
        } else {
            erg.copyFrom( t1 );
        }

        if (error != 0) return -1; // error

    } // k

    return error;

}

void afcplxPolynom::setCoeff_NA(
    const int32_t idx,
    AffineComplex& A
) {
    if (
        (idx < 0) ||
        (idx >= MAXAFFINECPLXPOLYNOMGRAD)
    ) {
        LOGMSG2("\nerror. setting coefficient out of index rangfe %i\n",idx);
        exit(99);
    }

    if (A.real.belongs != A.imag.belongs) {
        LOGMSG("\nerror. coefficient has differing AffinePath\n");
        exit(99);
    }

    if (!commonbelongs) commonbelongs=A.real.belongs;

    if (
        (commonbelongs != A.real.belongs) ||
        (commonbelongs != A.imag.belongs)
    ) {
        LOGMSG("\nerror. all coefficients must belong to the same AffinePath\n");
        exit(99);
    }

    used[idx]=1;
    coeff[idx].copyFrom( A );

    if (idx > GRAD) GRAD=idx;

}

void afcplxPolynom::ausgabe(FILE* f) {

    for(int32_t k=GRAD;k>=0;k--) {
        if (used[k] <= 0) continue;

        if (k == 0) fprintf(f,"+(");
        else if (k == 1) fprintf(f,"+z*(");
        else fprintf(f,"+z^%i*(",k);

        coeff[k].ausgabe(f);
        fprintf(f,")");
    } // k

}

int32_t afcplxPolynom::setFromCustPolynom_AP(
    custcplxPolynom& cpoly,
    AffinePath& apath
) {
    int32_t error=0;

    clearCoeff();

    //printf("\n\nafpoly: cpoly: "); cpoly.ausgabe(stdout);

    for(int32_t k=0;k<=cpoly.GRAD;k++) {

        if (custcplxpredZero( cpoly.coeff[k] ) > 0) continue;

        //printf("\nk:%i: cpoly: ",k);
        //cpoly.coeff[k].ausgabe(stdout);

        AffineComplex t1;
        error += setFromCustInt_TAP(t1.real,cpoly.coeff[k].re,apath);
        error += setFromCustInt_TAP(t1.imag,cpoly.coeff[k].im,apath);

        //printf("\n  afset(err%i): ",error);
        //t1.ausgabe(stdout);

        setCoeff_NA(k,t1);

        //printf("\ntgtgrad: %i",GRAD);

    } // k

    return error;

}

int32_t AffineForm::getAbsWidthNonCentral_T(CustInterval& erg) {

    int32_t error=0;

    // sum |xi| EXCLUDING x0
    erg.setToZero();

    for(int32_t k=1;k<=highestusedsymbol;k++) {
        if (used[k] <= 0) continue;

        CustInterval t1,ergcopy;
        if (xi[k] < 0.0) t1.left=t1.right=-xi[k]; // exact operation
        else t1.left=t1.right=xi[k];

        ergcopy.copyFrom(erg);
        error += custAdd_ZAB(erg,ergcopy,t1);

    } // k

    return error;

}

#endif

