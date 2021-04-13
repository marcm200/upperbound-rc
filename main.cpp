/*

	Guaranteed Upper Area Bound for the Mandelbrot Set

	reliable orbit constructiuon using
	inerval arithmetics (number type custComplex)
	and affine arithmetics (AffineComplex)

	Marc Meidlinger
	April 2021

*/

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "time.h"
#include "fenv.h"
#include "custint.h"
#include "stdint.h"


// consts

enum {
	FUNC_Z2C=0,

	FUNCANZ
};

// for use in areaCombine
enum {
    FORMAT_RAW=1,
    FORMAT_COMPRESSED
};

const int16_t MAXCOMPRESSEDPER=30000;
const int32_t COMBINEEOL=-1;
const int32_t COMBINEEOF=-2;
const int32_t INT32MAX=(int32_t)1 << 30;
const int64_t INT48=(int64_t)1 << 48;
const int64_t INT61=(int64_t)1 << 61;

const char funcname[][32] = {
	"Z2C"
};

const int8_t INP12_ERROR=-1;
const int8_t INP12_MAYBE=0;
const int8_t INP12_INSIDE=1;
const int8_t INP12_OUTSIDE=2;

int32_t TILELEVEL=15;
const int32_t MAXTILELEVEL=18;
const int32_t MAXREFINEMENTM=30;

// for Mset, at most 63
// not all are used in this project
const int32_t COLOR_UNDEF=-1;
const uint8_t MTILE_GRAY=0b00;
const uint8_t MTILE_WHITE=0b01; // for future compatibility
const uint8_t MTILE_BLACK=0b10;
const uint8_t MTILE_CHECKFORBLACK=4;
const uint8_t MTILE_CHECKFORWHITE=MTILE_CHECKFORBLACK;
const uint8_t MTILE_GRAY_DONE=6;
const uint8_t MTILE_POTENTIALBLACK=7;
const uint8_t MTILE_POTENTIALHC=8;
const uint8_t MTILE_TEMP_TOFOLLOW=9;
const uint8_t MTILE_TEMP_FOLLOWED=10;
const uint8_t MTILE_NUMERICALCYCLE=11;
const int32_t RGBANZ=12;

const uint8_t DICT_GRAY=MTILE_GRAY;
const uint8_t DICT_NUMERICALCYCLE=MTILE_NUMERICALCYCLE;
const uint8_t DICT_WHITE=MTILE_WHITE;
const uint8_t DICT_BLACK=MTILE_BLACK;

// color RGB values
// must be signed int32_t as correctColor subtracts
// a radius value
const int32_t rgbs[][3] = {
	/* 0 */ {127,127,127}, 	/* MTILE_GRAY */
	/* 1 */ {255,255,255}, 	/* white */
	/* 2 */ {0,0,0},		/* black */
	/* 3 */ {255,0,0}, 		/* not used */
	/* 4 */ {0,0,255}, 		/* check for black */
	/* 5 */ {255,0,0}, 		/* not used */
	/* 6 */ {255,255,0}, 	/* done */
	/* 7 */ {63,63,63},     /* potential black */
	/* 8 */ {255,255,0},    /* potential black and hc */
	/* 9 */ {63,63,63},     /* temporary color */
	/* 10 */ {255,127,127}, /* temporary color */
	/* 11 */ {255,0,0}     /* numerical cycle */
};


// structs

// corner coordinates of an Mtile or Jtile
// values are accurately represented as only dyadic fractions
// are assigned
struct Calcs {
    int64_t ctrw,ctrg,ctrb,ctrnum;

    void setToZero(void);
};

// pixel coordinates of a virtual screen
struct ScreenRect64 {
	int64_t x0,x1,y0,y1;

	void copyFrom(ScreenRect64&);
	void ausgabe(FILE*);
};

// chunk size dependend on operating system
//#define _CHUNK512

const int32_t MAXPTR=2048;

#ifdef _CHUNK512
const uint64_t CHARMAPCHUNKSIZE=( (uint64_t)1 << 27 );
#else
// win64 => 1 GB
const uint64_t CHARMAPCHUNKSIZE=( (uint64_t)1 << 30 );
#endif

typedef unsigned char BYTE;

struct RGB3 {
	BYTE R,G,B;
};

typedef BYTE *PBYTE;

// memory manager for byte allocations
struct ArrayByteMgr {
	BYTE* current;
	int32_t allocatedIdx,freeFromIdx,allocatePerBlockIdx;
	PBYTE ptr[MAXPTR];
	int32_t anzptr;

	ArrayByteMgr();
	void FreeAll(void);
	virtual ~ArrayByteMgr();
	PBYTE getMemory(const int32_t);
};

// stores a full Jset or the part in a virtual screen where the
// active region of an Mset hyperbolic component resides
// data structure is equivalent to an 256-color bitmap
// data can be viewed externally
// colors used are: black = interior, gray = untested
// yellow = tested and judged gray at current resolution
struct TiledPlane {
	int64_t xlen,ylen;
	int64_t memused;
	PBYTE* cmpY;
	RGB3 palette[256];

	TiledPlane();
	virtual ~TiledPlane();

	void setlenxy(const int64_t,const int64_t);
	void saveFormatBase(const int32_t,const char*);
	void saveBitmap(const char*);
	void saveTwice_FXY(const int32_t,const int64_t,const int64_t);
	int8_t saveCompressed(const char*);
	int8_t loadCompressed(const char*);
	int8_t loadFormatBase(const int32_t,const char*);
	void lineHorizVert(const int32_t,const int32_t,const int32_t,const int32_t,const BYTE);
	int8_t loadBitmap(const char*);
	void fillrect(const int32_t,const int32_t,const int32_t,const int32_t,const BYTE);
	void setPaletteRGB(const int32_t,const BYTE,const BYTE,const BYTE);
	void setPunkt(const int32_t,const int32_t,const BYTE);
	BYTE getPunkt(const int32_t,const int32_t);
	int8_t copyFrom( TiledPlane& );
};


// globals 1
TiledPlane dict;

// affine arithemtics
int64_t GLOBALHIGHESTUSEDSYMBOL=0;
int64_t AFFINEFOLDBACK=-1;

// counter
int64_t ctrexterior=0;
int64_t ctrfoldback=0;

// control flow
int8_t USEGRIDFILL=1;
int8_t USEAFFINE=1;
int8_t P123ON=0;
int8_t USEDICTIONARY=1;
int8_t SAVETWICE=1;
int8_t CREATESTART=0;

// heuristics
int32_t GRIDLOOPS=8;
int64_t GRIDFILLJUMP=9; // sth non-dyadic
int64_t SKIPIFLESS=-1; // off

// dictionary
int64_t DICTLEVEL=0;
int64_t DICTLEN=0;

// general
int64_t GRIDWIDTH=17;
int32_t MAXIT=-1;
int64_t TILELEN=(int64_t)1 << TILELEVEL;
int8_t SYMMETRY_TO_XAXIS=0;
int32_t LOADFORMAT=FORMAT_COMPRESSED;
int32_t SAVEFORMAT=FORMAT_COMPRESSED;
int8_t SYMMETRY_TO_YAXIS=1;
int32_t FUNC=FUNC_Z2C;
FILE *flog=NULL;
ArrayByteMgr planemgr;
TiledPlane md2; // for refinement
int64_t tilepixelsint=0;

int64_t time_at_enter_main=0; // in CLOCKS_PER_SEC
int64_t endtime=-1; // in CLOCKS_PER_SEC

// Mset
int32_t REFINEMENTM=10;
double scaleRangePerPixelM,scalePixelPerRangeM;
int32_t SCREENWIDTHM;
int32_t RANGEM0=-2,RANGEM1=2;
double COMPLETEM0,COMPLETEM1;


#define SETCPLX(TT,RR,II) \
{\
    (TT).re.left=(TT).re.right=RR;\
    (TT).im.left=(TT).im.right=II;\
}

#define LOGMSG4(TT,AA,BB,CC) \
{\
	fprintf(flog,TT,AA,BB,CC); fflush(flog);\
	printf(TT,AA,BB,CC);\
}

#define LOGMSG5(TT,AA,BB,CC,DD) \
{\
	fprintf(flog,TT,AA,BB,CC,DD); fflush(flog);\
	printf(TT,AA,BB,CC,DD);\
}

#define LOGMSG(TT) \
{\
	fprintf(flog,TT); fflush(flog);\
	printf(TT);\
}

#define LOGMSG2(TT,AA) \
{\
	fprintf(flog,TT,AA); fflush(flog);\
	printf(TT,AA);\
}

#define LOGMSG3(TT,AA,BB) \
{\
	fprintf(flog,TT,AA,BB); fflush(flog);\
	printf(TT,AA,BB);\
}

// include number type support for
// interval and affine arithmetics

#include "custint_local.cpp"
#include "affine.cpp"


// globals for number types used
AffinePath globalafpath;
CustInterval scaleRangePerPixelMcust;
CustInterval COMPLETEM0cust,COMPLETEM1cust;
CustComplex cplx4,cplx7,cplx8,cplxm4;

// declarations
void write2(FILE*,const uint8_t,const uint8_t);
void write4(FILE*,const uint8_t,const uint8_t,const uint8_t,const uint8_t);


// forward
void saveTwice_BFXY(TiledPlane&,const int32_t,const int64_t,const int64_t);
int32_t sum_int32t(const int32_t a,const int32_t b);
char* getFnBase_TSXY(char*,const int64_t,const int64_t,const int64_t);
int8_t predP1_z2c(CustComplex&);
int8_t predP2_z2c(CustComplex&);
void saveDictionaryAndTwice(void);


// routines
inline int8_t rectsOverlap(
    CustComplex& A,
	const double bx0,const double bx1,
	const double by0,const double by1
) {
	// it is not to assume any of the rectangles
	// is smaller than the other

	if (
		(A.re.right < bx0) || /* A fully left b */
		(A.re.left  > bx1) || /* A fully right b */
		(A.im.right < by0) || /* A below b */
		(A.im.left  > by1) /* A above b */
	) {
		return 0;
	}

	return 1;
}

char* upper(char* s) {
	if (!s) return NULL;

	for(int32_t i=0;i<(int32_t)strlen(s);i++) {
		if ((s[i]>='a')&&(s[i]<='z')) s[i]=s[i]-'a'+'A';
	}

	return s;
}

void setPaletteTo(TiledPlane& md) {
	// image colors
	for(int32_t i=0;i<256;i++) md.setPaletteRGB(i,255,0,0);
	for(int32_t i=0;i<RGBANZ;i++) {
        md.setPaletteRGB(i,
            rgbs[i][0],
            rgbs[i][1],
            rgbs[i][2]
        );
    }
}

// struct TiledPlane
void write2(FILE *f,const BYTE a,const BYTE b) {
	fwrite(&a,1,sizeof(a),f);
	fwrite(&b,1,sizeof(b),f);
}

void write4(FILE *f,const BYTE a,const BYTE b,const BYTE c,const BYTE d) {
	fwrite(&a,1,sizeof(a),f);
	fwrite(&b,1,sizeof(b),f);
	fwrite(&c,1,sizeof(c),f);
	fwrite(&d,1,sizeof(d),f);
}

void saveTwice_BFXY(
    TiledPlane& md,
    const int32_t sf,
    const int64_t filex,
    const int64_t filey
) {
    // creates a tempory image in md2
    // of each quadrant

    int64_t SC2=(SCREENWIDTHM << 1);
    int64_t HALFLEN=TILELEN >> 1;

    md2.setlenxy(TILELEN,TILELEN);
    setPaletteTo(md2);

    for(int64_t dy=0;dy<=1;dy++) {
        int64_t y0=dy*HALFLEN,y1=(dy+1)*HALFLEN;

        for(int64_t dx=0;dx<=1;dx++) {
            int64_t x0=dx*HALFLEN,x1=(dx+1)*HALFLEN;

            int64_t tgty=0;
            for(int64_t srcy=y0;srcy<y1;srcy++) {
                int64_t tgtx=0;

                for(int64_t srcx=x0;srcx<x1;srcx++) {
                    int32_t f=md.getPunkt(srcx,srcy);

                    // only BLACK,WHITE,GRAY in TWICE image
                    if (
                        (f == MTILE_WHITE) ||
                        (f == MTILE_BLACK) ||
                        (f == MTILE_NUMERICALCYCLE)
                    ) {
                        // keep color
                    } else {
                        // set everything else tzo gray
                        f=MTILE_GRAY;
                    }

                    md2.setPunkt(tgtx  ,tgty  ,f);
                    md2.setPunkt(tgtx+1,tgty  ,f);
                    md2.setPunkt(tgtx  ,tgty+1,f);
                    md2.setPunkt(tgtx+1,tgty+1,f);

                    tgtx += 2;

                } // srcx

                tgty += 2;

            } // srcy

            char fn[2048];
            getFnBase_TSXY(fn,SC2,
                (filex << 1) + dx*TILELEN,
                (filey << 1) + dy*TILELEN
            );

            md2.saveFormatBase(SAVEFORMAT,fn);

        } // dx

    } // dy

}

int8_t TiledPlane::copyFrom( TiledPlane& src ) {
    setlenxy( src.xlen,src.ylen );
    if (
        (xlen != src.xlen) ||
        (ylen != src.ylen)
    ) return 0;

    for(int32_t i=0;i<=255;i++) {
        palette[i].R=src.palette[i].R;
        palette[i].G=src.palette[i].G;
        palette[i].B=src.palette[i].B;
    }

    for(int64_t y=0;y<src.ylen;y++) {
        for(int64_t x=0;x<src.xlen;x++) {
            setPunkt(x,y,src.getPunkt(x,y));
        } // x
    } // y

    return 1;

}

int8_t TiledPlane::loadCompressed(const char* afn) {

    #define READF(TYP,WW) \
    {\
        TYP tmpw;\
        if (fread(&tmpw,1,sizeof(tmpw),f) != sizeof(tmpw)) {\
            fclose(f);\
            return 0;\
        }\
        WW=tmpw;\
    }

    FILE *f=fopen(afn,"rb");
    if (!f) return 0;

    int x64,y64;
    READF(int64_t,x64)
    READF(int64_t,y64)

    setlenxy(x64,y64);
    setPaletteTo(*this);

    if ( (xlen != x64) || (ylen != y64) ) {
        fclose(f);
        return 0; // error, probably memory
    }

    int64_t y=0,x=0;

    while (1) {
        uint8_t color;
        int16_t anz;

        READF(uint8_t,color)
        READF(int16_t,anz)

        if (anz == COMBINEEOF) {
            // EOF
            break;
        }

        if (anz == COMBINEEOL) {
            // EOLine
            if (x != xlen) {
                // error
                printf("\nerror. reading line %I64d\n",y);
                fclose(f);
                exit(99);
            }
            x=0;
            y++;
            continue;
        }

        if (anz == 0) {
            // error
            printf("\nstreak 0 not allowed\n");
            exit(99);
        }

        int64_t xe=x+anz-1;
        lineHorizVert(x,y,xe,y,color);
        x=xe + 1;

    } // while

    if (y != ylen) { // or y64, the same
        // error
        fclose(f);
        return 0;
    }

    fclose(f);

    return 1;
}

int8_t TiledPlane::saveCompressed(const char* afn) {
    // stores the imnage data in a compressed form
    /*

        File-format:
        int64_t xlen
        int64_t ylen
        streaks: {color uint8_t,anz int16_t}
        End-of-line: {0,-1}
        End-of-File: {0,-2}

    */

    #define WRITEF(TYP,WW) \
    {\
        TYP tmpw=WW;\
        if (fwrite(&tmpw,1,sizeof(tmpw),f) != sizeof(tmpw)) {\
            fclose(f);\
            return 0;\
        }\
    }

    #define SAVESTREAK_ANZF(ANZ,FF) \
    {\
        int64_t tmpa=ANZ;\
        while (tmpa >= MAXCOMPRESSEDPER) {\
            WRITEF(uint8_t,FF)\
            WRITEF(int16_t,MAXCOMPRESSEDPER)\
            tmpa -= MAXCOMPRESSEDPER;\
        }\
        if (tmpa > 0) {\
            WRITEF(uint8_t,FF)\
            WRITEF(int16_t,tmpa)\
        }\
    }

    FILE *f=fopen(afn,"wb");
    if (!f) return 0;

    WRITEF(int64_t,xlen)
    WRITEF(int64_t,ylen)

    for(int64_t y=0;y<ylen;y++) {
        int64_t anz=0;
        int32_t color=-1; // some INVALID color

        for(int64_t x=0;x<xlen;x++) {
            int32_t lokalf=getPunkt(x,y);

            if (color < 0) {
                color=lokalf;
                anz=1;
                continue;
            }

            if (lokalf == color) {
                anz++;
                continue;
            }

            SAVESTREAK_ANZF(anz,color)

            color=lokalf;
            anz=1;
        } // x

        if (anz > 0) {
            SAVESTREAK_ANZF(anz,color)
        }

        // EOL
        WRITEF(uint8_t,0)
        WRITEF(int16_t,COMBINEEOL)

    } // y

    // EOF
    WRITEF(uint8_t,0)
    WRITEF(int16_t,COMBINEEOF)

    fclose(f);

    return 1;
}

int8_t TiledPlane::loadBitmap(const char* afn) {
	// 8 bit Bitmap

	FILE *fbmp=fopen(afn,"rb");
	if (!fbmp) return 0;

	BYTE dummypuffer[4096];

	#define DUMMYREAD(NR) \
		fread(dummypuffer,NR,sizeof(BYTE),fbmp);

	DUMMYREAD(2)

	uint32_t off;
	uint32_t filelen;

	fread(&filelen,1,sizeof(filelen),fbmp);
	DUMMYREAD(4)
	fread(&off,1,sizeof(off),fbmp); // offset, ab da beginnen die PIXEL
	DUMMYREAD(4)

	uint32_t wx,wy;
	fread(&wx,sizeof(wx),1,fbmp);
	fread(&wy,sizeof(wy),1,fbmp);
	setlenxy(wx,wy);
	DUMMYREAD(2);
	uint16_t bitsperpixel;
	fread(&bitsperpixel,sizeof(bitsperpixel),1,fbmp);
	if (bitsperpixel != 8) {
		printf("Fehler. TiledPlane::read BitsPerPixel MUSS 8 sein.\n");
		exit(99);
	}

	DUMMYREAD(24)
	BYTE puffer[4];
	for(int32_t i=0;i<256;i++) {
		fread(puffer,4,sizeof(BYTE),fbmp);
		palette[i].B=puffer[0];
		palette[i].G=puffer[1];
		palette[i].R=puffer[2];
	}

	for(int32_t y=0;y<ylen;y++) {
		//cmpY[y]=planemgr.getMemory(xlen);
		fread(cmpY[y],xlen,sizeof(BYTE),fbmp);
	}

	fclose(fbmp);

	return 1;
}

void TiledPlane::saveFormatBase(const int32_t FORMAT,const char* afnbase) {
    char fn[2048];
    if (FORMAT == FORMAT_RAW) {
        sprintf(fn,"%s.data",afnbase);
        saveBitmap(fn);

        return;
    }

    sprintf(fn,"%s.cpr",afnbase);
    saveCompressed(fn);
}

int8_t TiledPlane::loadFormatBase(const int32_t FORMAT,const char* afnbase) {
    char fn[2048];
    if (FORMAT == FORMAT_RAW) {
        sprintf(fn,"%s.data",afnbase);
        return loadBitmap(fn);
    }

    sprintf(fn,"%s.cpr",afnbase);
    return loadCompressed(fn);
}

void TiledPlane::saveBitmap(const char* afn) {
	// the data is stored in the format of an 8 bit Bitmap

	FILE *fbmp=fopen(afn,"wb");
	write2(fbmp,66,77); // BM

	uint32_t off
		=		14 // FILEHeader
			+	40 // Bitmapheader
			+	256*4 // ColorPalette
		;

	// filelen will overflow if image width is too large
	// but external viewers
	// can display the image nonetheless
	uint32_t filelen
			=	off
			+	(ylen*xlen);
		;

	fwrite(&filelen,1,sizeof(filelen),fbmp);
	write4(fbmp,0,0,0,0);
	fwrite(&off,1,sizeof(off),fbmp); // offset, ab da beginnen die PIXEL
	write4(fbmp,40,0,0,0);

	uint32_t w = xlen;
	fwrite(&w,sizeof(w),1,fbmp);
	w = ylen;
	fwrite(&w,sizeof(w),1,fbmp);
	write2(fbmp,1,0);
	write2(fbmp,8,0); // 8 bits per pixel
	write4(fbmp,0,0,0,0);
	write4(fbmp,0,0,0,0);
	write4(fbmp,19,10,0,0);
	write4(fbmp,19,10,0,0);
	write4(fbmp,0,1,0,0); // number of palette entries
	write4(fbmp,0,0,0,0);
	BYTE puffer[4];
	for(int32_t i=0;i<256;i++) {
		puffer[0]=palette[i].B;
		puffer[1]=palette[i].G;
		puffer[2]=palette[i].R;
		puffer[3]=0;
		fwrite(puffer,4,sizeof(BYTE),fbmp);
	}

	for(int32_t y=0;y<ylen;y++) {
		fwrite(cmpY[y],xlen,sizeof(BYTE),fbmp);
	}

	fclose(fbmp);
}

void TiledPlane::lineHorizVert(const int32_t ax,const int32_t ay,const int32_t bx,const int32_t by,const BYTE awert) {
	if (!cmpY) return;

	if (ax == bx) {
		// vertical line
		int32_t y0,y1;
		if (ay < by) { y0=ay; y1=by; } else { y0=by; y1=ay; }
		for(int32_t y=y0;y<=y1;y++) {
			cmpY[y][ax]=awert;
		} // y
	} else if (ay == by) {
		// horizontal line
		int32_t x0,x1;
		if (ax < bx) { x0=ax; x1=bx; } else { x0=bx; x1=ax; }
		for(int32_t x=x0;x<=x1;x++) {
			cmpY[ay][x]=awert;
		} // y
	} else {
		// geht nicht
		printf("TiledPlane: Diagonal line not implemented.\n");
		return;
	}
}

void TiledPlane::fillrect(const int32_t ax,const int32_t ay,const int32_t bx,const int32_t by,const BYTE aw) {
	int32_t lx=ax,rx=bx;
	if (ax > bx) { lx=bx; rx=ax; }
	int32_t ly=ay,ry=by;
	if (ay > by) { ly=by; ry=ay; }

	for(int32_t y=ly;y<=ry;y++) {
		for(int32_t x=lx;x<=rx;x++) {
			cmpY[y][x]=aw;
		}
	}
}

void TiledPlane::setlenxy(const int64_t ax,const int64_t ay) {

    if (
        (cmpY) &&
        (ax == xlen) &&
        (ay == ylen)
    ) return; // used already allocated memory

	// memory itself is not deallocated here
	// as setlenxy is not called more than once

	if ((xlen>0)&&(cmpY)) {
		delete[] cmpY;
	}

	cmpY=new PBYTE[ay];
	xlen=ax;
	ylen=ay;
	for(int32_t y=0;y<ylen;y++) {
		cmpY[y]=planemgr.getMemory(xlen);
	}
}

TiledPlane::TiledPlane() {
	xlen=ylen=0;
	memused=0;
	cmpY=NULL;
}

TiledPlane::~TiledPlane() {
	// memory itself is deallocated when destroying the memory managaer struct
	if ((xlen>0)&&(cmpY)) {
		delete[] cmpY;
	}
}

void TiledPlane::setPaletteRGB(const int32_t pos,const BYTE ar,const BYTE ag,const BYTE ab) {
	if ((pos<0)||(pos>255)) return;
	palette[pos].R=ar;
	palette[pos].G=ag;
	palette[pos].B=ab;
}

void TiledPlane::setPunkt(const int32_t ax,const int32_t ay,const BYTE awert) {
	if (!cmpY) return;
	if (
		(ax<0) || (ax >= xlen) || (ay<0) || (ay>=ylen)
	) return;

	cmpY[ay][ax]=awert;
}

BYTE TiledPlane::getPunkt(const int32_t ax,const int32_t ay) {
	if (
		(!cmpY) ||
		(ax<0) || (ax >= xlen) || (ay<0) || (ay>=ylen)
	) {
		LOGMSG("Error. getPunkt for a pixel not in memory\n");
		exit(99);
	}

	return cmpY[ay][ax];
}

// struct ArrayDDByteManager

// alloctes memory in chunks for the TiledPlane
// objects to keep memory fragmentation low

ArrayByteMgr::ArrayByteMgr() {
	current=NULL;
	allocatedIdx=0;
	freeFromIdx=-1;
	anzptr=0;
	double d=CHARMAPCHUNKSIZE; d /= sizeof(BYTE);
	allocatePerBlockIdx=(int)floor(d);
}

void ArrayByteMgr::FreeAll(void) {
	for(int32_t i=0;i<anzptr;i++) {
		delete[] ptr[i];
	}
	anzptr=0;
	current=NULL;
}

ArrayByteMgr::~ArrayByteMgr() {
	FreeAll();
}

PBYTE ArrayByteMgr::getMemory(const int32_t aanz) {
	if (anzptr >= (MAXPTR-8)) {
		LOGMSG("Error, memory. ArrayByteMgr/1");
		exit(99);
	}
	if (
		(!current) ||
		((freeFromIdx + aanz + 2) >= allocatedIdx)
	) {
		ptr[anzptr]=current=new BYTE[allocatePerBlockIdx];
		anzptr++;
		if (!current) {
			LOGMSG("Error, memory. ArrayByteMgr/2");
			exit(99);
		}
		freeFromIdx=0;
		allocatedIdx=allocatePerBlockIdx;
	}

	PBYTE p=&current[freeFromIdx];
	freeFromIdx += aanz;
	return p;
}

void setFunc(void) {
    SYMMETRY_TO_XAXIS=1;
    SYMMETRY_TO_YAXIS=0;
    RANGEM0=-2;
    RANGEM1=2;
}

char* getFnBase_TSXY(char* erg,const int64_t ALEN,const int64_t ax,const int64_t ay) {
    sprintf(erg,"_M%010I64d_y%010I64d_x%010I64d",
        ALEN,ay,ax);
    return erg;
}

int8_t cust_escaping_C_affine(CustComplex& c) {
    // return: > 0 => escaping
    // = 0 => non-escpaing (so far)
    // < 0 => error

    // affine arithmetics
    // convert c to an affine form
    // follow the orbit (no 2-sqzare-cover fast test)
    // clear used noise symboles as those are
    // only valid here temporarily
    int64_t restore_nextfreesymbol=globalafpath.nextfreesymbol;

    #define RETURN(WW) \
    {\
        globalafpath.nextfreesymbol=restore_nextfreesymbol;\
        return (WW);\
    }

    int32_t error=0;

    // convert into affine-form
    AffineComplex afc;
    error += setFromCustInt_TAP(afc.real,c.re,globalafpath);
    error += setFromCustInt_TAP(afc.imag,c.im,globalafpath);
    if (error != 0) { RETURN(-1); } ; // error

    // critical point as start
    AffineComplex afz;
    afz.setToZero();
    afz.real.setBelongs(&globalafpath);
    afz.imag.setBelongs(&globalafpath);

    // currently no fast-2-square-covering-.test
    for(int32_t k=0;k<MAXIT;k++) {

        // fast pre-test: check if central value is outside 2-square
        if (
            (afz.real.xi[0] < COMPLETEM0cust.left ) ||
            (afz.real.xi[0] > COMPLETEM1cust.right) ||
            (afz.imag.xi[0] < COMPLETEM0cust.left ) ||
            (afz.imag.xi[0] > COMPLETEM1cust.right)
        ) {
            // if so => convert afz to CustComplex
            // and check whetehr the whole rectangle is
            // putside => escaping
            CustComplex z;
            error += convertToCustInt_TA(z.re,afz.real);
            error += convertToCustInt_TA(z.im,afz.imag);
            if (error != 0) { RETURN(-1); } // error, cannot judge

            // does this fully lie outside the 2-square
            if (
                (z.re.right < COMPLETEM0cust.left ) ||
                (z.re.left  > COMPLETEM1cust.right) ||
                (z.im.right < COMPLETEM0cust.left ) ||
                (z.im.left  > COMPLETEM1cust.right)
            ) {
                RETURN(1) // ESCAPING
            }

            // not fully poutside, continue iteration

        } // central value outside 2-square

        // iteration
        // afz <- afz*afz+c
        AffineComplex t1;
        error += afcplxMul_TAB_const(t1,afz,afz);
        // afz <- t1+c
        error += afcplxAdd_TAB_const(afz,t1,afc);

        if (error != 0) { RETURN(-1); }

        // TODO: if number of noise symbols gets too large
        // (define this), convert the affineform into an interval
        // and back, so starting with only one noise symbol
        // in sum, this still should lead to a much smaller
        // orbit point but might enhance chances of detecting
        // escaping

        if (
            (afz.real.highestusedsymbol > AFFINEFOLDBACK) ||
            (afz.imag.highestusedsymbol > AFFINEFOLDBACK)
        ) {
            // as real,imag can share noise symbols, BOTH have to be
            // converted back and forth
            CustComplex ztmp;
            error += convertToCustInt_TA(ztmp.re,afz.real);
            error += convertToCustInt_TA(ztmp.im,afz.imag);
            if (error != 0) { RETURN(-1); } // error

            // now the affinePath will be reset to the
            // startin value
            globalafpath.nextfreesymbol=restore_nextfreesymbol;
            // now afz is INVALID

            // and ztmp will be converted back to an affine form
            // and stored into afz
            error += setFromCustInt_TAP(afz.real,ztmp.re,globalafpath);
            error += setFromCustInt_TAP(afz.imag,ztmp.im,globalafpath);

            ctrfoldback++;

        } // folding event

    } // k

    if (error != 0) {
        RETURN(-1);
    }

    RETURN(0) // non-escaping

}

int8_t cust_escaping_C(CustComplex& c) {
    // return: 0 => non-escaping
    // 1 => escaping
    // < 0 => error

    CustComplex z;
    z.setToZero(); // critical point

    int32_t error=0;

    for(int32_t k=0;k<MAXIT;k++) {
        // escaping =
        if (z.re.right < COMPLETEM0cust.left) return 1;
        if (z.re.left  > COMPLETEM1cust.right) return 1;
        if (z.im.right < COMPLETEM0cust.left) return 1;
        if (z.im.left  > COMPLETEM1cust.right) return 1;

        // if z contains 0 and the whole 2-square
        // it can never be judged as escaping
        // return at 0
        // see: https://fractalforums.org/index.php?topic=2725.msg21221#msg21221

        if (
            (z.re.left  < COMPLETEM0cust.left ) &&
            (z.re.right > COMPLETEM1cust.right) &&
            (z.im.left  < COMPLETEM0cust.left ) &&
            (z.im.right > COMPLETEM1cust.right)
        ) return 0;

        // iteration
        CustComplex z2,zcopy;
        zcopy.copyFrom(z);
        error += custCplxMul_ZAB(z2,z,zcopy);
        error += custCplxAdd_ZAB(z,z2,c);
        if (error != 0) return -1;

    } // k

    if (error != 0) return -1;

    return 0; // not escaping
}

int32_t coordinate_TABC(
    CustInterval& erg,
    const int64_t A,
    CustInterval& B,
    CustInterval& C
) {
    // computes: A*B+C
    int32_t error=0;
    CustInterval tA;
    char tt[256];
    sprintf(tt,"%I64d",A);

    error += tA.setStr(tt);
    // erg=tA*B+c
    CustInterval t2;
    error += custMul_ZAB(t2,tA,B);
    // erg=t2+c
    error += custAdd_ZAB(erg,t2,C);

    return error;

}

// fast analysis : take many pixels together and analyze those
// (only used at creation level)
void analyzeSpecial_BXYA(
    TiledPlane& md,
    const int64_t filex,
    const int64_t filey,
    const int64_t DELTA
) {
    int64_t noch=1,noch0=(md.ylen >> 2) / DELTA;

    CustComplex c;
    int32_t error=0;

    for(int64_t y=0;y<md.ylen;y+=DELTA) {
        if ((--noch)<=0) {
            printf(".");
            noch=noch0;
        }

        // y*scaleRangePerPixelMia + COMPLETEM0ia
        CustInterval im1,im2;

        error=0;
        error += coordinate_TABC(im1,filey+y,scaleRangePerPixelMcust,COMPLETEM0cust);
        error += coordinate_TABC(im2,filey+y+DELTA,scaleRangePerPixelMcust,COMPLETEM0cust);

        c.im.left=im1.left;
        c.im.right=im2.right;

        if (error != 0) continue; // no return error as pixels are not changted

        for(int64_t x=0;x<md.xlen;x+=DELTA) {
            if (md.getPunkt(x,y) != MTILE_GRAY) continue;

            md.setPunkt(x,y,MTILE_GRAY_DONE);
            // will be set to WHITE if escaping

            CustInterval re1,re2;

            error=0;
            error += coordinate_TABC(re1,filex+x,scaleRangePerPixelMcust,COMPLETEM0cust);
            error += coordinate_TABC(re2,filex+x+DELTA,scaleRangePerPixelMcust,COMPLETEM0cust);

            c.re.left=re1.left;
            c.re.right=re2.right;

            if (error != 0) continue; // no return error as pixels are not changted

            int8_t esc=cust_escaping_C(c);

            if (
                (USEAFFINE > 0) &&
                (esc <= 0)
            ) {
                // if classic IA cannot detect escaping, try affine
                esc=cust_escaping_C_affine(c);
            }

            // < 0 => error, will be trated as NO CHANGE to image
            // > 0 => reliably escaping escaping
            if (esc > 0) {
                md.fillrect(x,y,x+DELTA-1,y+DELTA-1,MTILE_WHITE);
            }

        } // x
    } // y

    // remove everything non-black, non-white and set
    // to gray
    for(int64_t y=0;y<md.ylen;y++) {
        for(int64_t x=0;x<md.xlen;x++) {
            int32_t f=md.getPunkt(x,y);

            if (
                (f != MTILE_WHITE) &&
                (f != MTILE_BLACK)
            ) {
                md.setPunkt(x,y,MTILE_GRAY);
            }

        } // x
    } // y


}

void analyzeFileToCheck_cust_BWXYN_T(
    TiledPlane& md,
    ScreenRect64& inrect,
    const int64_t filex,
    const int64_t filey,
    const int64_t ctr,
    ScreenRect64& newwhite
) {
    newwhite.x0=-1; // none found

    CustComplex c;
    int32_t error=0;

    for(int64_t y=inrect.y0;y<=inrect.y1;y++) {
        // y*scaleRangePerPixelMia + COMPLETEM0ia
        CustInterval iay,t1,im1,im2;
        char tt[256];
        sprintf(tt,"%I64d",filey+y);
        iay.setStr(tt);

        error=0;

        error += custMul_ZAB(t1,iay,scaleRangePerPixelMcust);
        error += custAdd_ZAB(im1,t1,COMPLETEM0cust);
        error += custAdd_ZAB(im2,im1,scaleRangePerPixelMcust);
        c.im.left=im1.left;
        c.im.right=im2.right;

        if (error != 0) continue; // no return error as pixels are not changted

        for(int64_t x=inrect.x0;x<=inrect.x1;x++) {
            if (md.getPunkt(x,y) != MTILE_CHECKFORWHITE) continue;

            md.setPunkt(x,y,MTILE_GRAY_DONE);
            // will be set to WHITE if escaping

            CustInterval iax,t2,re1,re2;
            char tt[256];
            sprintf(tt,"%I64d",filex+x);
            iax.setStr(tt);

            error += custMul_ZAB(t2,iax,scaleRangePerPixelMcust);
            error += custAdd_ZAB(re1,t2,COMPLETEM0cust);
            error += custAdd_ZAB(re2,re1,scaleRangePerPixelMcust);
            c.re.left=re1.left;
            c.re.right=re2.right;

            if (error != 0) continue; // no return error as pixels are not changted

            int8_t esc=cust_escaping_C(c);
            // < 0 => error, will be trated as NO CHANGE to image
            // > 0 => reliably escaping escaping

            if (
                (USEAFFINE > 0) &&
                (esc <= 0)
            ) {
                // try affine arithmetics
            if (
                (USEAFFINE > 0) &&
                (esc <= 0)
            ) {
                // if classic IA cannot detect escaping, try affine
                esc=cust_escaping_C_affine(c);
            }

            } //

            if (esc > 0) {
                md.setPunkt(x,y,MTILE_WHITE);
                if (newwhite.x0 < 0) {
                    newwhite.x0=newwhite.x1=x;
                    newwhite.y0=newwhite.y1=y;
                } else {
                    if (x < newwhite.x0) newwhite.x0=x;
                    if (x > newwhite.x1) newwhite.x1=x;
                    if (y < newwhite.y0) newwhite.y0=y;
                    if (y > newwhite.y1) newwhite.y1=y;
                }
            }

        } // x
    } // y

}

void rectFill_BA(
    TiledPlane& md,
    const int64_t GJ
) {
    // md does NEVER contain the full Mset, so
    // any exterior-cirumferenced region is itself
    // exterior by simple connectivity

    for(int64_t ygrid=GJ;ygrid<md.ylen;ygrid+=GJ) {

        for(int64_t xgrid=GJ;xgrid<md.xlen;xgrid+=GJ) {

            int8_t white=1;
            for(int64_t g=0;g<=GJ;g++) {
                if (
                    (md.getPunkt(xgrid-g ,ygrid   ) != MTILE_WHITE) ||
                    (md.getPunkt(xgrid-g ,ygrid-GJ) != MTILE_WHITE) ||
                    (md.getPunkt(xgrid   ,ygrid-g ) != MTILE_WHITE) ||
                    (md.getPunkt(xgrid-GJ,ygrid-g ) != MTILE_WHITE)
                ) {
                    white=0;
                    break;
                }

            } // g

            if (white <= 0) continue; // next

            int8_t gn=0;
            for(int64_t dy=1;dy<GJ;dy++) {
                for(int64_t dx=1;dx<GJ;dx++) {
                    if (md.getPunkt(xgrid-dx,ygrid-dy) != MTILE_WHITE) {
                        gn=1;
                        break;
                    }
                } // dx

                if (gn > 0) break;
            } // dy

            if (gn <= 0) continue; // nothing to fill

            md.fillrect(xgrid,ygrid,xgrid-GJ,ygrid-GJ,MTILE_WHITE);

        } // xgtrid

    } // ygird

}

void analyzeFile(
    TiledPlane& md,
    const int64_t filex,
    const int64_t filey,
    const char* fnbase,
    TiledPlane& mdhelp
) {

    // first analyze by P123 in cust
    if (P123ON > 0) {
        printf(" p123 ... ");
        int64_t noch=1,noch0=md.ylen >> 3;
        CustComplex c;

        for(int64_t y=0;y<md.ylen;y++) {
            if ((--noch)<=0) {
                printf(".");
                noch=noch0;
            }
            // y*scaleRangePerPixelMia + COMPLETEM0ia
            CustInterval iay,t1,im1,im2;
            char tt[256];
            sprintf(tt,"%I64d",filey+y);
            iay.setStr(tt);

            int32_t error=0;

            error += custMul_ZAB(t1,iay,scaleRangePerPixelMcust);
            error += custAdd_ZAB(im1,t1,COMPLETEM0cust);
            error += custAdd_ZAB(im2,im1,scaleRangePerPixelMcust);
            c.im.left=im1.left;
            c.im.right=im2.right;

            if (error != 0) continue; // no return error as pixels are not changted

            for(int64_t x=0;x<md.xlen;x++) {
                if (md.getPunkt(x,y) != MTILE_GRAY) continue;

                CustInterval iax,t2,re1,re2;
                char tt[256];
                sprintf(tt,"%I64d",filex+x);
                iax.setStr(tt);

                error += custMul_ZAB(t2,iax,scaleRangePerPixelMcust);
                error += custAdd_ZAB(re1,t2,COMPLETEM0cust);
                error += custAdd_ZAB(re2,re1,scaleRangePerPixelMcust);
                c.re.left=re1.left;
                c.re.right=re2.right;

                if (error != 0) continue; // no return error as pixels are not changted

                // check via reliable p1,p2 to set some
                // gray pixels to interior
                // (not necessary for exterior test, but
                // reduces number of pixels for a subsequentz
                // guide_ps orbit analysis and reduces torage

                int8_t set123=0;
                if (set123 <= 0) if (predP2_z2c(c) == INP12_INSIDE) set123=1;
                if (set123 <= 0) if (predP1_z2c(c) == INP12_INSIDE) set123=1;

                if (set123 > 0) {
                    md.setPunkt(x,y,MTILE_BLACK);
                    continue;
                }
            } // x
        } // y

    } // P123

    // mark in a grid some pixels to be
    // reliably analyzed: very sparse and some
    // on the image border
    // this givces NEW exterior in case a fully gray image is
    // analyzed, or gives some new white blobs that are
    // (NOT YET) connected to the rest of the Mset exterior

    for(int64_t y=0;y<md.ylen;y+=GRIDWIDTH) {
        if (md.getPunkt(0,y) == MTILE_GRAY) {
            md.setPunkt(0,y,MTILE_CHECKFORWHITE);
        }

        if (md.getPunkt(md.xlen-1,y) == MTILE_GRAY) {
            md.setPunkt(md.xlen-1,y,MTILE_CHECKFORWHITE);
        }
    } // y

    for(int64_t x=0;x<md.xlen;x+=GRIDWIDTH) {
        if (md.getPunkt(x,0) == MTILE_GRAY) {
            md.setPunkt(x,0,MTILE_CHECKFORWHITE);
        }

        if (md.getPunkt(x,md.ylen-1) == MTILE_GRAY) {
            md.setPunkt(x,md.ylen-1,MTILE_CHECKFORWHITE);
        }
    } // y

    int8_t changed=0;
    int64_t ctrcheck=0;
    int8_t first=1;

    ScreenRect64 fulltile;

    fulltile.x0=0;
    fulltile.x1=TILELEN-1;
    fulltile.y0=0;
    fulltile.y1=TILELEN-1;

    ScreenRect64 newhite;
    changed=1;
    // 1st method: only analyze vertical and horizobntal
    // grid lines of a specific disrtance

    #define MARKIF(XX,YY) \
    {\
        int32_t f=md.getPunkt(XX,YY);\
        \
        if (f == MTILE_GRAY) {\
            int8_t wn=0;\
            for(int64_t dy=-1;dy<=1;dy++) {\
                int64_t y2=(YY)+dy;\
                if ( (y2<0) || (y2 >= md.ylen) ) continue;\
                \
                for(int64_t dx=-1;dx<=1;dx++) {\
                    int64_t x2=(XX)+dx;\
                    if ( (x2<0) || (x2 >= md.xlen) ) continue;\
                    \
                    if (md.getPunkt(x2,y2) == MTILE_WHITE) {\
                        wn=1;\
                        break;\
                    }\
                \
                }\
                \
                if (wn > 0) break;\
                \
            }\
            \
            if (wn > 0) {\
                ctrcheck++;\
                md.setPunkt(XX,YY,MTILE_CHECKFORWHITE);\
            }\
        } else if (f == MTILE_CHECKFORWHITE) {\
                    ctrcheck++;\
        }\
    }

    if (USEGRIDFILL > 0) {
        printf(" grid ... ");
        // check grid points, if two adjacent are at least one
        // white => mark everything tray on that line as
        // TO CHECK
        // then do oner flood-fill attempt

        for(int32_t gloop=0;gloop<GRIDLOOPS;gloop++) {
            ctrcheck=0;

            for(int64_t y=GRIDFILLJUMP;y<md.ylen;y+=GRIDFILLJUMP) {

                for(int64_t x=GRIDFILLJUMP;x<md.xlen;x+=GRIDFILLJUMP) {

                    if (
                        (md.getPunkt(x,y) == MTILE_WHITE) ||
                        (md.getPunkt(x-GRIDFILLJUMP,y) == MTILE_WHITE)
                    ) {
                        // horizontal left
                        for(int64_t x2=(x-GRIDFILLJUMP);x2<=x;x2++) {
                            int32_t f=md.getPunkt(x2,y);

                            if (f == MTILE_CHECKFORWHITE) ctrcheck++;
                            else if (f == MTILE_GRAY) {
                                md.setPunkt(x2,y,MTILE_CHECKFORWHITE);
                                ctrcheck++;
                            }
                        } // x2
                    }

                    if (
                        (md.getPunkt(x,y) == MTILE_WHITE) ||
                        (md.getPunkt(x,y-GRIDFILLJUMP) == MTILE_WHITE)
                    ) {
                        // vertical donw
                        for(int64_t y2=(y-GRIDFILLJUMP);y2<=y;y2++) {
                            int32_t f=md.getPunkt(x,y2);

                            if (f == MTILE_CHECKFORWHITE) ctrcheck++;
                            else if (f == MTILE_GRAY) {
                                md.setPunkt(x,y2,MTILE_CHECKFORWHITE);
                                ctrcheck++;
                            }

                        } // y2
                    }

                } // x

            } // y

            // now connect two CHECKFROWHITE grid points
            for(int64_t y=GRIDFILLJUMP;y<md.ylen;y+=GRIDFILLJUMP) {

                for(int64_t x=GRIDFILLJUMP;x<md.xlen;x+=GRIDFILLJUMP) {

                    if (
                        (md.getPunkt(x,y) == MTILE_CHECKFORWHITE) &&
                        (md.getPunkt(x-GRIDFILLJUMP,y) == MTILE_CHECKFORWHITE)
                    ) {
                        // horizontal left
                        for(int64_t x2=(x-GRIDFILLJUMP);x2<=x;x2++) {
                            int32_t f=md.getPunkt(x2,y);

                            if (f == MTILE_CHECKFORWHITE) ctrcheck++;
                            else if (f == MTILE_GRAY) {
                                md.setPunkt(x2,y,MTILE_CHECKFORWHITE);
                                ctrcheck++;
                            }
                        } // x2
                    }

                    if (
                        (md.getPunkt(x,y) == MTILE_CHECKFORWHITE) &&
                        (md.getPunkt(x,y-GRIDFILLJUMP) == MTILE_CHECKFORWHITE)
                    ) {
                        // vertical donw
                        for(int64_t y2=(y-GRIDFILLJUMP);y2<=y;y2++) {
                            int32_t f=md.getPunkt(x,y2);

                            if (f == MTILE_CHECKFORWHITE) ctrcheck++;
                            else if (f == MTILE_GRAY) {
                                md.setPunkt(x,y2,MTILE_CHECKFORWHITE);
                                ctrcheck++;
                            }

                        } // y2
                    }

                } // x

            } // y

            if (ctrcheck > 0) {
                printf("%I64d ",ctrcheck);

                char tt[1024];

                /*
                sprintf(tt,"_gctr%06i_before_grid_analysis.bmp",gctr);
                md.saveBitmap(tt);
                gctr++;
                */


                ScreenRect64 check;
                check.copyFrom(fulltile);

                analyzeFileToCheck_cust_BWXYN_T(md,check,filex,filey,ctrcheck,newhite);

                rectFill_BA(md,GRIDFILLJUMP);

            } // analyze

        } // gloop

    } // use gridfill

    ScreenRect64 toa; // rect to check for gray adj white
    toa.copyFrom(fulltile);

    changed=1;
    int64_t ctrkum=0,ctrportions=0;
    int64_t PRINT=(int64_t)1 << 20;

    while (changed > 0) {
        changed=0;

        // go over image and look for TOCHECK
        // and GRAY with adjacent WHITE
        ctrcheck=0;
        newhite.x0=-1; // empty

        toa.copyFrom(fulltile);

        ScreenRect64 check; // where to look
        check.x0=TILELEN;
        check.x1=0;
        check.y0=TILELEN;
        check.y1=0;

        for(int64_t y=toa.y0;y<=toa.y1;y++) {
            for(int64_t x=toa.x0;x<=toa.x1;x++) {
                int32_t f=md.getPunkt(x,y);

                if (f == MTILE_GRAY) {
                    // white neighbour ?
                    int8_t wn=0;
                    for(int64_t dy=-1;dy<=1;dy++) {
                        int64_t y2=y+dy;
                        if ( (y2<0) || (y2 >= md.ylen) ) continue;

                        for(int64_t dx=-1;dx<=1;dx++) {
                            int64_t x2=x+dx;
                            if ( (x2<0) || (x2 >= md.xlen) ) continue;

                            if (md.getPunkt(x2,y2) == MTILE_WHITE) {
                                wn=1;
                                break;
                            }

                        } // dx

                        if (wn > 0) break;

                    } // dy

                    if (wn > 0) {
                        ctrcheck++;
                        md.setPunkt(x,y,MTILE_CHECKFORWHITE);

                        if (x < check.x0) check.x0=x;
                        if (x > check.x1) check.x1=x;
                        if (y < check.y0) check.y0=y;
                        if (y > check.y1) check.y1=y;
                    }

                } else if (f == MTILE_CHECKFORWHITE) {
                    ctrcheck++;
                    if (x < check.x0) check.x0=x;
                    if (x > check.x1) check.x1=x;
                    if (y < check.y0) check.y0=y;
                    if (y > check.y1) check.y1=y;
                }
            } // x
        } // y

        if (ctrcheck <= 0) break;

        if (first > 0) {
            printf(" check %I64d ",ctrcheck);
            first=0;
        } else {
            ctrkum += ctrcheck;
            ctrportions++;

            if (ctrkum > PRINT) {
                printf(" %I64d(/%I64d) ",ctrkum,ctrportions);
                if (ctrportions > 100) {
                    if (PRINT > 10000) PRINT >>= 1;
                }
                ctrkum=0;
                ctrportions=0;
            }
        }

        if (
            (SKIPIFLESS > 0) &&
            (ctrcheck > 0)
        ) {
            // if too few => next file
            if (ctrcheck < SKIPIFLESS) {
                changed=0;
                printf("  skipped as below %I64d ",SKIPIFLESS);

                break;
            }
        }

        if (ctrcheck > 0) {

            analyzeFileToCheck_cust_BWXYN_T(md,check,filex,filey,ctrcheck,newhite);
            changed=1;

            if (newhite.x0 >= 0) {
                // new white set
                toa.copyFrom(newhite);
                // extand by 1 as GRAY will be analyzed that has
                // a white neighbour
                if (toa.x0 > 0) toa.x0--;
                if (toa.x1 < (TILELEN-1)) toa.x1=TILELEN-1;
                if (toa.y0 > 0) toa.y0--;
                if (toa.y1 < (TILELEN-1)) toa.y1=TILELEN-1;
            } else {
                // whole TILE
                toa.copyFrom(fulltile);
            }

        }

    } // while

}

void setDict_file_XYF(
    const int64_t filex,
    const int64_t filey,
    const int32_t af
) {
    if (USEDICTIONARY <= 0) return;

    int64_t dictx=filex >> TILELEVEL;
    int64_t dicty=filey >> TILELEVEL;

    if (
        (dictx >= dict.xlen) ||
        (dicty >= dict.ylen)
    ) {
        LOGMSG3("\nimplementation error. dict access outside (%I64d,%I64d)\n",dictx,dicty);
        exit(99);
    }

    dict.setPunkt(dictx,dicty,af);

    //printf("\nsaving dictionary ... ");
    saveDictionaryAndTwice();
    //printf("done. save to terminate program after next file start\n");

}

int32_t getDictColor_fileXY(
    const int64_t filex,
    const int64_t filey
) {
    if (USEDICTIONARY <= 0) return DICT_GRAY;

    int64_t dictx=filex >> TILELEVEL;
    int64_t dicty=filey >> TILELEVEL;

    if (
        (dictx >= dict.xlen) ||
        (dicty >= dict.ylen)
    ) {
        LOGMSG3("\nimplementation error. dict access get outside (%I64d,%I64d)\n",dictx,dicty);
        exit(99);
    }

    return dict.getPunkt(dictx,dicty);

}

void calcFile_TB(
    Calcs& ergs,
    TiledPlane& md
) {
    ergs.setToZero();

    for(int64_t y=0;y<md.ylen;y++) {
        for(int64_t x=0;x<md.xlen;x++) {
            int32_t f=md.getPunkt(x,y);

            if (f == MTILE_WHITE) ergs.ctrw++;
            else if (f == MTILE_BLACK) ergs.ctrb++;
            else if (f == MTILE_NUMERICALCYCLE) ergs.ctrnum++;
            else ergs.ctrg++;
        } // x
    } // y

}

// main function
void analyze(void) {
    ctrexterior=0;

    int64_t YEND=SCREENWIDTHM >> 1;
    TiledPlane md,mdhelp;
    md.setlenxy(TILELEN,TILELEN);
    setPaletteTo(md);
    mdhelp.setlenxy(TILELEN,TILELEN);
    setPaletteTo(mdhelp);

    for(int64_t filey=0;filey<YEND;filey+=TILELEN) {
        for(int64_t filex=0;filex<SCREENWIDTHM;filex+=TILELEN) {

            int32_t dictf=getDictColor_fileXY(filex,filey);

            if (
                (dictf == DICT_BLACK) ||
                (dictf == DICT_NUMERICALCYCLE)
            ) continue; // no analysis, no saving

            if (dictf == DICT_WHITE) {
                ctrexterior += (int64_t)TILELEN*TILELEN;
                continue; // no analysuis, no twice

            }

            char fn[2048];
            getFnBase_TSXY(fn,SCREENWIDTHM,filex,filey);
            if (md.loadFormatBase(LOADFORMAT,fn) <= 0) {
                // if dictionary entry is GRAY
                // => error, there must be a file
                if (dictf == DICT_GRAY) {
                    LOGMSG("\nconsistency error. dictionary entry is gray, but no file exists.\n");
                    LOGMSG2("%s\N",fn);
                    exit(99);
                }

                continue;
            }

            printf("\n[y%I64d]",filey);
            printf(" <x%I64d> ",filex);

            // file exists
            analyzeFile(md,filex,filey,fn,mdhelp);

            Calcs erg;
            calcFile_TB(erg,md);
            ctrexterior += erg.ctrw;

            if (
                (erg.ctrw == tilepixelsint) &&
                (erg.ctrg == 0) &&
                (erg.ctrnum == 0) &&
                (erg.ctrb == 0)
            ) {
                // full white
                setDict_file_XYF(filex,filey,DICT_WHITE);
                continue;
            } else
            if (
                (erg.ctrnum == tilepixelsint) &&
                (erg.ctrg == 0) &&
                (erg.ctrw == 0) &&
                (erg.ctrb == 0)
            ) {
                // if all are NUMERICALCYCLE
                // discard as it is maybe the inteiror of
                // a hyperbolic component => save stoarge

                setDict_file_XYF(filex,filey,DICT_NUMERICALCYCLE);

                continue;
            }

            // not used so far

            if (
                (erg.ctrw == 0) &&
                (erg.ctrg == 0) &&
                (erg.ctrnum == 0) &&
                (erg.ctrb == tilepixelsint)
            ) {
                // full white
                setDict_file_XYF(filex,filey,DICT_BLACK);
                continue;
            }

            // and store again
            md.saveFormatBase(SAVEFORMAT,fn);

            if (SAVETWICE > 0) {
                saveTwice_BFXY(md,SAVEFORMAT,filex,filey);
            }

        } // x
    } // y
}

void helppage(void) {
    printf("\n\ncommand line parameters:\n------------------------\n\n");

    printf("working example: start the following two commands:\n");
    printf("    upperbound-rc anew lenm=12 lent=11 load=raw save=raw\n");
    printf("    upperbound-rc lenm=12 lent=11 load=raw save=raw\n");
    printf("\nThe stored files in format .DATA are 8-bit bitmaps and can\n");
    printf("be viewed after temporary renaming them to .BMP\n");
    printf("\nThe files dict_XXX store the dictionary, a 2^LENT-fold reliably\n");
    printf("downscaled version to only refine files for higher levels if the\n");
    printf("dictionary entry is gray.\n\n");

    printf("ANEW      : mandatory first step and parameter to start a new\n");
    printf("            computation in an otherwise empty directory.\n\n");

    printf("LENM=n    : the refinement level. Denotes a full image of\n");
    printf("            2^n x 2^n pixels covering the complex 2-square.\n");
    printf("            Only the lower half is analyzed and can be\n");
    printf("            mirrored due to symmetry.\n");

    printf("LENT=m    : the total image's lower half is partitioned into\n");
    printf("            files of size 2^m x 2^m pixels. If not specified,\n");
    printf("            the value 15 is set.\n");

    printf("P123      : flag to emply a test for membership of a pixel to\n");
    printf("            the period-1 or period-2 component of the quadratic\n");
    printf("            Mandelbrot set. Used only to reduce disk space usage\n");
    printf("            as fully interior chunks are not stored as files,\n");
    printf("            but solely in the dictionary.\n");

    printf("LOAD=RAW  : specifies the format of data is RAW, i.e. 8-bit bitmap\n");
    printf("            with file extension .DATA (and can be temporarily renamed).\n");
    printf("            If not specified, loads data as compress format .CPR\n");

    printf("SAVE=RAW  : stores data as bitmaps; if not specified as compressed..\n");

    printf("SKIPLESS=N: repeated sweeps over an image stop, if fewer\n");
    printf("            than N gray, unjudged pixels have an exterior\n");
    printf("            neighbour.\n");

    printf("MAXIT=n   : maximum iteration count for escape time analysis.\n");
    printf("            If not specified, the number \"3*LENM\" is used.\n");

}

void ScreenRect64::ausgabe(FILE* f) {
    fprintf(f,"[%I64d..%I64d,%I64d..%I64d]",
        x0,x1,y0,y1);
}

void ScreenRect64::copyFrom(ScreenRect64& A) {
    x0=A.x0;
    x1=A.x1;
    y0=A.y0;
    y1=A.y1;
}

int32_t sum_int32t(
	const int32_t a,
	const int32_t b
) {
	int64_t sum=(int64_t)a + (int64_t)b;
	if (
		(sum < (-INT32MAX)) ||
		(sum >   INT32MAX)
	) {
		LOGMSG("\nError. Overflow by int32 addition\n");
		exit(99);
	}

	return (int32_t)sum; // lower half
}

char* getDictFn_TW(char* erg,const int64_t alevel) {
    if (alevel > 50) {
        printf("\nerror. dict level wrong at save\n");
        exit(99);
    }

    sprintf(erg,"_dictL%03I64d.data",alevel);

    return erg;
}

void createanew(void) {
    // at exact one level higher than tilelevel
    int64_t YE=SCREENWIDTHM >> 1;

    printf(" ... ");

    TiledPlane md;
    md.setlenxy(TILELEN,TILELEN);
    setPaletteTo(md);
    md.fillrect(0,0,md.xlen-1,md.ylen-1,MTILE_GRAY);

    int32_t ST=SAVETWICE;
    SAVETWICE=0;

    for(int64_t filey=0;filey<YE;filey+=TILELEN) {

        for(int64_t filex=0;filex<SCREENWIDTHM;filex+=TILELEN) {
            char fn[2048];
            getFnBase_TSXY(fn,SCREENWIDTHM,filex,filey);

            md.fillrect(0,0,md.xlen-1,md.ylen-1,MTILE_GRAY);
            // to analyze very large portions for the fast
            // escaping parts
            analyzeSpecial_BXYA(md,filex,filey,64);

            md.saveFormatBase(SAVEFORMAT,fn);

        } // filex

    } // filey

    SAVETWICE=ST;

    // dictionary
    dict.setlenxy(DICTLEN,DICTLEN);
    setPaletteTo(dict);
    dict.fillrect(0,0,dict.xlen-1,dict.ylen-1,DICT_GRAY);
    char tt[2048];
    getDictFn_TW(tt,DICTLEVEL);
    dict.saveBitmap(tt);

    printf("done\n\n");
}

void loadDictionary(void) {
    char fn[2048];
    getDictFn_TW(fn,DICTLEVEL);
    dict.setlenxy(DICTLEN,DICTLEN);
    setPaletteTo(dict);

    if (dict.loadBitmap(fn) <= 0) {
        // no longer usable
        USEDICTIONARY=0;
        return;
    }

    USEDICTIONARY=1;
}

void saveDictionaryAndTwice(void) {
    char fn1[2048],fn2[2048];
    int64_t D2=DICTLEN << 1;
    getDictFn_TW(fn1,DICTLEVEL);
    dict.saveBitmap(fn1);

    // twice image
    TiledPlane dict2;
    dict2.setlenxy(D2,D2);
    setPaletteTo(dict2);

    int64_t tgty=-2;
    for(int64_t srcy=0;srcy<dict.ylen;srcy++) {
        tgty += 2;

        int64_t tgtx=-2;
        for(int64_t srcx=0;srcx<dict.xlen;srcx++) {
            tgtx += 2;
            int32_t f=dict.getPunkt(srcx,srcy);

            dict2.setPunkt(tgtx  ,tgty  ,f);
            dict2.setPunkt(tgtx+1,tgty  ,f);
            dict2.setPunkt(tgtx  ,tgty+1,f);
            dict2.setPunkt(tgtx+1,tgty+1,f);

        } // x

    } // y

    getDictFn_TW(fn2,DICTLEVEL+1);
    dict2.saveBitmap(fn2);

}

void Calcs::setToZero(void) {
    ctrw=ctrg=ctrb=ctrnum=0;
}

int8_t predP1_z2c(CustComplex& A) {
    // return: constant of INP12_XXX

    /*
		https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/Mandelbrot_set/mandelbrot

		q=((x-0.25)^2+y^2)
		q*(q+(x-0.25)) < 0.25*y^2
		((x-0.25)^2+y^2)*(((x-0.25)^2+y^2)+(x-1/4)) - 1/4*y^2 < 0
    */

    // fast test if at all
    // main cardioid is fully contained therein
    /*
		quick check: rectangular region encompassing the main cardioid

		solve -3/256+x/8-3*x^2/8-3*y^2/8+x^4+2*x^2*y^2+y^4 < 0
		for x as input variable and as well for y as input variable
		to get (e.g. using WolframAlpha):

		main cardioid lies in rectangle [-3/4..3/8] x [-11/16..11/16]*i
    */

	if (A.re.right < (double)(-0.75)) return INP12_OUTSIDE;
	if (A.re.left  > (double)( 0.375)) return INP12_OUTSIDE;
	if (A.im.right < (double)(-0.6875)) return INP12_OUTSIDE;
	if (A.im.left  > (double)( 0.6875)) return INP12_OUTSIDE;

	// using interval arithmetics to evaluate
	// 256*x^4+512*x^2*y^2-96*x^2+32*x+256*y^4-96*y^2-3 < 0
	// A.re =: x, A.im =: y
	CustInterval x4,x2,y2,y4;
	CustInterval Are,Aim;
	Are.copyFrom(A.re);
	Aim.copyFrom(A.im);
	if (custMul_ZAB(x2,Are,A.re) != 0) return INP12_ERROR;
	if (custMul_ZAB(y2,Aim,A.im) != 0) return INP12_ERROR;
	CustInterval t1;
	t1.copyFrom(x2);
	if (custMul_ZAB(x4,x2,t1) != 0) return INP12_ERROR;
	t1.copyFrom(y2);
	if (custMul_ZAB(y4,y2,t1) != 0) return INP12_ERROR;

	CustInterval v256,v512,v96,v3,v32;
	v3  .left=v3  .right=3.0;
	v32 .left=v32 .right=32.0;
	v96 .left=v96 .right=96.0;
	v256.left=v256.right=256.0;
	v512.left=v512.right=512.0;

	// 256*x^4+512*x^2*y^2-96*x^2+32*x+256*y^4-96*y^2-3 < 0
	CustInterval w1;
	if (custMul_ZAB(w1,v256,x4) != 0) return INP12_ERROR;
	// w1+512*x^2*y^2-96*x^2+32*x+256*y^4-96*y^2-3 < 0
	CustInterval w2;
	if (custMul_ZAB(w2,v512,x2) != 0) return INP12_ERROR;
	// w1+w2*y^2-96*x^2+32*x+256*y^4-96*y^2-3 < 0
	CustInterval w3;
	if (custMul_ZAB(w3,w2,y2) != 0) return INP12_ERROR;
	// w1+w3-96*x^2+32*x+256*y^4-96*y^2-3 < 0
	CustInterval w4;
	if (custAdd_ZAB(w4,w1,w3) != 0) return INP12_ERROR;
	// w4-96*x^2+32*x+256*y^4-96*y^2-3 < 0
	CustInterval w5;
	if (custMul_ZAB(w5,v96,x2) != 0) return INP12_ERROR;
	// w4-w5+32*x+256*y^4-96*y^2-3 < 0
	CustInterval w6;
	if (custSub_ZAB(w6,w4,w5) != 0) return INP12_ERROR;
	// w6+32*x+256*y^4-96*y^2-3 < 0
	CustInterval w7;
	if (custMul_ZAB(w7,v32,A.re) != 0) return INP12_ERROR;
	// w6+w7+256*y^4-96*y^2-3 < 0
	CustInterval w8;
	if (custAdd_ZAB(w8,w6,w7) != 0) return INP12_ERROR;
	// w8+256*y^4-96*y^2-3 < 0
	CustInterval w9;
	if (custMul_ZAB(w9,v256,y4) != 0) return INP12_ERROR;
	// w8+w9-96*y^2-3 < 0
	CustInterval w10;
	if (custAdd_ZAB(w10,w8,w9) != 0) return INP12_ERROR;
	// w10-96*y^2-3 < 0
	CustInterval w11;
	if (custMul_ZAB(w11,v96,y2) != 0) return INP12_ERROR;
	// w10-w11-3 < 0
	CustInterval w12;
	if (custSub_ZAB(w12,w10,w11) != 0) return INP12_ERROR;
	// w12-3 < 0
	CustInterval test;
	if (custSub_ZAB(test,w12,v3) != 0) return INP12_ERROR;

	// < 0 => inside
	if (test.right < 0.0) return INP12_INSIDE;

	// fully outside
	if (test.left > 0.0) return INP12_OUTSIDE;

	// straddling => can be either
	return INP12_MAYBE;
}

int8_t predP2_z2c(CustComplex& A) {

    // return: constant of INP12_XXX

	// circle of radius 0.25 around -1+0i
	// first check rectangle for faster result

	if (A.re.right < (double)(-1.25)) return INP12_OUTSIDE;
	if (A.re.left  > (double)(-0.75)) return INP12_OUTSIDE;
	if (A.im.right < (double)(-0.25)) return INP12_OUTSIDE;
	if (A.im.left  > (double)( 0.25)) return INP12_OUTSIDE;

    // IA-distance of A to poiint (-1,0)
    CustComplex p2hc;
    SETCPLX(p2hc,-1,0)
    CustComplex dist;
    if (custCplxSub_ZAB(dist,A,p2hc) != 0) return INP12_ERROR;

    CustInterval normQ;
    if (custCplxNormQ_ZA(normQ,dist) != 0) return INP12_ERROR;

    // if squared distance's right interval end < (0.25)^2 => fully inside
    if (normQ.right < (double)(0.25*0.25)) return INP12_INSIDE;

    // if squared distance's left interval end > (0.25)^2 => fully outside
    if (normQ.left > (double)(0.25*0.25)) return INP12_OUTSIDE;

    // else straddling => could be inside, outside (outward rounding result)
	return INP12_MAYBE;
}

void setCustConstants(void) {
    SETCPLX(cplx4,4,0)
    SETCPLX(cplxm4,-4,0)
    SETCPLX(cplx7,7,0)
    SETCPLX(cplx8,8,0)

}


// main

int32_t main(int32_t argc,char** argv) {
    if (argc < 2) {
        helppage();
        return 0;
    }

	// /////////////////////

	// initialize number type prerequisites
    CustRoundingUpwards *rd=new CustRoundingUpwards;
    if (!rd) {
        printf("\nerror. cannot set upward rounding mode\n");
        exit(99);
    }

    initializeCustInt(); // for Taylor polynomials
    globalafpath.init();

    // ////////////////////////

	setCustConstants();

	time_at_enter_main=clock();

	flog=fopen("upperbound-rd.log.txt","at");
	fprintf(flog,"\n-----------------\n");
	printf("upperbound-rd\n");
	// standard values
	// Mset and Jsets are computed in the 2-square
	RANGEM0=-2;
	RANGEM1=2;
	COMPLETEM0=RANGEM0;
	COMPLETEM1=RANGEM1;
	// standard values in case non given
	REFINEMENTM=12;
	SCREENWIDTHM=(1 << REFINEMENTM);

	for(int32_t i=1;i<argc;i++) {
		upper(argv[i]);

		if (!strcmp(argv[i],"ANEW")) {
            CREATESTART=1;
		} else
		if (!strcmp(argv[i],"P123")) {
            P123ON=1;
		} else
		if (strstr(argv[i],"LOAD=RAW")) {
            LOADFORMAT=FORMAT_RAW;
		} else
		if (strstr(argv[i],"SAVE=RAW")) {
            SAVEFORMAT=FORMAT_RAW;
		} else
		if (strstr(argv[i],"MAXIT=")) {
			int32_t a;
			if (sscanf(&argv[i][6],"%i",&a) == 1) {
				MAXIT=a;
			}
		} else
		if (strstr(argv[i],"SKIPLESS=")) {
			int32_t a;
			if (sscanf(&argv[i][9],"%i",&a) == 1) {
                SKIPIFLESS=a;
			}
		} else
		if (strstr(argv[i],"LENM=")==argv[i]) {
			// refinement level for the Mandelbrot set, i.e. the active region
			int a;
			if (sscanf(&argv[i][5],"%i",&a) == 1) {
				if (a < 8) a=8;
				REFINEMENTM=a;
				SCREENWIDTHM=(1 << a);
			}
		} else
		if (strstr(argv[i],"LENT=")==argv[i]) {
			int a;
			if (sscanf(&argv[i][5],"%i",&a) == 1) {
				if (a < 8) a=8;
				TILELEVEL=a;
			}
		}
	} // i

	if (REFINEMENTM > 30) {
        LOGMSG("\nerror. only implemented up to M30\n");
        exit(99);
	}

	if (CREATESTART > 0) {
        REFINEMENTM=TILELEVEL+1;
	}

    if (
        (TILELEVEL > MAXTILELEVEL) ||
        (REFINEMENTM > MAXREFINEMENTM) ||
        (TILELEVEL >= REFINEMENTM)
    ) {
        printf("\nerror. tilelevel must be <= %i, lenm <= %i\n",
            MAXTILELEVEL,MAXREFINEMENTM);
        exit(99);
    }

    TILELEN=(int64_t)1 << TILELEVEL;

    if (MAXIT < 0) {
        MAXIT = 3*REFINEMENTM;
    }

    if (AFFINEFOLDBACK < 0) {
        AFFINEFOLDBACK = 8*REFINEMENTM;
    }

    if (AFFINEFOLDBACK > MAXSYMBOLSPERAA) {
        LOGMSG3("\nerror. affine forms do not provide enough noise symbols (%I64d <-> %I64d).\n",
            (int64_t)MAXSYMBOLSPERAA,(int64_t)AFFINEFOLDBACK );
        exit(99);
    }

    if (TILELEVEL > 30) {
        printf("\nerror. cannot compute pixel count\N");
        exit(99);
    }

    tilepixelsint=TILELEN;
    tilepixelsint *= TILELEN;

	DICTLEVEL=REFINEMENTM - TILELEVEL;
	DICTLEN=(int64_t)1 << DICTLEVEL;
    loadDictionary();

    if (USEDICTIONARY <= 0) {
        if (CREATESTART <= 0) {
            LOGMSG("\nerror. no dictionary found.\n");
            exit(99);
        }
    }

	// and save directly, so one can prematurely
	// terminate the calculation wiuthout losing
	// data
	saveDictionaryAndTwice();

	SCREENWIDTHM=(int64_t)1 << REFINEMENTM;
	GRIDWIDTH=(TILELEN >> 4) + 3; // exact value is not relevant

	if (REFINEMENTM <= TILELEVEL) {
        LOGMSG2("\nerror. refinement must be larger than TileLevel %i\n",TILELEVEL);
        exit(99);
	}

	LOGMSG("\nSTATUS:\n");
	LOGMSG("  for levels up to M30: using CustInterval\n");
	if (USEAFFINE > 0) {
        LOGMSG("    using affine arithmetics if cust failed for orbit construction\n");
        LOGMSG2("    at %I64d noise symbols, af forms are folded back\n",AFFINEFOLDBACK);

	}

    if (P123ON > 0) {
        LOGMSG("  check for p1,p2 membership to exclude from further refinement\n");
    }

    if (USEGRIDFILL > 0) {
        LOGMSG3("  using gridfill at %i loops and jump size %I64d\n",GRIDLOOPS,GRIDFILLJUMP);
    }

    if (SKIPIFLESS > 0) {
        LOGMSG2("  if fewer than %I64d to check => file terminated",SKIPIFLESS);
    }

    if (USEAFFINE) {
        if (GRIDLOOPS < 16) GRIDLOOPS=16;
    }

	setFunc();
	LOGMSG4("  function=%s sym xaxis=%i sym yaxis=%i\n",
		funcname[FUNC],SYMMETRY_TO_XAXIS,SYMMETRY_TO_YAXIS);
    LOGMSG2("  max it %i\n",MAXIT);

	RANGEM0=-RANGEM1;
	// set AFTER setFunc call
	COMPLETEM0=-RANGEM1;
	COMPLETEM1=RANGEM1;

	char tt[2048];
	sprintf(tt,"%i",-RANGEM1);
	COMPLETEM0cust.setStr(tt);
	sprintf(tt,"%i",RANGEM1);
	COMPLETEM1cust.setStr(tt);

	CustInterval acust,bcust;
	sprintf(tt,"%i",(RANGEM1-RANGEM0));
	acust.setStr(tt);
	sprintf(tt,"%I64d",SCREENWIDTHM);
	bcust.setStr(tt);

	if (custDiv_ZAB(scaleRangePerPixelMcust,acust,bcust) != 0) {
        LOGMSG("\nerror. cannot compute scale.\n");
        exit(99);
	}

	if (scaleRangePerPixelMcust.left != scaleRangePerPixelMcust.right) {
        printf("\nerror. cust scale must be a degenrate interval.\n");
        exit(99);
	}

	LOGMSG2("  Mset %i-square\n",RANGEM1);
	LOGMSG("  scale ");
	scaleRangePerPixelMcust.ausgabe(stdout);
	scaleRangePerPixelMcust.ausgabe(flog);

	// calculate scale values
	// valid up to M48 for -.2..+2 as range
	double w=(RANGEM1-RANGEM0) / (double)SCREENWIDTHM;
	scaleRangePerPixelM=w; // dyadic fraction

	LOGMSG3("\n  M%02i T%02i\n",REFINEMENTM,TILELEVEL);

	// needed for twice save
	md2.setlenxy(TILELEN,TILELEN);
	setPaletteTo(md2);

	printf("\nanalysis ...");

	// //////////////////////////////////////
	//
	// main routine
	//

	if (CREATESTART > 0) {
        printf("\ncreating initial data ... ");
        createanew();
	} else {
        // main routine
        analyze();

        if (SYMMETRY_TO_XAXIS > 0) {
            if (ctrexterior < INT61) {
                ctrexterior *= 2;
            } else {
                LOGMSG("\nerror. int64_t not sufficient.\n");
                exit(99);
            }
        }

        if (SYMMETRY_TO_YAXIS > 0) {
            if (ctrexterior < INT61) {
                ctrexterior *= 2;
            } else {
                LOGMSG("\nerror. int64_t not sufficient.\n");
                exit(99);
            }
        }

    }

	saveDictionaryAndTwice();

	//
	/////////////////////////////

	if (CREATESTART <= 0) {
        LOGMSG("\n\n----------------\nRESULT:\n");
        LOGMSG("  (reliable as pixels represent complex intervals and the fate of\n");
        LOGMSG("   ALL underlying complex numbers):\n");
        LOGMSG3("\n  M(z^2+c) <= 16 - %I64d * 4 / (2^%i)^2\n",
            ctrexterior,REFINEMENTM);
    }

	int32_t c1=clock();\
	LOGMSG2("\nduration %.0lf sec\n",(double)(c1-time_at_enter_main)/CLOCKS_PER_SEC);\

	fclose(flog);

    return 0;
}

