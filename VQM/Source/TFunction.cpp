// ===========================================================================
//	TFunction.cp
//	Project: Visual Quantum Mechanics.
//  Supervision: Bernd Thaller.
//  Programming: Manfred Liebmann 1998
//
// 	Changes by Wolfgang Thaller on 15.08.98:
//
//	-		Added function: IsFunction1D
//
// ===========================================================================
// 
//	Class for functions

//RM2018
#include <iostream>
using namespace std;

#include "MathLinkUtilities.h"
#include "TFunction.h"
#include "TOperator.h"

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>


// ---------------------------------------------------------------------------
//		$ TFunction
// ---------------------------------------------------------------------------
//	Constructor

TFunction::TFunction( void )
{
	mArrayP = NULL;
	mCountP = NULL;
	mHeadH = NULL;
	mDepth = 0;
	mID = 0;
}



// ---------------------------------------------------------------------------
//		$ ~TFunction
// ---------------------------------------------------------------------------
//	Destructor

TFunction::~TFunction( void )
{
	delete [] mArrayP;
	delete [] mCountP;
}



// ---------------------------------------------------------------------------
//		$ GetArray
// ---------------------------------------------------------------------------
//	Read float array from MathLink

//JPK:
//int n;
//complex c;
//If(MLGetFunction(stdlink,"Complex",&n) && 2==n){
//double re,im;
//MLGetDouble(stdlink,&re);
//MLGetDouble(stdlink,&im);
//c=complex(re,im);
//}


Int32	TFunction::GetArray( void )
{
//	MLPrint(stdlink, "................... in GetArray");
//	MLPrintReal(stdlink, (float) mDepth);
	return MLGetRealArray2(stdlink, mArrayP, mCountP, mDepth);
}



// ---------------------------------------------------------------------------
//		$ PutArray
// ---------------------------------------------------------------------------
//	Return float array to MathLink

Int32	TFunction::PutArray( void )
{	
 //RM2018	return MLPutFloatArray(stdlink, mArrayP, mCountP, mHeadH, mDepth);
//	return MLPutReal32Array(stdlink, mArrayP, (int *)mCountP, (const char **)mHeadH, mDepth);
// RM2018: Only in MacOSX this does not work :
//	return MLPutReal32Array(stdlink, (float *)mArrayP, (int *)mCountP, (const char **)mHeadH, mDepth);

	    Int32 i;
		Int32 size = 1;
		for(i=0;i<mDepth;i++)
		{
			size *= mCountP[i];
		}
		return MLPutReal32List(stdlink, (float *)mArrayP, size);
}

// RM: commmented out 2007-06-28
// ---------------------------------------------------------------------------
//		$ PutToSocket
// ---------------------------------------------------------------------------
//	Send float array over socket
//
//void	TFunction::PutToSocket( int sock )
//{
//	Int32 tmp;
//	int i;
//	Int32 size = 1;
//	
//	tmp = htonl(mDepth);
//	write2(sock,&tmp,4);
//	
//	for(i=0;i<mDepth;i++)
//	{
//		size *= mCountP[i];
//		tmp = htonl(mCountP[i]);
//		write2(sock,&tmp,4);
//	}
//	
//	if(htonl(0xDEADBEEF) == 0xDEADBEEF)
//		write2(sock,mArrayP,4*size);//###
//	else
//		for(i=0;i<size;i++)
//		{
//			float tmpF = mArrayP[i];
//			tmp = htonl(*(Int32*) &tmpF);
//			write2(sock,&tmp,4);
//		}
//}


// ---------------------------------------------------------------------------
//		$ PutInfo
// ---------------------------------------------------------------------------
//	Return list of array dimensions

Int32	TFunction::PutInfo( void )
{
	Int32	i;
	
	MLPutFunction(stdlink, "List", mDepth);
	for( i = 0; i < mDepth; i++ )
//		MLPutLongInteger(stdlink, mCountP[mDepth-1-i]);
//RM2018
		MLPutInteger32(stdlink, mCountP[mDepth-1-i]);
	
	return eOK;
}



// ---------------------------------------------------------------------------
//		$ PutColor
// ---------------------------------------------------------------------------
//	Return RGB color array

Int32	TFunction::PutColor( void )
{
	Int32	i, size;	
	Float	*realP, *imagP, *colorP;
//RM20180422 add const:
	Int8	const *headH[3] = {"List", "List", "RGBColor"};
	Int32	countP[3] = {0, 0, 3};
	Float	re, im, r, s, s0, phi, phi0, red, green, blue, temp;
	Float	zero = 0.0, one = 1.0, two = 2.0;
	Float	k3pi = 3.0 / M_PI, k4pi = 4.0 / M_PI;


	if( !IsFunction2D(countP[1], countP[0]) )
		return MLErrorReport(stdlink, "two-dimensional function expected");
	if( !IsFunctionC(realP, imagP) )
		return MLErrorReport(stdlink, "complex function expected");

	size = countP[0] * countP[1];
	colorP = new Float[3 * size];
	if( !colorP )
		return MLErrorReport(stdlink, "out of memory");
	if( eError == MLCheckMemoryReserve(stdlink) ) return eError;

	for( i = 0; i < size; i++ ) {
		re = *realP++;
		im = *imagP++;

		r = hypot( re, im );
		s = k4pi * atan( r );
		phi = k3pi * atan2( im, re );
		phi0 = fabs( phi );

		if( phi0 < one ) { red = one; green = phi0; blue = zero; }
		else {
			if( phi0 < two ) { red = two - phi0; green = one; blue = zero; }
			else { red = zero; green = one; blue = phi0 - two; }
		}

		if( phi < zero ) { temp = green; green = blue; blue = temp; }
		
		if( s < one ) { red *= s; green *= s; blue *= s; }
		else {
			s = two - s; s0 = one - s;
			red = red * s + s0; green = green * s + s0; blue = blue * s + s0;
		}

		if( !isfinite(red) ) red = zero;
		if( !isfinite(green) ) green = zero;
		if( !isfinite(blue) ) blue = zero;		
		*colorP++ = red;
		*colorP++ = green;
		*colorP++ = blue;
	}
	colorP -= 3*size;


//RM2018	MLPutFloatArray(stdlink, colorP, countP, headH, 3);
//	MLPutReal32Array(stdlink, colorP, (int *)countP, (const char **)headH, 3);
MLPutReal32List(stdlink, (float *)colorP, 3*size);
	delete [] colorP;

	return eOK;
}



// ---------------------------------------------------------------------------
//		$ PutGray
// ---------------------------------------------------------------------------
//	Return gray array

Int32	TFunction::PutGray( void )
{
	Int32	i, size;
	Float	*realP, *imagP, *grayP;
//RM20180420 add const:
	Int8	const *headH[3] = {"List", "List", "GrayLevel"};
	Int32	countP[3]={0, 0, 1};
	Float	re, im, r, gray;
	Float	zero = 0.0, one = 1.0;
	Float	k2pi = 2.0 / M_PI;


	if( !IsFunction2D(countP[1], countP[0]) )
		return MLErrorReport(stdlink, "two-dimensional function expected");
	if( !IsFunctionC(realP, imagP) )
		return MLErrorReport(stdlink, "complex function expected");

	size = countP[0] * countP[1];
	grayP = new Float[size];
	if( !grayP )
		return MLErrorReport(stdlink, "out of memory");
	if( eError == MLCheckMemoryReserve(stdlink) ) return eError;

	for( i = 0; i < size; i++ ) {
		re = *realP++;
		im = *imagP++;
		r = re * re + im * im;
		gray = k2pi * atan( r );

		if( !isfinite(gray) ) gray = zero;
		*grayP++ = gray;
	}
grayP -= size;

//RM2018	MLPutFloatArray(stdlink, grayP, countP, headH, 3);
//	MLPutReal32Array(stdlink, grayP, (int *)countP, (const char **)headH, 3);
	MLPutReal32List(stdlink, (float *)grayP, size);

	delete [] grayP;
	return eOK;
}



// ---------------------------------------------------------------------------
//		$ PutRedBlue
// ---------------------------------------------------------------------------
//	Return red-blue color array

Int32	TFunction::PutRedBlue( void )
{
	Int32	i, size;	
	Float	*realP, *colorP;
//RM20180422 add const:
	Int8	const *headH[3] = {"List", "List", "RGBColor"};
	Int32	countP[3] = {0, 0, 3};
	Float	r, s, red, green, blue, temp;
	Float	zero = 0.0, one = 1.0;
	Float	k4pi = 4.0 / M_PI;


	if( !IsFunction2D(countP[1], countP[0]) )
	return MLErrorReport(stdlink, "two-dimensional function expected");
	if( !IsFunctionR(realP) )
		return MLErrorReport(stdlink, "real function expected");

	size = countP[0] * countP[1];
	colorP = new Float[3 * size];
	if( !colorP )
		return MLErrorReport(stdlink, "out of memory");
	if( eError == MLCheckMemoryReserve(stdlink) ) return eError;

	for( i = 0; i < size; i++ ) {
		r = *realP++;
		s = k4pi * atan( fabs( r ) );

		if( s < one ) { red = s; green = blue = zero; }
		else { red = one; green = blue = s - one; }
		
		if( r < zero ) { temp = red; red = blue; blue = temp; }
		
		if( !isfinite(red) ) red = zero;
		if( !isfinite(green) ) green = zero;
		if( !isfinite(blue) ) blue = zero;		
		*colorP++ = red;
		*colorP++ = green;
		*colorP++ = blue;
	}
	colorP -= 3*size;

//RM2018	MLPutFloatArray(stdlink, colorP, countP, headH, 3);
	MLPutReal32Array(stdlink, colorP, (int *)countP, (const char **)headH, 3);
	delete [] colorP;

	return eOK;
}



// ---------------------------------------------------------------------------
//		$ PutBlackWhite
// ---------------------------------------------------------------------------
//	Return black-white array


Int32	TFunction::PutBlackWhite( void )
{
	Int32	i, size;
 	Float	*realP, *grayP;
//RM20180422 add const
 	Int8	const *headH[3] = {"List", "List", "GrayLevel"};
 	Int32	countP[3] = {0, 0, 1};
 	Float	gray;
 	Float	zero = 0.0, one = 1.0;
 
	if( !IsFunction2D(countP[1], countP[0]) )
		return MLErrorReport(stdlink, "two-dimensional function expected");
	if( !IsFunctionR(realP) )
 		return MLErrorReport(stdlink, "real function expected");
 
 	size = countP[0] * countP[1];
 	grayP = new Float[size];
 	if( !grayP )
 		return MLErrorReport(stdlink, "out of memory");
 	if( eError == MLCheckMemoryReserve(stdlink) ) return eError;
 
 	for( i = 0; i < size; i++ ) {
 		if( *realP++ < zero ) gray = one;
 		else gray = zero;
 		
 		*grayP++ = gray;
 	}
 	grayP -= size;
 
 //RM2018	MLPutFloatArray(stdlink, grayP, countP, headH, 3);
 	MLPutReal32Array(stdlink, grayP, (int *)countP, (const char **)headH, 3);
 	delete [] grayP;
 
	return eOK;
}



// ---------------------------------------------------------------------------
//		$ PutAbs
// ---------------------------------------------------------------------------
//	Return abs array
//RM2018 https://stackoverflow.com/a/33630813/887505

Int32	TFunction::PutAbs( void )
{	
	    Int32 i;
		Int32 size = 1;
	    Int32 countP[2] = {0, 0};
 	    Float	*realP, *imagP, *absP, re, im, r;
	 	Float	zero = 0.0;

//RM2018:
		if( !IsFunction2D(countP[1], countP[0]) )
		return MLErrorReport(stdlink, "two-dimensional function expected");

 	    if( !IsFunctionC(realP, imagP) )
		return MLErrorReport(stdlink, "complex function expected");

		size = countP[0] * countP[1];

 	absP = new Float[size];

 	if( !absP )
 		return MLErrorReport(stdlink, "out of memory");
 	if( eError == MLCheckMemoryReserve(stdlink) ) return eError;
 	
 	for( i = 0; i < size; i++ ) {
		re = *realP++;
		im = *imagP++;
 		r = hypot( re, im );
 
 		if(!isfinite(r) ) r = zero;
 		*absP++ = r;
 	}

//ASK: why ?
 	absP -= size;


   MLPutReal32List(stdlink, (float *)absP, size);
 	delete [] absP;
	return eOK;
}


// ---------------------------------------------------------------------------
//		$ IsFunction1D
// ---------------------------------------------------------------------------
//	Check if function is 1D

bool	TFunction::IsFunction1D(	Int32 &ioNi )
{
	if( mDepth != 2 ) return false;
	
	if( ioNi == 0) {
		ioNi = mCountP[1];
	}
	else if(ioNi != mCountP[1] ) return false;

	return true;
}

// ---------------------------------------------------------------------------
//		$ IsFunction2D
// ---------------------------------------------------------------------------
//	Check if function is 2D

bool	TFunction::IsFunction2D(	Int32 &ioNi,
									Int32 &ioNj )
{
//	MLPrint(stdlink, "mDepth = ");
//	MLPrintReal(stdlink, mDepth);
	if( mDepth != 3 ) return false;
	
	if( ioNi == 0 && ioNj == 0 ) {
		ioNj = mCountP[1];
		ioNi = mCountP[2];
	}
	else if( ioNj != mCountP[1] || ioNi != mCountP[2] ) return false;

//	MLPrint(stdlink, "in IsFunction2D ");
//	MLPrintReal(stdlink, mDepth);

	return true;
}



// ---------------------------------------------------------------------------
//		$ IsFunction3D
// ---------------------------------------------------------------------------
//	Check if function is 3D

bool	TFunction::IsFunction3D(	Int32 &ioNi,
									Int32 &ioNj,
									Int32 &ioNk )
{
	if( mDepth != 4 ) return false;
	
	if( ioNi == 0 && ioNj == 0 && ioNk == 0 ) {
		ioNk = mCountP[1];
		ioNj = mCountP[2];
		ioNi = mCountP[3];
	}
	else if( ioNk != mCountP[1] || ioNj != mCountP[2] || ioNi != mCountP[3] ) return false;

	return true;
}



// ---------------------------------------------------------------------------
//		$ IsFunctionR
// ---------------------------------------------------------------------------
//	Check if function is R-valued

bool	TFunction::IsFunctionR( Float* &outSP )
{
	if( mCountP[0] != 1 ) return false;
	
	outSP = mArrayP;
	return true;
}



// ---------------------------------------------------------------------------
//		$ IsFunctionR2
// ---------------------------------------------------------------------------
//	Check if function is R2-valued

bool	TFunction::IsFunctionR2(	Float* &outS1P,
									Float* &outS2P )
{
	return IsFunctionC(outS1P, outS2P);
}



// ---------------------------------------------------------------------------
//		$ IsFunctionR3
// ---------------------------------------------------------------------------
//	Check if function is R3-valued

bool	TFunction::IsFunctionR3(	Float* &outS1P,
									Float* &outS2P,
									Float* &outS3P )
{
	Int32	i, offset = 1;


	if( mCountP[0] != 3 ) return false;

	for( i = 1; i < mDepth; i++ ) offset *= mCountP[i];
	outS1P = mArrayP;
	outS2P = mArrayP + offset;
	outS3P = mArrayP + 2*offset;
	return true;

}



// ---------------------------------------------------------------------------
//		$ IsFunctionR4
// ---------------------------------------------------------------------------
//	Check if function is R4-valued

bool	TFunction::IsFunctionR4(	Float* &outS1P,
									Float* &outS2P,
									Float* &outS3P,
									Float* &outS4P )
{
	return IsFunctionC2(outS1P, outS2P, outS3P, outS4P);
}



// ---------------------------------------------------------------------------
//		$ IsFunctionC
// ---------------------------------------------------------------------------
//	Check if function is C-valued

bool	TFunction::IsFunctionC(	Float* &outS1P,
								Float* &outS2P )
{
	Int32	i, offset = 1;


//RM2018 ASK MIST MIST
	if( mCountP[0] != 2 ) return false;

	for( i = 1; i < mDepth; i++ )
	{
		offset *= mCountP[i];
		outS1P = mArrayP;
		outS2P = mArrayP + offset;
	}

//MLPrint(stdlink, "in IsFunctionC, offset = ");
//MLPrintReal(stdlink, (float) offset);

	return true;
}



// ---------------------------------------------------------------------------
//		$ IsFunctionC2
// ---------------------------------------------------------------------------
//	Check if function is C2-valued

bool	TFunction::IsFunctionC2(	Float* &outS1P,
									Float* &outS2P,
									Float* &outS3P,
									Float* &outS4P )
{
	Int32	i, offset = 1;


	if( mCountP[0] != 4 ) return false;

	for( i = 1; i < mDepth; i++ ) offset *= mCountP[i];
	outS1P = mArrayP;
	outS2P = mArrayP + offset;
	outS3P = mArrayP + 2*offset;
	outS4P = mArrayP + 3*offset;
	return true;
}



// ---------------------------------------------------------------------------
//		$ IsFunctionC4
// ---------------------------------------------------------------------------
//	Check if function is C4-valued

bool	TFunction::IsFunctionC4(	Float* &outS1P,
									Float* &outS2P,
									Float* &outS3P,
									Float* &outS4P,
									Float* &outS5P,
									Float* &outS6P,
									Float* &outS7P,
									Float* &outS8P )
{
	Int32	i, offset = 1;


	if( mCountP[0] != 8 ) return false;

	for( i = 1; i < mDepth; i++ ) offset *= mCountP[i];
	outS1P = mArrayP;
	outS2P = mArrayP + offset;
	outS3P = mArrayP + 2*offset;
	outS4P = mArrayP + 3*offset;
	outS5P = mArrayP + 4*offset;
	outS6P = mArrayP + 5*offset;
	outS7P = mArrayP + 6*offset;
	outS8P = mArrayP + 7*offset;
	return true;
}



// ---------------------------------------------------------------------------
//		$ UpdateWindow
// ---------------------------------------------------------------------------
//	Update window

Int32	TFunction::UpdateWindow( void )
{
#ifdef QUANTUM_KERNEL_UI
	TWindow*	theWindow;
	Int32	size, ID;

	size = gWindowList->GetSize();
	for( ID = 0; ID < size; ID++ ) {
		theWindow = gWindowList->Fetch(ID);
		if( theWindow ) {
			if( theWindow->HasCorrectID(mID) ) {
				if( eError == theWindow->Draw(this) ) return eError;
			}
		}
	}
	
#endif
	return eOK;
}



// ---------------------------------------------------------------------------
//		$ Copy
// ---------------------------------------------------------------------------
//	Copy function object

Int32	TFunction::Copy( TFunction* inFunction )
{
	Int32	i, size;

		
	size = 1;
	for( i = 0; i < inFunction->mDepth; i++ ) size *= inFunction->mCountP[i];

	delete [] mArrayP;
	mArrayP = new Float[size];
	if( !mArrayP ) {
		return MLErrorReport(stdlink, "out of memory");
	}
	for( i = 0; i < size; i++ ) mArrayP[i] = inFunction->mArrayP[i];

	delete [] mCountP;
	mCountP = new Int32[inFunction->mDepth];
	if( !mCountP ) {
		return MLErrorReport(stdlink, "out of memory");
	}
	for( i = 0; i < inFunction->mDepth; i++ ) mCountP[i] = inFunction->mCountP[i];
	
	mDepth = inFunction->mDepth;

	return eOK;
}



// ---------------------------------------------------------------------------
//		$ SwapArrayPointers
// ---------------------------------------------------------------------------
//	Swap array pointers

void	TFunction::SwapArrayPointers( TFunction* inFunction )
{
	Float	*theArrayP;

	theArrayP = inFunction->mArrayP;
	inFunction->mArrayP = mArrayP;
	mArrayP = theArrayP;
}
