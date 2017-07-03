// ===========================================================================
//	TDomain.cp
//	Project: Visual Quantum Mechanics.
//  Supervision: Bernd Thaller.
//  Programming: Manfred Liebmann 1998
// ===========================================================================
//
//	Class for operator-domains

#include "MathLinkUtilities.h"
#include "TDomain.h"



// ---------------------------------------------------------------------------
//		$ TDomain
// ---------------------------------------------------------------------------
//	Constructor

TDomain::TDomain( void )
{
	mDomainP = NULL;
}



// ---------------------------------------------------------------------------
//		$ ~TDomain
// ---------------------------------------------------------------------------
//	Destructor

TDomain::~TDomain( void )
{	
	delete [] mDomainP;
}



// ---------------------------------------------------------------------------
//		$ InitPlain2D
// ---------------------------------------------------------------------------
//	Create 2D domain array, plain format

Int32	TDomain::InitPlain2D(	Int8* &outDomP,
								Float* inD0P,
								Int32 inNi,
								Int32 inNj )
{
	Float	zero = 0.0;
	Int32	i, j;
	Int8	*domP, type;


	domP = new Int8[inNi*inNj];
	if( !domP )
		return MLErrorReport(stdlink, "out of memory");
	outDomP = mDomainP = domP;

	for( j = 0; j < inNj; j++ ) {
		for( i = 0; i < inNi; i++ ) {
			type = 0xff;
			if( !inD0P || *inD0P++ > zero ) {
				if( i == inNi-1 || i == 0 || j == inNj-1 || j == 0 ) type = ~type;
			}
			else type = ~type;
			*domP++ = type;
		}
	}

	return eOK;
}



// ---------------------------------------------------------------------------
//		$ InitPlain3D
// ---------------------------------------------------------------------------
//	Create 3D domain array, plain format

Int32	TDomain::InitPlain3D(	Int8* &outDomP,
								Float* inD0P,
								Int32 inNi,
								Int32 inNj,
								Int32 inNk )
{
	Float	zero = 0.0;
	Int32	i, j, k;
	Int8	*domP, type;


	domP = new Int8[inNi*inNj*inNk];
	if( !domP )
		return MLErrorReport(stdlink, "out of memory");
	outDomP = mDomainP = domP;

	for( k = 0; k < inNk; k++ ) {
		for( j = 0; j < inNj; j++ ) {
			for( i = 0; i < inNi; i++ ) {
				type = 0xff;
				if( !inD0P || *inD0P++ > zero ) {
					if( i == inNi-1 || i == 0 ||
						j == inNj-1 || j == 0 ||
						k == inNk-1 || k == 0 )
						type = ~type;
				}
				else type = ~type;
				*domP++ = type;
			}
		}
	}

	return eOK;
}



// ---------------------------------------------------------------------------
//		$ ClearC
// ---------------------------------------------------------------------------
//	Clear boundary points

void	TDomain::ClearC(	Float* inS1P,
							Float* inS2P,
							Int32 inLen )
{
	Float	zero = 0.0;
	Int32	i;
	Int8	*domP;

	domP = mDomainP;
	for( i = 0; i < inLen; i++ ) {
		if( !domP[i] ) inS1P[i] = inS2P[i] = zero;
	}

}



// ---------------------------------------------------------------------------
//		$ ClearC2
// ---------------------------------------------------------------------------
//	Clear boundary points

void	TDomain::ClearC2(	Float* inS1P,
							Float* inS2P,
							Float* inS3P,
							Float* inS4P,
							Int32 inLen )
{
	Float	zero = 0.0;
	Int32	i;
	Int8	*domP;

	domP = mDomainP;
	for( i = 0; i < inLen; i++ ) {
		if( !domP[i] ) inS1P[i] = inS2P[i] = inS3P[i] = inS4P[i] = zero;
	}


}

// ---------------------------------------------------------------------------
//		$ ClearC4
// ---------------------------------------------------------------------------
//	Clear boundary points

void	TDomain::ClearC4(	Float* inS1P,
							Float* inS2P,
							Float* inS3P,
							Float* inS4P,
							Float* inS5P,
							Float* inS6P,
							Float* inS7P,
							Float* inS8P,
							Int32 inLen )
{
	Float	zero = 0.0;
	Int32	i;
	Int8	*domP;

	domP = mDomainP;
	for( i = 0; i < inLen; i++ ) {
		if( !domP[i] )
			inS1P[i] = inS2P[i] = inS3P[i] = inS4P[i] =
			inS5P[i] = inS6P[i] = inS7P[i] = inS8P[i] = zero;
	}

}

