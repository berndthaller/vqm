// ===========================================================================
//	TOperator.cp
//	Project: Visual Quantum Mechanics.
//  Supervision: Bernd Thaller.
//  Programming: Manfred Liebmann 1998
// ===========================================================================
//
//	Abstract class for matrix-operators

#include "MathLinkUtilities.h"
#include "TOperator.h"
#include <math.h>



// ---------------------------------------------------------------------------
//		$ GetFractalNumber
// ---------------------------------------------------------------------------
//	Calculate fractal numbers

Int32	TOperator::GetFractalNumber(	Int32 inFractal,
										Int32 inIndex,
										Float &outReal,
										Float &outImag )
{
	Int32	i;
	Float	reP, imP, reT, imT, twom, temp;
	Float	zero = 0.0, half = 0.5, one = 1.0, two = 2.0;


	if( inFractal >= 32 ) return eError;
	
	reP = one;
	imP = zero;

	twom = two;
	for( i = 1; i < inFractal; i++ ) {
		twom += two;
		temp = tan( M_PI / twom );

		if( inIndex & 1 ) temp = -temp;
		inIndex >>= 1;
		reT = reP - imP * temp;
		imT = imP + reP * temp;
		reP = reT * half;
		imP = imT * half;
	}

	outReal = reP;
	outImag = imP;

	return eOK;
}

