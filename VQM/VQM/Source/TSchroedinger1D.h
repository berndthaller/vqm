// ===========================================================================
//	TSchroedinger2D.h
//	Project: Visual Quantum Mechanics.
//  Supervision: Bernd Thaller.
//  Programming: Wolfgang Thaller 1998
// ===========================================================================

#ifndef _H_TSchroedinger1D
#define _H_TSchroedinger1D
#pragma once

#include "TypeDefinition.h"
#include "TOperator.h"


class	TSchroedinger1D : public TOperator
{
	public:
				TSchroedinger1D(	Int32 inScalarID, 
									Float inMass,
									Float inUnits );
				~TSchroedinger1D( void );
		Int32	TimeEvolution(	TFunction* inFunction,
								Float inTimeStep,
								Int32 inFractal,
								Int32 inSteps );
		Int32	PutInfo( void );
private:
		Int32	mScalarID;
		Float	mB,mUnits;
};

#endif
