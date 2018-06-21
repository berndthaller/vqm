// ===========================================================================
//	TPauli2D.h
//	Project: Visual Quantum Mechanics.
//  Supervision: Bernd Thaller.
//  Programming: Manfred Liebmann 1998
// ===========================================================================

#ifndef _H_TPauli2D
#define _H_TPauli2D
#pragma once

#include "TypeDefinition.h"
#include "TOperator.h"


class	TPauli2D : public TOperator
{
	public:
				TPauli2D(	Int32 inScalarID,
							Int32 inVectorID,
							Int32 inDomainID,
							Float inMass,
							Float inCharge,
							Float inUnits );
				~TPauli2D( void );
		Int32	TimeEvolution(	TFunction* inFunction,
								Float inTimeStep,
								Int32 inFractal,
								Int32 inSteps );
		Int32	PutInfo( void );

	private:
		void	Kernel(	Float* rePsiP,
						Float* imPsiP,
						Float* rePhiP,
						Float* imPhiP,
						Float* v0P,
						Float* w0P,
						Float* a1P,
						Float* a2P,
						Int8* domP,
						Float re,
						Float im,
						Int32 ni,
						Int32 nj );

		Int32	mScalarID;
		Int32	mVectorID;
		Int32	mDomainID;
		Float	mCharge;
		Float	mMass;
		Float	mUnits;

};

#endif
