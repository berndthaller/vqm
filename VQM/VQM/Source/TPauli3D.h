// ===========================================================================
//	TPauli3D.h
//	Project: Visual Quantum Mechanics.
//  Supervision: Bernd Thaller.
//  Programming: Manfred Liebmann 1998
// ===========================================================================

#ifndef _H_TPauli3D
#define _H_TPauli3D
#pragma once

#include "TypeDefinition.h"
#include "TOperator.h"


class	TPauli3D : public TOperator
{
	public:
				TPauli3D(	Int32 inScalarID,
							Int32 inVectorID,
							Int32 inDomainID,
							Float inMass,
							Float inCharge,
							Float inUnits );
				~TPauli3D( void );
		Int32	TimeEvolution(	TFunction* inFunction,
								Float inTimeStep,
								Int32 inFractal,
								Int32 inSteps );
		Int32	PutInfo( void );

	private:
		void	Kernel(	Float* rePsi1P,
						Float* imPsi1P,
						Float* rePsi2P,
						Float* imPsi2P,
						Float* rePhi1P,
						Float* imPhi1P,
						Float* rePhi2P,
						Float* imPhi2P,
						Float* v0P,
						Float* w0P,
						Float* a1P,
						Float* a2P,
						Float* a3P,
						Int8* domP,
						Float reZ,
						Float imZ,
						Int32 ni,
						Int32 nj,
						Int32 nk );

		Int32	mScalarID;
		Int32	mVectorID;
		Int32	mDomainID;
		Float	mCharge;
		Float	mMass;
		Float	mUnits;

};

#endif
