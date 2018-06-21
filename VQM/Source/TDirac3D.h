// ===========================================================================
//	TDirac3D.h
//	Project: Visual Quantum Mechanics.
//  Supervision: Bernd Thaller.
//  Programming: Manfred Liebmann 1998
// ===========================================================================

#ifndef _H_TDirac3D
#define _H_TDirac3D
#pragma once

#include "TypeDefinition.h"
#include "TOperator.h"


class	TDirac3D : public TOperator
{
	public:
				TDirac3D(	Int32 inScalarID,
							Int32 inVectorID,
							Int32 inDomainID,
							Float inMass,
							Float inCharge,
							Float inUnits );
				~TDirac3D( void );
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
						Float* rePsi3P,
						Float* imPsi3P,
						Float* rePsi4P,
						Float* imPsi4P,
						Float* rePhi1P,
						Float* imPhi1P,
						Float* rePhi2P,
						Float* imPhi2P,
						Float* rePhi3P,
						Float* imPhi3P,
						Float* rePhi4P,
						Float* imPhi4P,
						Float* v0P,
						Float* w0P,
						Float* a1P,
						Float* a2P,
						Float* a3P,
						Float* a4P,
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
