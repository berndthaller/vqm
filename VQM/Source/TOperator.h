// ===========================================================================
//	TOperator.h
//	Project: Visual Quantum Mechanics.
//  Supervision: Bernd Thaller.
//  Programming: Manfred Liebmann 1998
// ===========================================================================
//
// Wolfgang Thaller: Made TimeEvolution and PutInfo PURE VIRTUAL funcitons
//

#ifndef _H_TOperator
#define _H_TOperator
#pragma once

#include "TypeDefinition.h"
#include "TFunction.h"
#include "TList.h"


class	TOperator
{
	public:
		virtual Int32	TimeEvolution(	TFunction* inFunction,
										Float inTimeStep,
										Int32 inFractal,
										Int32 inSteps ) = 0;
		virtual Int32	PutInfo( void ) = 0;

		Int32	mID;

		virtual			~TOperator() {}
	protected:
		Int32	GetFractalNumber(	Int32 inFractal,
									Int32 inIndex,
									Float &outReal,
									Float &outImag );

};

extern TList<TOperator> *gOperatorList;

enum
{
	kScalar = 1 << 0,
	kVector = 1 << 1
};

#endif
