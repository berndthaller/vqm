// ===========================================================================
//	TDomain.h
//	Project: Visual Quantum Mechanics.
//  Supervision: Bernd Thaller.
//  Programming: Manfred Liebmann 1998
// ===========================================================================

#ifndef _H_TDomain
#define _H_TDomain
#pragma once

#include "TypeDefinition.h"


class	TDomain
{
	public:
				TDomain( void );
				~TDomain( void );
		Int32	InitPlain2D(	Int8* &outDomP,
								Float* inD0P,
								Int32 inNi,
								Int32 inNj );
		Int32	InitPlain3D(	Int8* &outDomP,
								Float* inD0P,
								Int32 inNi,
								Int32 inNj,
								Int32 inNk );
		void	ClearC(	Float* inS1P,
						Float* inS2P,
						Int32 inLen );
		void	ClearC2(	Float* inS1P,
							Float* inS2P,
							Float* inS3P,
							Float* inS4P,
							Int32 inLen );
		void	ClearC4(	Float* inS1P,
							Float* inS2P,
							Float* inS3P,
							Float* inS4P,
							Float* inS5P,
							Float* inS6P,
							Float* inS7P,
							Float* inS8P,
							Int32 inLen );

	private:

		Int8*	mDomainP;

};

#endif
