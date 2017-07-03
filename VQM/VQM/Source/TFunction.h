// ===========================================================================
//	TFunction.h
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

#ifndef _H_TFunction
#define _H_TFunction
#pragma once

#include "TypeDefinition.h"
#include "TList.h"


class	TFunction
{
	public:
				TFunction( void );
				~TFunction( void );
		Int32	GetArray( void );
		Int32	PutInfo( void );
		Int32	PutArray( void );
		void	PutToSocket( int sock );

		Int32	PutColor( void );
		Int32	PutGray( void );
		Int32	PutRedBlue( void );
		Int32	PutBlackWhite( void );
		Int32	PutAbs( void );
		
		bool	IsFunction1D(	Int32 &ioNi );
		bool	IsFunction2D(	Int32 &ioNi,
								Int32 &ioNj );
		bool	IsFunction3D(	Int32 &ioNi,
								Int32 &ioNj,
								Int32 &ioNk );
		bool	IsFunctionR( Float* &outSP );
		bool	IsFunctionR2(	Float* &outS1P,
								Float* &outS2P );
		bool	IsFunctionR3(	Float* &outS1P,
								Float* &outS2P,
								Float* &outS3P );
		bool	IsFunctionR4(	Float* &outS1P,
								Float* &outS2P,
								Float* &outS3P,
								Float* &outS4P );
		bool	IsFunctionC(	Float* &outS1P,
								Float* &outS2P );
		bool	IsFunctionC2(	Float* &outS1P,
								Float* &outS2P,
								Float* &outS3P,
								Float* &outS4P );
		bool	IsFunctionC4(	Float* &outS1P,
								Float* &outS2P,
								Float* &outS3P,
								Float* &outS4P,
								Float* &outS5P,
								Float* &outS6P,
								Float* &outS7P,
								Float* &outS8P );
		Int32	UpdateWindow( void );
		Int32	Copy( TFunction* inFunction );
		void	SwapArrayPointers( TFunction* inFunction );

		Int32	mID;

	private:
		Float*	mArrayP;
		Int32*	mCountP;
		Int8**	mHeadH;
		Int32	mDepth;
};

extern TList<TFunction> *gFunctionList;

#endif
