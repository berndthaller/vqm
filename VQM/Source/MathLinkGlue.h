// ===========================================================================
//	MathLinkGlue.h
//	Project: Visual Quantum Mechanics.
//  Supervision: Bernd Thaller.
//  Programming: Manfred Liebmann 1998
//
// Changes by Wolfgang Thaller on 15.08.98:
//
//		- Added prototype for QSchroedinger1D
// ===========================================================================
//

#ifndef _H_MathLinkGlue
#define _H_MathLinkGlue
#pragma once

#include "TypeDefinition.h"

//	Prototypes

void	QNewFunction( void );
void	QDisposeFunction( void );
//RM2018 void	QGetArray( void );
//RM2018
void	QGetAbsList( void );
void	QGetList( void );
void	QGetFunctionInfo( void );
void	QGetColorArray( void );
void	QGetGrayArray( void );
void	QGetRedBlueArray( void );
void	QGetBlackWhiteArray( void );
void	QGetAbsArray( void );
void	QInfo( void );

void	QSchroedinger1D( void );
void	QSchroedinger2D( void );
void	QSchroedinger3D( void );
void	QPauli2D( void );
void	QPauli3D( void );
void	QDirac2D( void );
void	QDirac3D( void );
void	QDisposeOperator( void );
void	QTimeEvolution( void );
void	QGetOperatorInfo( void );

void	QShowWindow( void );
void	QHideWindow( void );
void	QGetWindowInfo( void );
void	QBeginMovie( void );
void	QEndMovie( void );

#endif
