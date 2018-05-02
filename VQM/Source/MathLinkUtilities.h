// ===========================================================================
//	MathLinkUtilities.h
//	Project: Visual Quantum Mechanics.
//  Supervision: Bernd Thaller.
//  Programming: Manfred Liebmann 1998
// ===========================================================================

#ifndef _H_MathLinkUtilities
#define _H_MathLinkUtilities

#include <stdlib.h>
#include "mathlink.h"
#include "TypeDefinition.h"

#ifndef isfinite
#define isfinite(x) ((x) < HUGE_VAL && (x) > -HUGE_VAL)
#endif


// if set to true, MLErrorReport does nothing:
extern bool mlSuppressErrorReport;

//	Prototypes

//Int32	MLErrorReport(	MLINK inLink, Int8* inMessage );
//RM2018
Int32	MLPrint( MLINK inLink, const char *inMessage );

Int32	MLPrintReal( MLINK inLink, Float number);

Int32	MLErrorReport(	MLINK inLink, const char *inMessage );

Int32	MLGetRealArray2(	MLINK inLink,
							Float* &outArrayP,
							Int32* &outCountP,
							Int32 &outDepth	);
Int32	MLReadList(	MLINK inLink,
					Float* inArrayP,
					Int32* inCountP,
					Int32 inDepth,
					Int32 &ioIteration );

Int32	MLCheckMemoryReserve( MLINK inLink );
Int32	MLGetFunctionObject(	MLINK inLink,
								Int32 &outID );
Int32	MLGetOperatorObject(	MLINK inLink,
								Int32 &outID );
Int32	MLGetWindowObject(	MLINK inLink,
							Int32 &outID );
Int32	MLEventLoop( void );

//	Errorcodes

enum {
	eError = 1, eOK = 0
};

#endif
