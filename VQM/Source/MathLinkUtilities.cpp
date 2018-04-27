// ===========================================================================
//	MathLinkUtilities.c
//	Project: Visual Quantum Mechanics.
//  Supervision: Bernd Thaller.
//  Programming: Manfred Liebmann 1998
// ===========================================================================
//
//	MathLink utility functions

#include "MathLinkUtilities.h"
#include <string.h>
#include <stdio.h>


	// if set to true, MLErrorReport does nothing:
bool mlSuppressErrorReport = false;

// ---------------------------------------------------------------------------
//		$ MLErrorReport
// ---------------------------------------------------------------------------
//	Send error message

Int32	MLErrorReport(	MLINK inLink,
						Int8* inMessage )
{
	if(mlSuppressErrorReport)
		return eError;
		
	char	errMsg[256]; 	

	
	sprintf(errMsg, "%s\"%.192s\"%s","Message[QuantumKernel::err,",inMessage,"]");
	MLClearError(inLink);
	MLNewPacket(inLink);
	MLEvaluate(inLink, errMsg);
	while( MLNextPacket(inLink) != RETURNPKT ) MLNewPacket(inLink);
	MLNewPacket(inLink);
	MLPutSymbol(inLink, "$Failed");

	return eError;		//we have an error!
}



// ---------------------------------------------------------------------------
//		$ MLGetRealArray2
// ---------------------------------------------------------------------------
//	Substitution for MLGetRealArray

Int32	MLGetRealArray2(	MLINK inLink,
							Float* &outArrayP,
							Int32* &outCountP,
							Int32 &outDepth	)
{
	Int32	i, len, size, ioIteration;
	const Int32	kIterationLimit = 16;	

	outCountP = new Int32[kIterationLimit];
	if( outCountP == NULL )
		return MLErrorReport(inLink, "out of memory");

	outDepth = 0;
	while( MLGetType(inLink) == MLTKFUNC ) {
		MLCheckFunction(inLink, "List", &len);
		if( MLError(inLink) || len == 0 )
			return MLErrorReport(inLink, "out of sequence");
		outCountP[outDepth] = len;
		outDepth++;
		if( outDepth >= kIterationLimit )
			return MLErrorReport(inLink, "iteration limit reached");		
	}
	if( outDepth == 0 ) return MLErrorReport(inLink, "out of sequence");
	
	size = 1;
	for( i = 0; i < outDepth; i++ ) size *= outCountP[i];
	outArrayP = new Float[size];
	if( outArrayP == NULL )
		return MLErrorReport(inLink, "out of memory");
	if( eError == MLCheckMemoryReserve(inLink) ) return eError;

	if( eError == MLReadList(inLink, outArrayP, outCountP, outDepth, ioIteration = 0) )
		return eError;

	return eOK;
}



// ---------------------------------------------------------------------------
//		$ MLReadList
// ---------------------------------------------------------------------------
//	Iterator for MLGetRealArray2

Int32	MLReadList(	MLINK inLink,
					Float* inArrayP,
					Int32* inCountP,
					Int32 inDepth,
					Int32 &ioIteration )
{
	static Float*	sArrayP;
	Int32	i, len;


	if( ioIteration == 0 ) sArrayP = inArrayP;

	if( sArrayP == inArrayP ) len = inCountP[ioIteration];
	else {
		MLCheckFunction(inLink, "List", &len);
		if(	MLError(inLink) || len != inCountP[ioIteration] )
			return MLErrorReport(inLink, "out of sequence");
	}

	ioIteration++;
	for( i = 0; i < len; i++) {
		if( ioIteration == inDepth ) MLGetFloat(inLink, sArrayP++);
		else {
			if( eError == MLReadList(inLink, inArrayP, inCountP, inDepth, ioIteration) )
				return eError;
		}
	}
	ioIteration--;

	if( MLError(inLink) )
		return MLErrorReport(inLink, "out of sequence");

	return eOK;
}



// ---------------------------------------------------------------------------
//		$ MLCheckMemoryReserve
// ---------------------------------------------------------------------------
//	Check memory reserve

Int32	MLCheckMemoryReserve( MLINK inLink )
{
	Int8	*reserveP[32];
	Int32	i, j;
	const Int32	kMemoryReserve = 32768;		//32 Kbytes memory reserve

	
	for( i = 0; i < 32; i++ ) {					//32 small memory blocks
		reserveP[i] = new Int8[kMemoryReserve/32];
		if( reserveP[i] == NULL ) {
			for( j = 0; j < i; j++ ) delete [] reserveP[j];
			return MLErrorReport(inLink, "low memory situation");
		}
	}
	for( i = 0; i < 32; i++ ) delete [] reserveP[i];
	
	reserveP[0] = new Int8[kMemoryReserve];		//one large memory block
	if( reserveP[0] == NULL )
		return MLErrorReport(inLink, "low memory situation");
	delete [] reserveP[0];

	return eOK;
}



// ---------------------------------------------------------------------------
//		$ MLGetFunctionObject
// ---------------------------------------------------------------------------
//	Get function identification

Int32	MLGetFunctionObject(	MLINK inLink,
								Int32 &outID )
{
	Int32	len, result, ID;
	const char*	symbolP;
	

	switch( MLGetType(inLink) ) {
	case MLTKFUNC:
		MLCheckFunction(inLink, "QFunctionObject", &len);	
		if( MLError(inLink) || len != 1 )
			return MLErrorReport(inLink, "Function object expected");
		if( !MLGetLongInteger(inLink, &ID) )
			return MLErrorReport(inLink, "integer ID expected");
		outID = ID;
		return eOK;
		break;
	case MLTKSYM:
		MLGetSymbol(inLink, &symbolP);
		result = strcmp( symbolP, "None");
//		MLDisownSymbol(inLink, symbolP);
		MLReleaseSymbol(inLink, symbolP);
		if( result == 0 ) {
			outID = 0;
			return eOK;
		}
	default:
		return MLErrorReport(inLink, "Function object or \"None\" expected");
	}
	
}



// ---------------------------------------------------------------------------
//		$ MLGetOperatorObject
// ---------------------------------------------------------------------------
//	Get operator identification

Int32	MLGetOperatorObject(	MLINK inLink,
								Int32 &outID )
{
	Int32	len, ID;


	MLCheckFunction(inLink, "QOperatorObject", &len);	
	if( MLError(inLink) || len != 1 )
		return MLErrorReport(inLink, "Operator object expected");
	if( !MLGetLongInteger(inLink, &ID) )
		return MLErrorReport(inLink, "integer ID expected");

	outID = ID;
	return eOK;
}



