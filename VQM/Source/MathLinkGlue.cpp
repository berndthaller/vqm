// ===========================================================================
//	MathLinkGlue.c
//	Project: Visual Quantum Mechanics.
//  Supervision: Bernd Thaller.
//  Programming: Manfred Liebmann 1998
//  RM, 2007-06-28 : MLGetFloat -> MLGetReal32
//
//  Changes by Wolfgang Thaller on 15.08.98:
//
//		- Added QSchroedinger1D glue function
// ===========================================================================
//
//	MathLink glue functions

#include "MathLinkUtilities.h"
#include "MathLinkGlue.h"
#include "TFunction.h"
#include "TOperator.h"
#include "TSchroedinger1D.h"
#include "TSchroedinger2D.h"
#include "TSchroedinger3D.h"
#include "TPauli2D.h"
#include "TPauli3D.h"
#include "TDirac2D.h"
#include "TDirac3D.h"

// ---------------------------------------------------------------------------
//		$ main
// ---------------------------------------------------------------------------
//	Start QuantumKernel

TList<TFunction>	*gFunctionList;
TList<TOperator>	*gOperatorList;


int main(int argc,char *argv[])
{
	const Int32	kListSize = 1024;


	gFunctionList = new TList<TFunction>();
	if( !gFunctionList ) return 1;
	if( !gFunctionList->New(kListSize) ) return 1;

	gOperatorList = new TList<TOperator>();
	if( !gOperatorList ) return 1;
	if( !gOperatorList->New(kListSize) ) return 1;

	return MLMain(argc,argv);
}


// ---------------------------------------------------------------------------
//		$ QNewFunction
// ---------------------------------------------------------------------------
//	Create function object

void	QNewFunction( void )
{
	Int32	ID;
	TFunction	*theFunction;


	theFunction = new TFunction();
	if ( theFunction == NULL ) {
		MLErrorReport(stdlink, "out of memory !");
		return;
	}
	
	if( eError == theFunction->GetArray() ) {
		delete theFunction;
		return;
	}

	ID = gFunctionList->Insert(theFunction);
	if( ID == 0 ) {
		delete theFunction;
		MLErrorReport(stdlink, "ID list is full");
		return;
	}
	theFunction->mID = ID;

	MLPutFunction(stdlink, "QFunctionObject", 1);
	MLPutLongInteger(stdlink, ID);
}



// ---------------------------------------------------------------------------
//		$ QDisposeFunction
// ---------------------------------------------------------------------------
//	Remove function object

void	QDisposeFunction( void )
{
	TFunction	*theFunction;
	Int32	ID;


	if( eError == MLGetFunctionObject(stdlink, ID) ) return;

	theFunction = gFunctionList->Remove(ID);
	
	if( !theFunction ) {
		MLErrorReport(stdlink, "invalid function ID");
		return;
	}

	delete theFunction;
	MLPutSymbol(stdlink, "Null");
}



// ---------------------------------------------------------------------------
//  RM2018:	$ QGetArray changed to QGetList
//  define QGetArray now in QuantumKernel.m
// ---------------------------------------------------------------------------
//	Return flat array

void	QGetList( void )
{
	TFunction	*theFunction;
 	Int32	ID;
 
 	
 	if( eError == MLGetFunctionObject(stdlink, ID) ) return;
 
 	theFunction = gFunctionList->Fetch(ID);
 	if( !theFunction ) {
 		MLErrorReport(stdlink, "invalid function ID");
 		return;
 	}
 	
//RM2018: This is now a flat array ...
 	theFunction->PutArray();
}



// ---------------------------------------------------------------------------
//		$ QGetFunctionInfo
// ---------------------------------------------------------------------------
//	Return function info

void	QGetFunctionInfo( void )
{
	TFunction	*theFunction;
	Int32	ID;

	
	if( eError == MLGetFunctionObject(stdlink, ID) ) return;

	theFunction = gFunctionList->Fetch(ID);
	if( !theFunction ) {
		MLErrorReport(stdlink, "invalid function ID");
		return;
	}

	theFunction->PutInfo();
}



// ---------------------------------------------------------------------------
//		$ QGetColorArray
// ---------------------------------------------------------------------------
//	Return RGB color array
//
void	QGetColorArray( void )
{
 	TFunction	*theFunction;
 	Int32	ID;
 
 	
 	if( eError == MLGetFunctionObject(stdlink, ID) ) return;
 
 	theFunction = gFunctionList->Fetch(ID);
 	if( !theFunction ) {
 		MLErrorReport(stdlink, "invalid function ID");
 		return;
 	}
 	
 	theFunction->PutColor();
}



// ---------------------------------------------------------------------------
//		$ QGetGrayArray
// ---------------------------------------------------------------------------
//	Return gray array
//RM2018COMMENT

void	QGetGrayArray( void )
{
	TFunction	*theFunction;
	Int32	ID;

	
	if( eError == MLGetFunctionObject(stdlink, ID) ) return;

	theFunction = gFunctionList->Fetch(ID);
	if( !theFunction ) {
		MLErrorReport(stdlink, "invalid function ID");
		return;
	}
	
	theFunction->PutGray();
}



// ---------------------------------------------------------------------------
//		$ QGetRedBlueArray
// ---------------------------------------------------------------------------
//	Return red-blue color array
//RM2018COMMENT

void	QGetRedBlueArray( void )
{
	TFunction	*theFunction;
 	Int32	ID;
 
 	
 	if( eError == MLGetFunctionObject(stdlink, ID) ) return;
 
 	theFunction = gFunctionList->Fetch(ID);
 	if( !theFunction ) {
 		MLErrorReport(stdlink, "invalid function ID");
 		return;
 	}
 	
 	theFunction->PutRedBlue();
}
 


// ---------------------------------------------------------------------------
//		$ QGetBlackWhiteArray
// ---------------------------------------------------------------------------
//	Return black-white array
//
void	QGetBlackWhiteArray( void )
{
 	TFunction	*theFunction;
 	Int32	ID;
 
 	
 	if( eError == MLGetFunctionObject(stdlink, ID) ) return;
 
 	theFunction = gFunctionList->Fetch(ID);
 	if( !theFunction ) {
 		MLErrorReport(stdlink, "invalid function ID");
 		return;
 	}
 	
 	theFunction->PutBlackWhite();
}



// ---------------------------------------------------------------------------
//		$ QGetAbsArray  defined in QuantumKernel.m
// ---------------------------------------------------------------------------
//	Return abs array

//RM2018
void	QGetAbsList( void )
{
	TFunction	*theFunction;
	Int32	ID;

	
	if( eError == MLGetFunctionObject(stdlink, ID) ) return;
	
	theFunction = gFunctionList->Fetch(ID);
	if( !theFunction ) {
		MLErrorReport(stdlink, "invalid function ID");
		return;
	}
	
	MLPrint(stdlink, "!!!!!!!!!!!!!! in QGetAbsList, vor PutAbs() ... ");

	theFunction->PutAbs();
}



// ---------------------------------------------------------------------------
//		$ QInfo
// ---------------------------------------------------------------------------
//	Return function, operator and window info

void	QInfo( void )
{
	Int32	len, size, ID;
	TFunction	*theFunction;
	TOperator	*theOperator;

	MLPutFunction(stdlink, "List", 3);

	size = gFunctionList->GetSize();
	len = 0;
	for( ID = 0; ID < size; ID++ ) {
		theFunction = gFunctionList->Fetch(ID);
		if( theFunction ) len++;
	}
	MLPutFunction(stdlink, "List", len);
	for( ID = 0; ID < size; ID++ ) {
		theFunction = gFunctionList->Fetch(ID);
		if( theFunction ) {
			MLPutFunction(stdlink, "List", 2);
			MLPutFunction(stdlink, "QFunctionObject", 1);
			MLPutLongInteger(stdlink, ID);
			theFunction->PutInfo();
		}
	}	

	size = gOperatorList->GetSize();
	len = 0;
	for( ID = 0; ID < size; ID++ ) {
		theOperator = gOperatorList->Fetch(ID);
		if( theOperator ) len++;
	}
	MLPutFunction(stdlink, "List", len);
	for( ID = 0; ID < size; ID++ ) {
		theOperator = gOperatorList->Fetch(ID);
		if( theOperator ) {
			MLPutFunction(stdlink, "List", 2);
			MLPutFunction(stdlink, "QOperatorObject", 1);
			MLPutLongInteger(stdlink, ID);
			theOperator->PutInfo();
		}
	}	

	MLPutFunction(stdlink, "List", 0);
}

/*

// ---------------------------------------------------------------------------
//		$ QSchroedinger1D
// ---------------------------------------------------------------------------
//	Create operator object

void	QSchroedinger1D( void )
{
	Int32	scalarID, ID;
	Float	mass, units;
	TOperator	*theOperator;


	if( eError == MLGetFunctionObject(stdlink, scalarID) ) return;
	
	MLGetReal32(stdlink, &mass);
	MLGetReal32(stdlink, &units);
	if( MLError(stdlink) ) {
		MLErrorReport(stdlink, "real numbers for mass and units expected");
		return;
	}

	theOperator = new TSchroedinger1D(scalarID, mass, units);
	if( theOperator == NULL ) {
		MLErrorReport(stdlink, "out of memory");
		return;
	}

	ID = gOperatorList->Insert(theOperator);
	if( ID == 0 ) {
		delete theOperator;
		MLErrorReport(stdlink, "ID list is full");
		return;
	}
	theOperator->mID = ID;

	MLPutFunction(stdlink, "QOperatorObject", 1);
	MLPutLongInteger(stdlink, ID);
}

// ---------------------------------------------------------------------------
//		$ QSchroedinger2D
// ---------------------------------------------------------------------------
//	Create operator object

void	QSchroedinger2D( void )
{
	Int32	scalarID, vectorID, domainID, ID;
	Float	mass, charge, units;
	TOperator	*theOperator;


	if( eError == MLGetFunctionObject(stdlink, scalarID) ) return;
	if( eError == MLGetFunctionObject(stdlink, vectorID) ) return;
	if( eError == MLGetFunctionObject(stdlink, domainID) ) return;
	
	MLGetReal32(stdlink, &mass);
	MLGetReal32(stdlink, &charge);
	MLGetReal32(stdlink, &units);
	if( MLError(stdlink) ) {
		MLErrorReport(stdlink, "real numbers for mass, charge and units expected");
		return;
	}

	theOperator = new TSchroedinger2D(scalarID, vectorID, domainID, mass, charge, units);
	if( theOperator == NULL ) {
		MLErrorReport(stdlink, "out of memory");
		return;
	}

	ID = gOperatorList->Insert(theOperator);
	if( ID == 0 ) {
		delete theOperator;
		MLErrorReport(stdlink, "ID list is full");
		return;
	}
	theOperator->mID = ID;

	MLPutFunction(stdlink, "QOperatorObject", 1);
	MLPutLongInteger(stdlink, ID);
}



// ---------------------------------------------------------------------------
//		$ QSchroedinger3D
// ---------------------------------------------------------------------------
//	Create operator object

void	QSchroedinger3D( void )
{
	Int32	scalarID, vectorID, domainID, ID;
	Float	mass, charge, units;
	TOperator	*theOperator;


	if( eError == MLGetFunctionObject(stdlink, scalarID) ) return;
	if( eError == MLGetFunctionObject(stdlink, vectorID) ) return;
	if( eError == MLGetFunctionObject(stdlink, domainID) ) return;
	
	MLGetReal32(stdlink, &mass);
	MLGetReal32(stdlink, &charge);
	MLGetReal32(stdlink, &units);
	if( MLError(stdlink) ) {
		MLErrorReport(stdlink, "real numbers for mass, charge and units expected");
		return;
	}

	theOperator = new TSchroedinger3D(scalarID, vectorID, domainID, mass, charge, units);
	if( theOperator == NULL ) {
		MLErrorReport(stdlink, "out of memory");
		return;
	}

	ID = gOperatorList->Insert(theOperator);
	if( ID == 0 ) {
		delete theOperator;
		MLErrorReport(stdlink, "ID list is full");
		return;
	}
	theOperator->mID = ID;

	MLPutFunction(stdlink, "QOperatorObject", 1);
	MLPutLongInteger(stdlink, ID);
}



// ---------------------------------------------------------------------------
//		$ QPauli2D
// ---------------------------------------------------------------------------
//	Create operator object

void	QPauli2D( void )
{
	Int32	scalarID, vectorID, domainID, ID;
	Float	mass, charge, units;
	TOperator	*theOperator;


	if( eError == MLGetFunctionObject(stdlink, scalarID) ) return;
	if( eError == MLGetFunctionObject(stdlink, vectorID) ) return;
	if( eError == MLGetFunctionObject(stdlink, domainID) ) return;
	
	MLGetReal32(stdlink, &mass);
	MLGetReal32(stdlink, &charge);
	MLGetReal32(stdlink, &units);
	if( MLError(stdlink) ) {
		MLErrorReport(stdlink, "real numbers for mass, charge and units expected");
		return;
	}

	theOperator = new TPauli2D(scalarID, vectorID, domainID, mass, charge, units);
	if( theOperator == NULL ) {
		MLErrorReport(stdlink, "out of memory");
		return;
	}

	ID = gOperatorList->Insert(theOperator);
	if( ID == 0 ) {
		delete theOperator;
		MLErrorReport(stdlink, "ID list is full");
		return;
	}
	theOperator->mID = ID;

	MLPutFunction(stdlink, "QOperatorObject", 1);
	MLPutLongInteger(stdlink, ID);
}



// ---------------------------------------------------------------------------
//		$ QPauli3D
// ---------------------------------------------------------------------------
//	Create operator object

void	QPauli3D( void )
{
	Int32	scalarID, vectorID, domainID, ID;
	Float	mass, charge, units;
	TOperator	*theOperator;


	if( eError == MLGetFunctionObject(stdlink, scalarID) ) return;
	if( eError == MLGetFunctionObject(stdlink, vectorID) ) return;
	if( eError == MLGetFunctionObject(stdlink, domainID) ) return;
	
	MLGetReal32(stdlink, &mass);
	MLGetReal32(stdlink, &charge);
	MLGetReal32(stdlink, &units);
	if( MLError(stdlink) ) {
		MLErrorReport(stdlink, "real numbers for mass, charge and units expected");
		return;
	}

	theOperator = new TPauli3D(scalarID, vectorID, domainID, mass, charge, units);
	if( theOperator == NULL ) {
		MLErrorReport(stdlink, "out of memory");
		return;
	}

	ID = gOperatorList->Insert(theOperator);
	if( ID == 0 ) {
		delete theOperator;
		MLErrorReport(stdlink, "ID list is full");
		return;
	}
	theOperator->mID = ID;

	MLPutFunction(stdlink, "QOperatorObject", 1);
	MLPutLongInteger(stdlink, ID);
}



// ---------------------------------------------------------------------------
//		$ QDirac2D
// ---------------------------------------------------------------------------
//	Create operator object

void	QDirac2D( void )
{
	Int32	scalarID, vectorID, domainID, ID;
	Float	mass, charge, units;
	TOperator	*theOperator;


	if( eError == MLGetFunctionObject(stdlink, scalarID) ) return;
	if( eError == MLGetFunctionObject(stdlink, vectorID) ) return;
	if( eError == MLGetFunctionObject(stdlink, domainID) ) return;
	
	MLGetReal32(stdlink, &mass);
	MLGetReal32(stdlink, &charge);
	MLGetReal32(stdlink, &units);

	if( MLError(stdlink) ) {
		MLErrorReport(stdlink, "real numbers for mass, charge and units expected");
		return;
	}

	theOperator = new TDirac2D(scalarID, vectorID, domainID, mass, charge, units);
	if( theOperator == NULL ) {
		MLErrorReport(stdlink, "out of memory");
		return;
	}

	ID = gOperatorList->Insert(theOperator);
	if( ID == 0 ) {
		delete theOperator;
		MLErrorReport(stdlink, "ID list is full");
		return;
	}
	theOperator->mID = ID;

	MLPutFunction(stdlink, "QOperatorObject", 1);
	MLPutLongInteger(stdlink, ID);
}



// ---------------------------------------------------------------------------
//		$ QDirac3D
// ---------------------------------------------------------------------------
//	Create operator object

void	QDirac3D( void )
{
	Int32	scalarID, vectorID, domainID, ID;
	Float	mass, charge, units;
	TOperator	*theOperator;


	if( eError == MLGetFunctionObject(stdlink, scalarID) ) return;
	if( eError == MLGetFunctionObject(stdlink, vectorID) ) return;
	if( eError == MLGetFunctionObject(stdlink, domainID) ) return;
	
	MLGetReal32(stdlink, &mass);
	MLGetReal32(stdlink, &charge);
	MLGetReal32(stdlink, &units);

	if( MLError(stdlink) ) {
		MLErrorReport(stdlink, "real numbers for mass, charge and units expected");
		return;
	}

	theOperator = new TDirac3D(scalarID, vectorID, domainID, mass, charge, units);
	if( theOperator == NULL ) {
		MLErrorReport(stdlink, "out of memory");
		return;
	}

	ID = gOperatorList->Insert(theOperator);
	if( ID == 0 ) {
		delete theOperator;
		MLErrorReport(stdlink, "ID list is full");
		return;
	}
	theOperator->mID = ID;

	MLPutFunction(stdlink, "QOperatorObject", 1);
	MLPutLongInteger(stdlink, ID);
}



// ---------------------------------------------------------------------------
//		$ QDisposeOperator
// ---------------------------------------------------------------------------
//	Remove operator object

void	QDisposeOperator( void )
{
	TOperator	*theOperator;
	Int32	ID;

	
	if( eError == MLGetOperatorObject(stdlink, ID) ) return;
	
	theOperator = gOperatorList->Remove(ID);
	if( !theOperator ) {
		MLErrorReport(stdlink, "invalid operator ID");
		return;
	}

	delete theOperator;
	MLPutSymbol(stdlink, "Null");
}



// ---------------------------------------------------------------------------
//		$ QTimeEvolution
// ---------------------------------------------------------------------------
//	Calculate time evolution

void	QTimeEvolution( void )
{
	Int32	operatorID, functionID, fractal, steps;
	Float	timeStep;
	TOperator	*theOperator;
	TFunction	*theFunction;
	

	if( eError == MLGetOperatorObject(stdlink, operatorID) ) return;
	if( eError == MLGetFunctionObject(stdlink, functionID) ) return;

	if( !MLGetReal32(stdlink, &timeStep) ) {
		MLErrorReport(stdlink, "real number for timestep expected");
		return;
	}
	if( !MLGetLongInteger(stdlink, &fractal) ) {
		MLErrorReport(stdlink, "integer for fractal expected");
		return;
	}
	if( !MLGetLongInteger(stdlink, &steps) ) {
		MLErrorReport(stdlink, "integer for steps expected");
		return;
	}

	theOperator = gOperatorList->Fetch(operatorID);
	if( !theOperator ) {
		MLErrorReport(stdlink, "invalid operator ID");
		return;
	}

	theFunction = gFunctionList->Fetch(functionID);
	if( !theFunction ) {
		MLErrorReport(stdlink, "invalid function ID");
		return;
	}

	theOperator->TimeEvolution(theFunction, timeStep, fractal, steps);
}

*/


// ---------------------------------------------------------------------------
//		$ QGetOperatorInfo
// ---------------------------------------------------------------------------
//	Return operator info

void	QGetOperatorInfo( void )
{
	TOperator	*theOperator;
	Int32	ID;


	if( eError == MLGetOperatorObject(stdlink, ID) ) return;

	theOperator = gOperatorList->Fetch(ID);
	if( !theOperator ) {
		MLErrorReport(stdlink, "invalid operator ID");
		return;
	}

	theOperator->PutInfo();
}
