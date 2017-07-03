// ===========================================================================
//	TDirac2D.cp
//	Project: Visual Quantum Mechanics.
//  Supervision: Bernd Thaller.
//  Programming: Manfred Liebmann 1998
//  RM: 2007-06-28: change MLPutDouble to MLPutReal32
// ===========================================================================
//
//	Class for Dirac operators

#include "MathLinkUtilities.h"
#include "TDirac2D.h"
#include "TDomain.h"



// ---------------------------------------------------------------------------
//		$ TDirac2D
// ---------------------------------------------------------------------------
//	Constructor

TDirac2D::TDirac2D(	Int32 inScalarID,
					Int32 inVectorID,
					Int32 inDomainID,
					Float inMass,
					Float inCharge,
					Float inUnits )
{
	mScalarID = inScalarID;
	mVectorID = inVectorID;
	mDomainID = inDomainID;
	mMass = inMass;
	mCharge = inCharge;
	mUnits = inUnits;
}



// ---------------------------------------------------------------------------
//		$ ~TDirac2D
// ---------------------------------------------------------------------------
//	Destructor

TDirac2D::~TDirac2D( void )
{

}



// ---------------------------------------------------------------------------
//		$ PutInfo
// ---------------------------------------------------------------------------
//	PutInfo

Int32	TDirac2D::PutInfo( void )
{	
	MLPutFunction(stdlink, "List", 7);
	MLPutFunction(stdlink, "Rule", 2);
	MLPutSymbol(stdlink, "Type");	
	MLPutSymbol(stdlink, "Dirac2D");

	MLPutFunction(stdlink, "Rule", 2);
	MLPutSymbol(stdlink, "ScalarPotential");
	if( !mScalarID ) MLPutSymbol(stdlink, "None");
	else {
		MLPutFunction(stdlink, "QFunctionObject", 1);
		MLPutInteger(stdlink, mScalarID);
	}
	
	MLPutFunction(stdlink, "Rule", 2);
	MLPutSymbol(stdlink, "VectorPotential");
	if( !mVectorID ) MLPutSymbol(stdlink, "None");
	else {
		MLPutFunction(stdlink, "QFunctionObject", 1);
		MLPutInteger(stdlink, mVectorID);
	}
	
	MLPutFunction(stdlink, "Rule", 2);
	MLPutSymbol(stdlink, "Domain");
	if( !mDomainID ) MLPutSymbol(stdlink, "None");
	else { 
		MLPutFunction(stdlink, "QFunctionObject", 1);
		MLPutInteger(stdlink, mDomainID);
	}
	
	MLPutFunction(stdlink, "Rule", 2);
	MLPutSymbol(stdlink, "Mass");
	MLPutReal32(stdlink, mMass);

	MLPutFunction(stdlink, "Rule", 2);
	MLPutSymbol(stdlink, "Charge");
	MLPutReal32(stdlink, mCharge);

	MLPutFunction(stdlink, "Rule", 2);
	MLPutSymbol(stdlink, "Units");
	MLPutReal32(stdlink, mUnits);

	return eOK;
}



// ---------------------------------------------------------------------------
//		$ TimeEvolution
// ---------------------------------------------------------------------------
//	TimeEvolution

Int32	TDirac2D::TimeEvolution(	TFunction* inFunction,
									Float inTimeStep,
									Int32 inFractal,
									Int32 inSteps )
{
	Int32	ni, nj, i, n, s;
	Float	*rePsi1P, *imPsi1P, *rePsi2P, *imPsi2P, *rePhi1P, *imPhi1P, *rePhi2P, *imPhi2P;
	Float	*v0P, *w0P, *a1P, *a2P, *a3P, *d0P, *tempP, re, im;
	TFunction	*scalarP, *vectorP, *domainP, tempFunction;
	TDomain	tempDomain;
	Int8	*domP;


	if( inFractal < 0 || inFractal >= 16 )
		return MLErrorReport(stdlink, "fractal order is out of range");

	ni = nj = 0;
	if( !inFunction->IsFunction2D(ni, nj) )
		return MLErrorReport(stdlink, "two-dimensional wavefunction expected");
	if( !inFunction->IsFunctionC2(rePsi1P, imPsi1P, rePsi2P, imPsi2P) )
		return MLErrorReport(stdlink, "spinor wavefunction expected");

	v0P = w0P = a1P = a2P = a3P = d0P = NULL;
	if( mScalarID ) {	//	scalar potential
		scalarP = gFunctionList->Fetch(mScalarID);
		if( !scalarP )
			return MLErrorReport(stdlink, "invalid function ID for scalar potential");
		if( !scalarP->IsFunction2D(ni, nj) )
			return MLErrorReport(stdlink, "scalar potential is not compatible");
		if( !scalarP->IsFunctionC(v0P, w0P) )
			return MLErrorReport(stdlink, "scalar potential is not C-valued");
	}
	if( mVectorID ) {	//	vector potential
		vectorP = gFunctionList->Fetch(mVectorID);
		if( !vectorP )
			return MLErrorReport(stdlink, "invalid function ID for vector potential");
		if( !vectorP->IsFunction2D(ni, nj) )
			return MLErrorReport(stdlink, "vector potential is not compatible");
		if( !vectorP->IsFunctionR3(a1P, a2P, a3P) )
			return MLErrorReport(stdlink, "vector potential is not R3-valued");
	}
	if( mDomainID ) {	//	domain function
		domainP = gFunctionList->Fetch(mDomainID);
		if( !domainP )
			return MLErrorReport(stdlink, "invalid function ID for domain function");
		if( !domainP->IsFunction2D(ni, nj) )
			return MLErrorReport(stdlink, "domain function is not compatible");
		if( !domainP->IsFunctionR(d0P) )
			return MLErrorReport(stdlink, "domain function is not R-valued");
	}

	domP = NULL;
	if( eError == tempDomain.InitPlain2D(domP, d0P, ni, nj) ) return eError;
	tempDomain.ClearC2(rePsi1P, imPsi1P, rePsi2P, imPsi2P, ni*nj);

	if( eError == tempFunction.Copy(inFunction) ) return eError;
	tempFunction.IsFunctionC2(rePhi1P, imPhi1P, rePhi2P, imPhi2P);

	for( s = 0; s < inSteps; s++ ) {
		n = 1 << inFractal >> 1;	//fractal iteration
		for( i = 0; i < n; i++ ) {
			GetFractalNumber( inFractal, i, re, im );				
			re *= inTimeStep;
			im *= inTimeStep;
			Kernel(	rePsi1P, imPsi1P, rePsi2P, imPsi2P, rePhi1P, imPhi1P, rePhi2P, imPhi2P,
					v0P, w0P, a1P, a2P, a3P, domP, re, im, ni, nj );

			tempP = rePsi1P; rePsi1P = rePhi1P; rePhi1P = tempP;
			tempP = imPsi1P; imPsi1P = imPhi1P; imPhi1P = tempP;
			tempP = rePsi2P; rePsi2P = rePhi2P; rePhi2P = tempP;
			tempP = imPsi2P; imPsi2P = imPhi2P; imPhi2P = tempP;
		}
		if( n == 1 ) inFunction->SwapArrayPointers(&tempFunction);

		if( eError == inFunction->UpdateWindow() ) return eError;

		if(MLYieldFunction(stdlink))
			MLCallYieldFunction(MLYieldFunction(stdlink), stdlink, (MLYieldParameters)0);
		if(MLAbort) {
			MLPutFunction(stdlink, "Abort", 0);
			return eError;
		}
	}

	MLPutSymbol(stdlink, "Null");
	return eOK;
}



// ---------------------------------------------------------------------------
//		$ Kernel
// ---------------------------------------------------------------------------

void	TDirac2D::Kernel(	Float* rePsi1P,
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
							Int32 nj )
{
	Float	rePsi1C, imPsi1C, reEta1C, imEta1C;
	Float	rePsi1R, imPsi1R, rePsi1L, imPsi1L;
	Float	rePsi1U, imPsi1U, rePsi1D, imPsi1D;
	Float	rePsi2C, imPsi2C, reEta2C, imEta2C;
	Float	rePsi2R, imPsi2R, rePsi2L, imPsi2L;
	Float	rePsi2U, imPsi2U, rePsi2D, imPsi2D;
	Float	v0C, w0C, a1C, a2C, a3C;
	Float	ch, m, e;
	Float	one = 1.0, two = 2.0;
	Int32	i, mode = 0;


	ch = one / ( two * mUnits );
	m = mMass;
	e = mCharge;

	if( mScalarID ) mode |= kScalar;
	if( mVectorID ) mode |= kVector;
	
	for( i = 0; i < ni*nj; i++ ) {
		if( *domP++ ) {
			rePsi1R = *(rePsi1P + 1);
			imPsi1R = *(imPsi1P + 1);
			rePsi1L = *(rePsi1P +-1);
			imPsi1L = *(imPsi1P +-1);
			rePsi1U = *(rePsi1P + ni);
			imPsi1U = *(imPsi1P + ni);
			rePsi1D = *(rePsi1P +-ni);
			imPsi1D = *(imPsi1P +-ni);
			reEta2C = ch * (imPsi1R - imPsi1L + rePsi1U - rePsi1D);
			imEta2C = ch * (rePsi1L - rePsi1R + imPsi1U - imPsi1D);
			rePsi2R = *(rePsi2P + 1);
			imPsi2R = *(imPsi2P + 1);
			rePsi2L = *(rePsi2P +-1);
			imPsi2L = *(imPsi2P +-1);
			rePsi2U = *(rePsi2P + ni);
			imPsi2U = *(imPsi2P + ni);
			rePsi2D = *(rePsi2P +-ni);
			imPsi2D = *(imPsi2P +-ni);
			reEta1C = ch * (imPsi2R - imPsi2L - rePsi2U + rePsi2D);
			imEta1C = ch * (rePsi2L - rePsi2R - imPsi2U + imPsi2D);
			rePsi1C = *rePsi1P++;
			imPsi1C = *imPsi1P++;
			rePsi2C = *rePsi2P++;
			imPsi2C = *imPsi2P++;
			if( mode & kVector ) {
				a1C = *a1P++;
				a2C = *a2P++;
				a3C = *a3P++;
				reEta1C -= e * a1C * rePsi2C;
				imEta1C -= e * a1C * imPsi2C;
				reEta2C -= e * a1C * rePsi1C;
				imEta2C -= e * a1C * imPsi1C;
				reEta1C -= e * a2C * imPsi2C;
				imEta1C += e * a2C * rePsi2C;
				reEta2C += e * a2C * imPsi1C;
				imEta2C -= e * a2C * rePsi1C;
				reEta1C -= e * a3C * rePsi1C;
				imEta1C -= e * a3C * imPsi1C;
				reEta2C += e * a3C * rePsi2C;
				imEta2C += e * a3C * imPsi2C;
			}
			reEta1C += m * rePsi1C;
			imEta1C += m * imPsi1C;
			reEta2C -= m * rePsi2C;
			imEta2C -= m * imPsi2C;
			if( mode & kScalar ) {
				v0C = *v0P++;
				w0C = *w0P++;
				reEta1C += e * v0C * rePsi1C;
				imEta1C += e * v0C * imPsi1C;
				reEta2C += e * v0C * rePsi2C;
				imEta2C += e * v0C * imPsi2C;
				reEta1C -= e * w0C * imPsi1C;
				imEta1C += e * w0C * rePsi1C;
				reEta2C -= e * w0C * imPsi2C;
				imEta2C += e * w0C * rePsi2C;
			}
			*rePhi1P++ = rePsi1C + imZ * reEta1C + reZ * imEta1C;
			*imPhi1P++ = imPsi1C - reZ * reEta1C + imZ * imEta1C;
			*rePhi2P++ = rePsi2C + imZ * reEta2C + reZ * imEta2C;
			*imPhi2P++ = imPsi2C - reZ * reEta2C + imZ * imEta2C;
		}
		else {
			rePsi1P++;
			imPsi1P++;
			rePsi2P++;
			imPsi2P++;
			if( mode & kScalar ) { v0P++; w0P++; }
			rePhi1P++;
			imPhi1P++;
			rePhi2P++;
			imPhi2P++;
			if( mode & kVector ) { a1P++; a2P++; a3P++; }
		}
	}

}


