// ===========================================================================
//	TDirac3D.cp
//	Project: Visual Quantum Mechanics.
//  Supervision: Bernd Thaller.
//  Programming: Manfred Liebmann 1998
//  RM: 2007-06-28: change MLPutDouble to MLPutReal32
// ===========================================================================
//
//	Class for Dirac operators

#include "MathLinkUtilities.h"
#include "TDirac3D.h"
#include "TDomain.h"



// ---------------------------------------------------------------------------
//		$ TDirac3D
// ---------------------------------------------------------------------------
//	Constructor

TDirac3D::TDirac3D(	Int32 inScalarID,
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
//		$ ~TDirac3D
// ---------------------------------------------------------------------------
//	Destructor

TDirac3D::~TDirac3D( void )
{

}



// ---------------------------------------------------------------------------
//		$ PutInfo
// ---------------------------------------------------------------------------
//	PutInfo

Int32	TDirac3D::PutInfo( void )
{	
	MLPutFunction(stdlink, "List", 7);
	MLPutFunction(stdlink, "Rule", 2);
	MLPutSymbol(stdlink, "Type");	
	MLPutSymbol(stdlink, "Dirac3D");

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

Int32	TDirac3D::TimeEvolution(	TFunction* inFunction,
									Float inTimeStep,
									Int32 inFractal,
									Int32 inSteps )
{
	Int32	ni, nj, nk, i, n, s;
	Float	*rePsi1P, *imPsi1P, *rePsi2P, *imPsi2P, *rePsi3P, *imPsi3P, *rePsi4P, *imPsi4P;
	Float	*rePhi1P, *imPhi1P, *rePhi2P, *imPhi2P, *rePhi3P, *imPhi3P, *rePhi4P, *imPhi4P;
	Float	*v0P, *w0P, *a1P, *a2P, *a3P, *a4P, *d0P, *tempP, re, im;
	TFunction	*scalarP, *vectorP, *domainP, tempFunction;
	TDomain	tempDomain;
	Int8	*domP;


	if( inFractal < 0 || inFractal >= 16 )
		return MLErrorReport(stdlink, "FRACTAL order is out of range");

	ni = nj = nk = 0;
	if( !inFunction->IsFunction3D(ni, nj, nk) )
		return MLErrorReport(stdlink, "three-dimensional wavefunction expected");
	if( !inFunction->IsFunctionC4(	rePsi1P, imPsi1P, rePsi2P, imPsi2P,
									rePsi3P, imPsi3P, rePsi4P, imPsi4P) )
		return MLErrorReport(stdlink, "bispinor wavefunction expected");

	v0P = w0P = a1P = a2P = a3P = a4P = d0P = NULL;
	if( mScalarID ) {	//	scalar potential
		scalarP = gFunctionList->Fetch(mScalarID);
		if( !scalarP )
			return MLErrorReport(stdlink, "invalid function ID for scalar potential");
		if( !scalarP->IsFunction3D(ni, nj, nk) )
			return MLErrorReport(stdlink, "scalar potential is not compatible");
		if( !scalarP->IsFunctionC(v0P, w0P) )
			return MLErrorReport(stdlink, "scalar potential is not C-valued");
	}
	if( mVectorID ) {	//	vector potential
		vectorP = gFunctionList->Fetch(mVectorID);
		if( !vectorP )
			return MLErrorReport(stdlink, "invalid function ID for vector potential");
		if( !vectorP->IsFunction3D(ni, nj, nk) )
			return MLErrorReport(stdlink, "vector potential is not compatible");
		if( !vectorP->IsFunctionR4(a1P, a2P, a3P, a4P) )
			return MLErrorReport(stdlink, "vector potential is not R4-valued");
	}
	if( mDomainID ) {	//	domain function
		domainP = gFunctionList->Fetch(mDomainID);
		if( !domainP )
			return MLErrorReport(stdlink, "invalid function ID for domain function");
		if( !domainP->IsFunction3D(ni, nj, nk) )
			return MLErrorReport(stdlink, "domain function is not compatible");
		if( !domainP->IsFunctionR(d0P) )
			return MLErrorReport(stdlink, "domain function is not R-valued");
	}

	domP = NULL;
	if( eError == tempDomain.InitPlain3D(domP, d0P, ni, nj, nk) ) return eError;
	tempDomain.ClearC4(	rePsi1P, imPsi1P, rePsi2P, imPsi2P,
						rePsi3P, imPsi3P, rePsi4P, imPsi4P, ni*nj*nk);

	if( eError == tempFunction.Copy(inFunction) ) return eError;
	tempFunction.IsFunctionC4(	rePhi1P, imPhi1P, rePhi2P, imPhi2P,
								rePhi3P, imPhi3P, rePhi4P, imPhi4P);

	for( s = 0; s < inSteps; s++ ) {
		n = 1 << inFractal >> 1;	//fractal iteration
		for( i = 0; i < n; i++ ) {
			GetFractalNumber( inFractal, i, re, im );				
			re *= inTimeStep;
			im *= inTimeStep;
			Kernel(	rePsi1P, imPsi1P, rePsi2P, imPsi2P, rePsi3P, imPsi3P, rePsi4P, imPsi4P,
					rePhi1P, imPhi1P, rePhi2P, imPhi2P, rePhi3P, imPhi3P, rePhi4P, imPhi4P,
					v0P, w0P, a1P, a2P, a3P, a4P, domP, re, im, ni, nj, nk );

			tempP = rePsi1P; rePsi1P = rePhi1P; rePhi1P = tempP;
			tempP = imPsi1P; imPsi1P = imPhi1P; imPhi1P = tempP;
			tempP = rePsi2P; rePsi2P = rePhi2P; rePhi2P = tempP;
			tempP = imPsi2P; imPsi2P = imPhi2P; imPhi2P = tempP;
			tempP = rePsi3P; rePsi3P = rePhi3P; rePhi3P = tempP;
			tempP = imPsi3P; imPsi3P = imPhi3P; imPhi3P = tempP;
			tempP = rePsi4P; rePsi4P = rePhi4P; rePhi4P = tempP;
			tempP = imPsi4P; imPsi4P = imPhi4P; imPhi4P = tempP;

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

void	TDirac3D::Kernel(	Float* rePsi1P,
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
							Int32 nk )
{
	Float	rePsi1C, imPsi1C, reEta1C, imEta1C;
	Float	rePsi1R, imPsi1R, rePsi1L, imPsi1L;
	Float	rePsi1U, imPsi1U, rePsi1D, imPsi1D;
	Float	rePsi1F, imPsi1F, rePsi1B, imPsi1B;
	Float	rePsi2C, imPsi2C, reEta2C, imEta2C;
	Float	rePsi2R, imPsi2R, rePsi2L, imPsi2L;
	Float	rePsi2U, imPsi2U, rePsi2D, imPsi2D;
	Float	rePsi2F, imPsi2F, rePsi2B, imPsi2B;
	Float	rePsi3C, imPsi3C, reEta3C, imEta3C;
	Float	rePsi3R, imPsi3R, rePsi3L, imPsi3L;
	Float	rePsi3U, imPsi3U, rePsi3D, imPsi3D;
	Float	rePsi3F, imPsi3F, rePsi3B, imPsi3B;
	Float	rePsi4C, imPsi4C, reEta4C, imEta4C;
	Float	rePsi4R, imPsi4R, rePsi4L, imPsi4L;
	Float	rePsi4U, imPsi4U, rePsi4D, imPsi4D;
	Float	rePsi4F, imPsi4F, rePsi4B, imPsi4B;
	Float	v0C, w0C, a1C, a2C, a3C, a4C;
	Float	ch, m, e;
	Float	one = 1.0, two = 2.0;
	Int32	i, mode = 0;


	ch = one / ( two * mUnits );
	m = mMass;
	e = mCharge;

	if( mScalarID ) mode |= kScalar;
	if( mVectorID ) mode |= kVector;

	for( i = 0; i < ni*nj*nk; i++ ) {
		if( *domP++ ) {
			rePsi1R = *(rePsi1P + 1);
			imPsi1R = *(imPsi1P + 1);
			rePsi1L = *(rePsi1P +-1);
			imPsi1L = *(imPsi1P +-1);
			rePsi1U = *(rePsi1P + ni);
			imPsi1U = *(imPsi1P + ni);
			rePsi1D = *(rePsi1P +-ni);
			imPsi1D = *(imPsi1P +-ni);
			rePsi1F = *(rePsi1P + ni*nj);
			imPsi1F = *(imPsi1P + ni*nj);
			rePsi1B = *(rePsi1P +-ni*nj);
			imPsi1B = *(imPsi1P +-ni*nj);
			rePsi2R = *(rePsi2P + 1);
			imPsi2R = *(imPsi2P + 1);
			rePsi2L = *(rePsi2P +-1);
			imPsi2L = *(imPsi2P +-1);
			rePsi2U = *(rePsi2P + ni);
			imPsi2U = *(imPsi2P + ni);
			rePsi2D = *(rePsi2P +-ni);
			imPsi2D = *(imPsi2P +-ni);
			rePsi2F = *(rePsi2P + ni*nj);
			imPsi2F = *(imPsi2P + ni*nj);
			rePsi2B = *(rePsi2P +-ni*nj);
			imPsi2B = *(imPsi2P +-ni*nj);
			reEta3C = ch * (imPsi2R - imPsi2L - rePsi2U + rePsi2D + imPsi1F - imPsi1B);
			imEta3C = ch * (rePsi2L - rePsi2R - imPsi2U + imPsi2D + rePsi1B - rePsi1F);
			reEta4C = ch * (imPsi1R - imPsi1L + rePsi1U - rePsi1D + imPsi2B - imPsi2F);
			imEta4C = ch * (rePsi1L - rePsi1R + imPsi1U - imPsi1D + rePsi2F - rePsi2B);
			rePsi3R = *(rePsi3P + 1);
			imPsi3R = *(imPsi3P + 1);
			rePsi3L = *(rePsi3P +-1);
			imPsi3L = *(imPsi3P +-1);
			rePsi3U = *(rePsi3P + ni);
			imPsi3U = *(imPsi3P + ni);
			rePsi3D = *(rePsi3P +-ni);
			imPsi3D = *(imPsi3P +-ni);
			rePsi3F = *(rePsi3P + ni*nj);
			imPsi3F = *(imPsi3P + ni*nj);
			rePsi3B = *(rePsi3P +-ni*nj);
			imPsi3B = *(imPsi3P +-ni*nj);
			rePsi4R = *(rePsi4P + 1);
			imPsi4R = *(imPsi4P + 1);
			rePsi4L = *(rePsi4P +-1);
			imPsi4L = *(imPsi4P +-1);
			rePsi4U = *(rePsi4P + ni);
			imPsi4U = *(imPsi4P + ni);
			rePsi4D = *(rePsi4P +-ni);
			imPsi4D = *(imPsi4P +-ni);
			rePsi4F = *(rePsi4P + ni*nj);
			imPsi4F = *(imPsi4P + ni*nj);
			rePsi4B = *(rePsi4P +-ni*nj);
			imPsi4B = *(imPsi4P +-ni*nj);
			reEta1C = ch * (imPsi4R - imPsi4L - rePsi4U + rePsi4D + imPsi3F - imPsi3B);
			imEta1C = ch * (rePsi4L - rePsi4R - imPsi4U + imPsi4D + rePsi3B - rePsi3F);
			reEta2C = ch * (imPsi3R - imPsi3L + rePsi3U - rePsi3D + imPsi4B - imPsi4F);
			imEta2C = ch * (rePsi3L - rePsi3R + imPsi3U - imPsi3D + rePsi4F - rePsi4B);
			rePsi1C = *rePsi1P++;
			imPsi1C = *imPsi1P++;
			rePsi2C = *rePsi2P++;
			imPsi2C = *imPsi2P++;
			rePsi3C = *rePsi3P++;
			imPsi3C = *imPsi3P++;
			rePsi4C = *rePsi4P++;
			imPsi4C = *imPsi4P++;
			if( mode & kVector ) {
				a1C = *a1P++;
				a2C = *a2P++;
				a3C = *a3P++;
				a4C = *a4P++;
				reEta1C -= e * a1C * rePsi4C;
				imEta1C -= e * a1C * imPsi4C;
				reEta2C -= e * a1C * rePsi3C;
				imEta2C -= e * a1C * imPsi3C;
				reEta3C -= e * a1C * rePsi2C;
				imEta3C -= e * a1C * imPsi2C;
				reEta4C -= e * a1C * rePsi1C;
				imEta4C -= e * a1C * imPsi1C;
				reEta1C -= e * a2C * imPsi4C;
				imEta1C += e * a2C * rePsi4C;
				reEta2C += e * a2C * imPsi3C;
				imEta2C -= e * a2C * rePsi3C;
				reEta3C -= e * a2C * imPsi2C;
				imEta3C += e * a2C * rePsi2C;
				reEta4C += e * a2C * imPsi1C;
				imEta4C -= e * a2C * rePsi1C;
				reEta1C -= e * a3C * rePsi3C;
				imEta1C -= e * a3C * imPsi3C;
				reEta2C += e * a3C * rePsi4C;
				imEta2C += e * a3C * imPsi4C;
				reEta3C -= e * a3C * rePsi1C;
				imEta3C -= e * a3C * imPsi1C;
				reEta4C += e * a3C * rePsi2C;
				imEta4C += e * a3C * imPsi2C;
				reEta1C -= e * a4C * rePsi1C;
				imEta1C -= e * a4C * imPsi1C;
				reEta2C -= e * a4C * rePsi2C;
				imEta2C -= e * a4C * imPsi2C;
				reEta3C += e * a4C * rePsi3C;
				imEta3C += e * a4C * imPsi3C;
				reEta4C += e * a4C * rePsi4C;
				imEta4C += e * a4C * imPsi4C;
			}
			reEta1C += m * rePsi1C;
			imEta1C += m * imPsi1C;
			reEta2C += m * rePsi2C;
			imEta2C += m * imPsi2C;
			reEta3C -= m * rePsi3C;
			imEta3C -= m * imPsi3C;
			reEta4C -= m * rePsi4C;
			imEta4C -= m * imPsi4C;
			if( mode & kScalar ) {
				v0C = *v0P++;
				w0C = *w0P++;
				reEta1C += e * v0C * rePsi1C;
				imEta1C += e * v0C * imPsi1C;
				reEta2C += e * v0C * rePsi2C;
				imEta2C += e * v0C * imPsi2C;
				reEta3C += e * v0C * rePsi3C;
				imEta3C += e * v0C * imPsi3C;
				reEta4C += e * v0C * rePsi4C;
				imEta4C += e * v0C * imPsi4C;
				reEta1C -= e * w0C * imPsi1C;
				imEta1C += e * w0C * rePsi1C;
				reEta2C -= e * w0C * imPsi2C;
				imEta2C += e * w0C * rePsi2C;
				reEta3C -= e * w0C * imPsi3C;
				imEta3C += e * w0C * rePsi3C;
				reEta4C -= e * w0C * imPsi4C;
				imEta4C += e * w0C * rePsi4C;
			}
			*rePhi1P++ = rePsi1C + imZ * reEta1C + reZ * imEta1C;
			*imPhi1P++ = imPsi1C - reZ * reEta1C + imZ * imEta1C;
			*rePhi2P++ = rePsi2C + imZ * reEta2C + reZ * imEta2C;
			*imPhi2P++ = imPsi2C - reZ * reEta2C + imZ * imEta2C;
			*rePhi3P++ = rePsi3C + imZ * reEta3C + reZ * imEta3C;
			*imPhi3P++ = imPsi3C - reZ * reEta3C + imZ * imEta3C;
			*rePhi4P++ = rePsi4C + imZ * reEta4C + reZ * imEta4C;
			*imPhi4P++ = imPsi4C - reZ * reEta4C + imZ * imEta4C;
		}
		else {			
			rePsi1P++;
			imPsi1P++;
			rePsi2P++;
			imPsi2P++;
			rePsi3P++;
			imPsi3P++;
			rePsi4P++;
			imPsi4P++;
			if( mode & kScalar ) { v0P++; w0P++; }
			rePhi1P++;
			imPhi1P++;
			rePhi2P++;
			imPhi2P++;
			rePhi3P++;
			imPhi3P++;
			rePhi4P++;
			imPhi4P++;
			if( mode & kVector ) { a1P++; a2P++; a3P++; a4P++; }
		}
	}

}

