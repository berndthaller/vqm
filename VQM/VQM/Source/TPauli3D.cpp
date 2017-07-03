// ===========================================================================
//	TPauli3D.cp
//	Project: Visual Quantum Mechanics.
//  Supervision: Bernd Thaller.
//  Programming: Manfred Liebmann 1998
// ===========================================================================
//
//	Class for Pauli operators

#include "MathLinkUtilities.h"
#include "TPauli3D.h"
#include "TDomain.h"



// ---------------------------------------------------------------------------
//		$ TPauli3D
// ---------------------------------------------------------------------------
//	Constructor

TPauli3D::TPauli3D(	Int32 inScalarID,
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
//		$ ~TPauli3D
// ---------------------------------------------------------------------------
//	Destructor

TPauli3D::~TPauli3D( void )
{

}



// ---------------------------------------------------------------------------
//		$ PutInfo
// ---------------------------------------------------------------------------
//	PutInfo

Int32	TPauli3D::PutInfo( void )
{	
	MLPutFunction(stdlink, "List", 7);
	MLPutFunction(stdlink, "Rule", 2);
	MLPutSymbol(stdlink, "Type");	
	MLPutSymbol(stdlink, "Pauli3D");

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
	MLPutDouble(stdlink, mMass);

	MLPutFunction(stdlink, "Rule", 2);
	MLPutSymbol(stdlink, "Charge");
	MLPutDouble(stdlink, mCharge);

	MLPutFunction(stdlink, "Rule", 2);
	MLPutSymbol(stdlink, "Units");
	MLPutDouble(stdlink, mUnits);

	return eOK;
}



// ---------------------------------------------------------------------------
//		$ TimeEvolution
// ---------------------------------------------------------------------------
//	TimeEvolution

Int32	TPauli3D::TimeEvolution(	TFunction* inFunction,
									Float inTimeStep,
									Int32 inFractal,
									Int32 inSteps )
{
	Int32	ni, nj, nk, i, n, s;
	Float	*rePsi1P, *imPsi1P, *rePsi2P, *imPsi2P, *rePhi1P, *imPhi1P, *rePhi2P, *imPhi2P;
	Float	*v0P, *w0P, *a1P, *a2P, *a3P, *d0P, *tempP, re, im;
	TFunction	*scalarP, *vectorP, *domainP, tempFunction;
	TDomain	tempDomain;
	Int8	*domP;


	if( inFractal < 0 || inFractal >= 16 )
		return MLErrorReport(stdlink, "fractal order is out of range");

	ni = nj = nk = 0;
	if( !inFunction->IsFunction3D(ni, nj, nk) )
		return MLErrorReport(stdlink, "three-dimensional wavefunction expected");
	if( !inFunction->IsFunctionC2(rePsi1P, imPsi1P, rePsi2P, imPsi2P) )
		return MLErrorReport(stdlink, "spinor wavefunction expected");

	v0P = w0P = a1P = a2P = a3P = d0P = NULL;
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
		if( !vectorP->IsFunctionR3(a1P, a2P, a3P) )
			return MLErrorReport(stdlink, "vector potential is not R3-valued");
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
	tempDomain.ClearC2(rePsi1P, imPsi1P, rePsi2P, imPsi2P, ni*nj*nk);

	if( eError == tempFunction.Copy(inFunction) ) return eError;
	tempFunction.IsFunctionC2(rePhi1P, imPhi1P, rePhi2P, imPhi2P);

	for( s = 0; s < inSteps; s++ ) {
		n = 1 << inFractal >> 1;	//fractal iteration
		for( i = 0; i < n; i++ ) {
			GetFractalNumber( inFractal, i, re, im );				
			re *= inTimeStep;
			im *= inTimeStep;
			Kernel(	rePsi1P,imPsi1P, rePsi2P, imPsi2P, rePhi1P, imPhi1P, rePhi2P, imPhi2P,
					v0P, w0P, a1P, a2P, a3P, domP, re, im, ni, nj, nk );

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

void	TPauli3D::Kernel(	Float* rePsi1P,
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
	Float	reT, imT, reS, imS;
	Float	v0C, w0C;
	Float	a1C, a1R, a1L, a1U, a1D, a1F, a1B;
	Float	a2C, a2R, a2L, a2U, a2D, a2F, a2B;
	Float	a3C, a3R, a3L, a3U, a3D, a3F, a3B;
	Float	chh, ceh, deh, cee, e;
	Float	one = 1.0, two = 2.0, six = 6.0;
	Int32	i, mode = 0;


	chh = one / (two * mMass * mUnits * mUnits);
	ceh = mCharge / (two * mMass * mUnits);
	deh = ceh / two;
	cee = mCharge * mCharge / (two * mMass);
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
			rePsi1C = *rePsi1P++;
			imPsi1C = *imPsi1P++;
			reT = rePsi1R + rePsi1L + rePsi1U + rePsi1D + rePsi1F + rePsi1B;
			imT = imPsi1R + imPsi1L + imPsi1U + imPsi1D + imPsi1F + imPsi1B;
			reEta1C = chh * (six * rePsi1C - reT);
			imEta1C = chh * (six * imPsi1C - imT);
			if( mode & kVector ) {
				a1C = *a1P;
				a2C = *a2P;
				a3C = *a3P;
				reEta1C -= ceh * a1C * (imPsi1R - imPsi1L);
				imEta1C += ceh * a1C * (rePsi1R - rePsi1L);
				reEta1C -= ceh * a2C * (imPsi1U - imPsi1D);
				imEta1C += ceh * a2C * (rePsi1U - rePsi1D);
				reEta1C -= ceh * a3C * (imPsi1F - imPsi1B);
				imEta1C += ceh * a3C * (rePsi1F - rePsi1B);
			}
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
			rePsi2C = *rePsi2P++;
			imPsi2C = *imPsi2P++;
			reT = rePsi2R + rePsi2L + rePsi2U + rePsi2D + rePsi2F + rePsi2B;
			imT = imPsi2R + imPsi2L + imPsi2U + imPsi2D + imPsi2F + imPsi2B;
			reEta2C = chh * (six * rePsi2C - reT);
			imEta2C = chh * (six * imPsi2C - imT);
			if( mode & kVector ) {
				a1C = *a1P;
				a2C = *a2P;
				a3C = *a3P;
				reEta2C -= ceh * a1C * (imPsi2R - imPsi2L);
				imEta2C += ceh * a1C * (rePsi2R - rePsi2L);
				reEta2C -= ceh * a2C * (imPsi2U - imPsi2D);
				imEta2C += ceh * a2C * (rePsi2U - rePsi2D);
				reEta2C -= ceh * a3C * (imPsi2F - imPsi2B);
				imEta2C += ceh * a3C * (rePsi2F - rePsi2B);
			}
			if( mode & kVector ) {
				a1R = *(a1P + 1);
				a2R = *(a2P + 1);
				a3R = *(a3P + 1);
				a1L = *(a1P +-1);
				a2L = *(a2P +-1);
				a3L = *(a3P +-1);
				a1U = *(a1P + ni);
				a2U = *(a2P + ni);
				a3U = *(a3P + ni);
				a1D = *(a1P +-ni);
				a2D = *(a2P +-ni);
				a3D = *(a3P +-ni);
				a1F = *(a1P + ni*nj);
				a2F = *(a2P + ni*nj);
				a3F = *(a3P + ni*nj);
				a1B = *(a1P +-ni*nj);
				a2B = *(a2P +-ni*nj);
				a3B = *(a3P +-ni*nj);
				a1C = *a1P++;
				a2C = *a2P++;
				a3C = *a3P++;
				reT = cee * (a1C * a1C + a2C * a2C + a3C * a3C);
				reEta1C += reT * rePsi1C;
				imEta1C += reT * imPsi1C;
				reEta2C += reT * rePsi2C;
				imEta2C += reT * imPsi2C;
				imT = deh * (a1R - a1L + a2U - a2D + a3F - a3B);
				reEta1C -= imT * imPsi1C;
				imEta1C += imT * rePsi1C;
				reEta2C -= imT * imPsi2C;
				imEta2C += imT * rePsi2C;
				reT = deh * (a1U - a1D - a2R + a2L);
				reEta1C += reT * rePsi1C;
				imEta1C += reT * imPsi1C;
				reEta2C -= reT * rePsi2C;
				imEta2C -= reT * imPsi2C;
				imS = deh * (a3R - a3L - a1F + a1B);
				reEta1C += imS * imPsi2C;
				imEta1C -= imS * rePsi2C;
				reEta2C -= imS * imPsi1C;
				imEta2C += imS * rePsi1C;
				reS = deh * (a2F - a2B - a3U + a3D);
				reEta1C += reS * rePsi2C;
				imEta1C += reS * imPsi2C;
				reEta2C += reS * rePsi1C;
				imEta2C += reS * imPsi1C;
			}
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


