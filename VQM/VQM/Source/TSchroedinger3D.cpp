// ===========================================================================
//	TSchroedinger3D.cp
//	Project: Visual Quantum Mechanics.
//  Supervision: Bernd Thaller.
//  Programming: Manfred Liebmann 1998
// ===========================================================================
//
//	Class for Schroedinger operators

#include "MathLinkUtilities.h"
#include "TSchroedinger3D.h"
#include "TDomain.h"



// ---------------------------------------------------------------------------
//		$ TSchroedinger3D
// ---------------------------------------------------------------------------
//	Constructor

TSchroedinger3D::TSchroedinger3D(	Int32 inScalarID,
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
//		$ ~TSchroedinger3D
// ---------------------------------------------------------------------------
//	Destructor

TSchroedinger3D::~TSchroedinger3D( void )
{

}



// ---------------------------------------------------------------------------
//		$ PutInfo
// ---------------------------------------------------------------------------
//	PutInfo

Int32	TSchroedinger3D::PutInfo( void )
{	
	MLPutFunction(stdlink, "List", 7);
	MLPutFunction(stdlink, "Rule", 2);
	MLPutSymbol(stdlink, "Type");	
	MLPutSymbol(stdlink, "Schroedinger3D");

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

Int32	TSchroedinger3D::TimeEvolution(	TFunction* inFunction,
										Float inTimeStep,
										Int32 inFractal,
										Int32 inSteps )
{
	Int32	ni, nj, nk, i, n, s;
	Float	*rePsiP, *imPsiP, *rePhiP, *imPhiP;
	Float	*v0P, *w0P, *a1P, *a2P, *a3P, *d0P, *tempP, re, im;
	TFunction	*scalarP, *vectorP, *domainP, tempFunction;
	TDomain	tempDomain;
	Int8	*domP;


	if( inFractal < 0 || inFractal >= 16 )
		return MLErrorReport(stdlink, "fractal order is out of range");

	ni = nj = nk = 0;
	if( !inFunction->IsFunction3D(ni, nj, nk) )
		return MLErrorReport(stdlink, "three-dimensional wavefunction expected");
	if( !inFunction->IsFunctionC(rePsiP, imPsiP) )
		return MLErrorReport(stdlink, "complex wavefunction expected");

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
	tempDomain.ClearC(rePsiP, imPsiP, ni*nj*nk);

	if( eError == tempFunction.Copy(inFunction) ) return eError;
	tempFunction.IsFunctionC(rePhiP, imPhiP);

	for( s = 0; s < inSteps; s++ ) {
		n = 1 << inFractal >> 1;	//fractal iteration
		for( i = 0; i < n; i++ ) {
			GetFractalNumber( inFractal, i, re, im );				
			re *= inTimeStep;
			im *= inTimeStep;
			Kernel(	rePsiP, imPsiP, rePhiP, imPhiP,
					v0P, w0P, a1P, a2P, a3P, domP, re, im, ni, nj, nk);

			tempP = rePsiP; rePsiP = rePhiP; rePhiP = tempP;
			tempP = imPsiP; imPsiP = imPhiP; imPhiP = tempP;
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

void	TSchroedinger3D::Kernel(	Float* rePsiP,
									Float* imPsiP,
									Float* rePhiP,
									Float* imPhiP,
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
	Float	rePsiC, imPsiC, reEtaC, imEtaC;
	Float	rePsiR, imPsiR, rePsiL, imPsiL;
	Float	rePsiU, imPsiU, rePsiD, imPsiD;
	Float	rePsiF, imPsiF, rePsiB, imPsiB;
	Float	reT, imT;
	Float	v0C, w0C;
	Float	a1C, a1R, a1L;
	Float	a2C, a2U, a2D;
	Float	a3C, a3F, a3B;
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
			rePsiR = *(rePsiP + 1);
			imPsiR = *(imPsiP + 1);
			rePsiL = *(rePsiP +-1);
			imPsiL = *(imPsiP +-1);
			rePsiU = *(rePsiP + ni);
			imPsiU = *(imPsiP + ni);
			rePsiD = *(rePsiP +-ni);
			imPsiD = *(imPsiP +-ni);
			rePsiF = *(rePsiP + ni*nj);
			imPsiF = *(imPsiP + ni*nj);
			rePsiB = *(rePsiP +-ni*nj);
			imPsiB = *(imPsiP +-ni*nj);
			rePsiC = *rePsiP++;
			imPsiC = *imPsiP++;
			reT = rePsiR + rePsiL + rePsiU + rePsiD + rePsiF + rePsiB;
			imT = imPsiR + imPsiL + imPsiU + imPsiD + imPsiF + imPsiB;
			reEtaC = chh * (six * rePsiC - reT);
			imEtaC = chh * (six * imPsiC - imT);
			if( mode & kVector ) {
				a1R = *(a1P + 1);
				a1L = *(a1P +-1);
				a2U = *(a2P + ni);
				a2D = *(a2P +-ni);
				a3F = *(a3P + ni*nj);
				a3B = *(a3P +-ni*nj);
				a1C = *a1P++;
				a2C = *a2P++;
				a3C = *a3P++;
				reEtaC -= ceh * a1C * (imPsiR - imPsiL);
				imEtaC += ceh * a1C * (rePsiR - rePsiL);
				reEtaC -= ceh * a2C * (imPsiU - imPsiD);
				imEtaC += ceh * a2C * (rePsiU - rePsiD);
				reEtaC -= ceh * a3C * (imPsiF - imPsiB);
				imEtaC += ceh * a3C * (rePsiF - rePsiB);
				reT = cee * (a1C * a1C + a2C * a2C + a3C * a3C);
				reEtaC += reT * rePsiC;
				imEtaC += reT * imPsiC;
				imT = deh * (a1R - a1L + a2U - a2D + a3F - a3B);					
				reEtaC -= imT * imPsiC;
				imEtaC += imT * rePsiC;
			}
			if( mode & kScalar ) {
				v0C = *v0P++;
				w0C = *w0P++;
				reEtaC += e * v0C * rePsiC;
				imEtaC += e * v0C * imPsiC;
				reEtaC -= e * w0C * imPsiC;
				imEtaC += e * w0C * rePsiC;
			}
			*rePhiP++ = rePsiC + imZ * reEtaC + reZ * imEtaC;
			*imPhiP++ = imPsiC - reZ * reEtaC + imZ * imEtaC;
		}
		else {			
			rePsiP++;
			imPsiP++;
			if( mode & kScalar ) { v0P++; *w0P++; }
			rePhiP++;
			imPhiP++;
			if( mode & kVector ) { a1P++; a2P++; a3P++; }
		}
	}

}


