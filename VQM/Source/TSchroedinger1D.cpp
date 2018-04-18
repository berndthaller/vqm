// ===========================================================================
//	TSchroedinger1D.cp
//	Project: Visual Quantum Mechanics.
//  Supervision: Bernd Thaller.
//  Programming: Wolfgang Thaller 1998
// ===========================================================================
//
//	Class for Schroedinger operators

#include "MathLinkUtilities.h"
#include "TSchroedinger1D.h"

//RM (2007/06/28) : minimal change to suppress compiler warnings:
//#include <complex.h>
#include <complex>
using std::complex;


// ---------------------------------------------------------------------------
//		$ TSchroedinger1D
// ---------------------------------------------------------------------------
//	Constructor

TSchroedinger1D::TSchroedinger1D(	Int32 inScalarID,
									Float inMass,
									Float inUnits )
{
	mScalarID = inScalarID;
	mB = 1/(2*inMass);
	mUnits = inUnits;
}



// ---------------------------------------------------------------------------
//		$ ~TSchroedinger1D
// ---------------------------------------------------------------------------
//	Destructor

TSchroedinger1D::~TSchroedinger1D( void )
{

}



// ---------------------------------------------------------------------------
//		$ PutInfo
// ---------------------------------------------------------------------------
//	PutInfo

Int32	TSchroedinger1D::PutInfo( void )
{	
	MLPutFunction(stdlink, "List", 4);
	
		MLPutFunction(stdlink, "Rule", 2);
			MLPutSymbol(stdlink, "Type");	
			MLPutSymbol(stdlink, "Schroedinger1D");

		MLPutFunction(stdlink, "Rule", 2);
			MLPutSymbol(stdlink, "ScalarPotential");
			if(!mScalarID)
				MLPutSymbol(stdlink, "None");
			else
			{
				MLPutFunction(stdlink, "QFunctionObject", 1);
				MLPutInteger(stdlink, mScalarID);
			}
			
		MLPutFunction(stdlink, "Rule", 2);
			MLPutSymbol(stdlink, "Mass");
			MLPutDouble(stdlink, 1/(2*mB));

		MLPutFunction(stdlink, "Rule", 2);
			MLPutSymbol(stdlink, "Units");
			MLPutDouble(stdlink, mUnits);

	return eOK;
}



// ---------------------------------------------------------------------------
//		$ TimeEvolution
// ---------------------------------------------------------------------------
//	TimeEvolution

Int32	TSchroedinger1D::TimeEvolution(	TFunction* inFunction,
										Float inTimeStep,
										Int32 inFractal,
										Int32 inSteps )
{
	Int32	ni;
	Float *rePsiP, *imPsiP;
	Float *reVP, *imVP;			// V at time j
	Float *reVP1, *imVP1;		// V at time j+1
	TFunction	*scalarP;
	
	ni = 0;
	if( !inFunction->IsFunction1D(ni) )
		return MLErrorReport(stdlink, "one-dimensional wavefunction expected");
	if( !inFunction->IsFunctionC(rePsiP, imPsiP) )
		return MLErrorReport(stdlink, "complex wavefunction expected");

	if( mScalarID ) {	//	scalar potential
		scalarP = gFunctionList->Fetch(mScalarID);
		if( !scalarP )
			return MLErrorReport(stdlink, "invalid function ID for scalar potential");
		if( !scalarP->IsFunction1D(ni) )
			return MLErrorReport(stdlink, "scalar potential is not compatible");
		if( !scalarP->IsFunctionC(reVP, imVP) )
			return MLErrorReport(stdlink, "scalar potential is not C-valued");
		reVP1 = reVP;	// for now, there is no change in V
		imVP1 = imVP;	
	}
	
	
	Float reU = mB/(mUnits*mUnits);
	Float imV = 1/inTimeStep;
	
	typedef complex<Float> fcomplex;
	fcomplex *a,*b,*d,*f;
	a = new fcomplex [ni];
	b = new fcomplex [ni];
	d = new fcomplex [ni];
	f = new fcomplex [ni];
	
	int n = ni-1;
	for(Int32 s = 0;s < inSteps; s++)
	{
			// Step 1
		for(int i=1;i<=n-1;i++)
		{
			Float rePOT,imPOT,rePOT1,imPOT1;
			
			if(mScalarID)
			{
				rePOT = reVP[i];
				imPOT = imVP[i];
				rePOT1 = reVP1[i];
				imPOT1 = imVP1[i];
			}
			else
				rePOT = imPOT = rePOT1 = imPOT1 = 0;
		
			fcomplex psiHere = fcomplex(rePsiP[i],imPsiP[i]);
			fcomplex psiRight = fcomplex(rePsiP[i+1],imPsiP[i+1]);
			fcomplex psiLeft = fcomplex(rePsiP[i-1],imPsiP[i-1]);

			a[i] = psiHere*fcomplex(2*reU+rePOT,2*imV+imPOT) -
				reU*(psiRight+psiLeft);
			d[i] = fcomplex(-2*reU-rePOT1,2*imV-imPOT1);
		}
		
			// Step 2
		f[1] = d[1];
		for(int i=2;i<=n-1;i++)
			f[i] = d[i] - (reU*reU)/f[i-1];
			
			// Step 3
		b[0] = 0;
		for(int i=1;i<=n-1;i++)
			b[i] = (a[i]-reU*b[i-1])/f[i];
			
			// Step 4
		
		fcomplex otherPsi;
		fcomplex newPsi;
		
		rePsiP[n] = imPsiP[n] = 0;
		
		for(int i=n-1;i>=1;i--)
		{
			otherPsi = fcomplex(rePsiP[i+1],imPsiP[i+1]);
			
			newPsi = b[i] - reU*otherPsi/f[i];
			
			rePsiP[i] = newPsi.real();
			imPsiP[i] = newPsi.imag();
		}
		
		rePsiP[0] = imPsiP[0] = 0;
	}
	
	delete[] a;
	delete[] b;
	delete[] d;
	delete[] f;
	
	MLPutSymbol(stdlink, "Null");
	return eOK;
}


