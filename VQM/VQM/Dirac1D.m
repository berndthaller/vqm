(* :Title:   Dirac1D *)

(* :Name:    VQM`Dirac1D` *)

(* :Copyright: Copyright 2007 Bernd Thaller *)

(* :Author:  Bernd Thaller,
             Institute of Mathematics,
             University of Graz,
             A-8010 Graz
             bernd.thaller@uni-graz.at
*)

(* :Source:
	Advanced Visual Quantum Mechanics
	Springer-Verlag New York, 2004
*)

(* :Summary:
This package provides code for solving the Dirac equation in
one dimension (a two-by-two matrix-differential equation).
For the free Dirac equation we use a fast Fourier transform.
For the Dirac equation with an external field we use the Crank-Nicolson
formula.
*)

(* :Date:    Jul 20, 2007 *)

(* :Package Version:        2.0 *)

(* :Mathematica Version:    6.0.1 *)

(* :Keywords:
    Dirac equation, Quantum mechanics, Wavefunction
*)

VQMmsgon  = Head[General::"spell"]  =!= $Off;
VQMmsgon1 = Head[General::"spell1"] =!= $Off;

Off[General::spell1,General::spell];

(*-----------------------------------*)
BeginPackage[
	"VQM`Dirac1D`",		(* package Context *)
	"VQM`Spinors`",
	"VQM`FastFourier`"
        ];
(*-----------------------------------*)

ClearAll[$c, $m];

VQM`Dirac1D`Private`Symbols=Hold[
	QRelativisticEnergy, QRelativisticVelocity, QSignedEnergy, QSignedMomentum,
	aplus, aminus,
	QBaseSpinorRight, QBaseSpinorLeft, QBaseSpinorPos, QBaseSpinorNeg,
	QDiracPlaneWaveRight, QDiracPlaneWaveLeft, QDiracPlaneWavePos, QDiracPlaneWaveNeg,
	QDiracMatrix1D, QProjectPositiveEnergy, QProjectNegativeEnergy,
	QFWMatrix, QInverseFWMatrix,
	QDiracEquation1D, QDiracOperator1D,
	QMomentumSpacePropagator,
	QTridiagonalBlockSolve, QDiracTimeStep1D,
	QInitializeDirac1D,
	QMomentumSpaceSpinorPos,QMomentumSpaceSpinorNeg,
	QMomentumSpaceSpinorRight,QMomentumSpaceSpinorLeft,
	QMomentumSpaceSpinorPosRight,QMomentumSpaceSpinorPosLeft,
	QMomentumSpaceSpinorNegRight,QMomentumSpaceSpinorNegLeft,
	QMomentumSpaceSpinor,QPositionSpaceSpinor,
	QPositionSpaceSpinorPos,QPositionSpaceSpinorNeg,
	QPositionSpaceSpinorRight,QPositionSpaceSpinorLeft,
	QPositionSpaceSpinorPosRight,QPositionSpaceSpinorPosLeft,
	QPositionSpaceSpinorNegRight,QPositionSpaceSpinorNegLeft,
	QPositionSpaceGrid,QPositionSpaceInterval,QPositionSpaceStep,
	QComputeMeanPosition, QInnerProductList, QSpinorFT, QInverseSpinorFT];

Unprotect @@ VQM`Dirac1D`Private`Symbols;
ClearAll @@ VQM`Dirac1D`Private`Symbols;


QRelativisticEnergy::usage = "QRelativisticEnergy[k] is the relativistic energy-
	momentum relation. Package: VQM`Dirac1D`.";

QRelativisticVelocity::usage = "QRelativisticVelocity[k] gives the velocity as a
	function of the momentum according to relativistic kinematics. Package: VQM`Dirac1D`.";

QSignedEnergy::usage = "QSignedEnergy[k]=Sign[k]*QRelativisticEnergy[k]. Package: VQM`Dirac1D`.";

QSignedMomentum::usage = "QSignedMomentum[en] relativistic relation between momentum
	and energy. Contains the factor Sign[en]. Package: VQM`Dirac1D`.";

aplus::usage = "aplus[en]. Auxiliary function. Part of the Foldy-Wouthuysen matrix.
	Package: VQM`Dirac1D`.";
aminus::usage = "aminus[en]. Auxiliary function. Part of the Foldy-Wouthuysen matrix.
	Package: VQM`Dirac1D`.";

QDiracMatrix1D::usage = "QDiracMatrix1D[k] is the one-dimensional Dirac operator
	in momentum space. Package: VQM`Dirac1D`.";

QProjectPositiveEnergy::usage = "QProjectPositiveEnergy[k] is the projection onto the
	eigenspace of QDiracMatrix1D[k] belonging to the positive eigenvalue. Package: VQM`Dirac1D`.";

QProjectNegativeEnergy::usage = "QProjectNegativeEnergy[k] is the projection onto the
	eigenspace of QDiracMatrix1D[k] belonging to the negative eigenvalue. Package: VQM`Dirac1D`.";

QMomentumSpacePropagator::usage = "QMomentumSpacePropagator[k,t] is the matrix
	giving the momentum-space representation of the unitary time evolution
	according to the free Dirac equation in one space-dimensions. Package: VQM`Dirac1D`.";

QBaseSpinorRight::usage = "QBaseSpinorRight[k] is the 'spinor-amplitude' of a
	right-moving plane wave of the Dirac equation in one dimension. Package: VQM`Dirac1D`.";

QBaseSpinorLeft::usage = "QBaseSpinorLeft[k] is the 'spinor-amplitude' of a
	left-moving plane wave of the Dirac equation in one dimension. Package: VQM`Dirac1D`.";

QBaseSpinorPos::usage = "QBaseSpinorPos[k] is the 'spinor-amplitude' of a
	positive-energy plane wave of the Dirac equation in one dimension. Package: VQM`Dirac1D`.";

QBaseSpinorNeg::usage = "QBaseSpinorNeg[k] is the 'spinor-amplitude' of a
	negative-energy plane wave of the Dirac equation in one dimension. Package: VQM`Dirac1D`.";

QDiracPlaneWaveRight::usage = "QDiracPlaneWaveRight[k,x] or QDiracPlaneWaveRight[k,x,t] gives
	a plane-wave solution of the one-dimensional Dirac equation, corresponding
	to a beam of particles with positive velocity. Package: VQM`Dirac1D`.";

QDiracPlaneWaveLeft::usage = "QDiracPlaneWaveLeft[k,x] or QDiracPlaneWaveLeft[k,x,t] gives
	a plane-wave solution of the one-dimensional Dirac equation, corresponding
	to a beam of particles with negative velocity. Package: VQM`Dirac1D`.";

QDiracPlaneWavePos::usage = "QDiracPlaneWavePos[k,x] or QDiracPlaneWavePos[k,x,t] gives
	a plane-wave solution of the one-dimensional Dirac equation, corresponding
	to particles with momentum k and positive energy. Package: VQM`Dirac1D`.";

QDiracPlaneWaveNeg::usage = "QDiracPlaneWaveNeg[k,x] or QDiracPlaneWaveNeg[k,x,t] gives
	a plane-wave solution of the one-dimensional Dirac equation, corresponding
	to particles with momentum k and negative energy. Package: VQM`Dirac1D`.";

QFWMatrix::usage = "QFWMatrix[k] and QInverseFWMatrix[k] is the unitary matrix 
	(resp. its inverse) transforming QDiracMatrix1D[k] to a diagonal form.
	Package: VQM`Dirac1D`.";

QInverseFWMatrix::usage = "QFWMatrix[k] and QInverseFWMatrix[k] is the unitary matrix 
	(resp. its inverse) transforming QDiracMatrix1D[k] to a diagonal form.
	Package: VQM`Dirac1D`.";

$c::usage = "Constant representing the velocity of light. Package: VQM`Dirac1D`.";

$m::usage = "Constant representing the mass of the particle. Package: VQM`Dirac1D`.";

QDiracEquation1D::usage = "QDiracEquation1D[psi[x,t],{x,t}]==0 is the
	one-dimensional Dirac equation for the two-component spinor function
	psi[x,t]. Package: VQM`Dirac1D`.";

QDiracOperator1D::usage = "QDiracOperator1D[psi[x],x]==0 gives the action
	one-dimensional Dirac operator on the two-component spinor function
	psi[x]. Package: VQM`Dirac1D`.";

QTridiagonalBlockSolve::usage = "QTridiagonalBlockSolve[a_List, b_List, c_List, r_List]
	is a variant of TridiagonalSolve in the standard package
	LinearAlgebra`Tridiagonal`. a,b,c are lists of 2 by 2 matrices, Length[b]=
	Length[c]+1=Length[a]+1, and r is
	a list of two-dimensional vectors, Length[r]=Length[b].
	If M[a,b,c] is the block-tridiagonal matrix with b as the diagonal, then
	QTridiagonalBlockSolve gives the solution of M.x=r as a list
	x of two-dimensional vectors. Package: VQM`Dirac1D`.";

QDiracTimeStep1D::usage = "QDiracTimeStep1D[psi_List, pot_List, dt, dx]
	performs a time-step dt of the time-evolution defined by the free,
	one-dimensional Dirac equation with potential-matrix pot.
	psi is given as a list of two-dimensional vectors that represent
	the values of the initial spinor on the points of a regular grid
	with grid-distance dx. pot is a list of 2 by 2 matrices representing
	the values of a potential matrix on the grid points. The algorithm
	uses a Crank-Nicolson formula. Package: VQM`Dirac1D`.";

QInitializeDirac1D::usage = "QInitializeDirac1D[Kpoints] must be called before executing
	any of the QMomentumSpaceSpinor- or QPositionSpaceSpinor- commands. It produces
	some lists (global variables) which are needed by all of these commands and which
	are pre-calculated once and for all for better performance.
	Kpoints is the grid of points in momentum space, where the initial spinor is
	given as a list of spinor-values. Package: VQM`Dirac1D`.";

QMomentumSpaceSpinorPos::usage = "QMomentumSpaceSpinorPos[t,Kvalues,Kpoints] computes
	the positive-energy part of QMomentumSpaceSpinor[t,Kvalues,Kpoints]. Package: VQM`Dirac1D`.";

QMomentumSpaceSpinorNeg::usage = "QMomentumSpaceSpinorNeg[t,Kvalues,Kpoints] computes
	the negative-energy part of QMomentumSpaceSpinor[t,Kvalues,Kpoints]. Package: VQM`Dirac1D`.";

QMomentumSpaceSpinorRight::usage = "QMomentumSpaceSpinorRight[t,Kvalues,Kpoints] computes
	the part with positive velocity of QMomentumSpaceSpinor[t,Kvalues,Kpoints]. Package: VQM`Dirac1D`.";

QMomentumSpaceSpinorLeft::usage = "QMomentumSpaceSpinorLeft[t,Kvalues,Kpoints] computes
	the part with negative velocity of QMomentumSpaceSpinor[t,Kvalues,Kpoints]. Package: VQM`Dirac1D`.";
	
QMomentumSpaceSpinorPosRight::usage = "QMomentumSpaceSpinorPosRight[t,Kvalues,Kpoints] computes
	the part with positive energy and positive velocity of
	QMomentumSpaceSpinor[t,Kvalues,Kpoints]. Package: VQM`Dirac1D`.";

QMomentumSpaceSpinorPosLeft::usage = "QMomentumSpaceSpinorPosLeft[t,Kvalues,Kpoints] computes
	the part with positive energy and negative velocity of
	QMomentumSpaceSpinor[t,Kvalues,Kpoints]. Package: VQM`Dirac1D`.";

QMomentumSpaceSpinorNegRight::usage = "QMomentumSpaceSpinorNegRight[t,Kvalues,Kpoints] computes
	the part with negative energy and positive velocity of
	QMomentumSpaceSpinor[t,Kvalues,Kpoints]. Package: VQM`Dirac1D`.";

QMomentumSpaceSpinorNegLeft::usage = "QMomentumSpaceSpinorNegLeft[t,Kvalues,Kpoints] computes
	the part with negative energy and negative velocity of
	QMomentumSpaceSpinor[t,Kvalues,Kpoints]. Package: VQM`Dirac1D`.";

QMomentumSpaceSpinor::usage = "QMomentumSpaceSpinor[t,Kvalues,Kpoints] computes
	the solution of the one-dimensional Dirac equation in momentum space at time t,
	given an initial spinor Kvalues on a grid Kpoints in momentum space. The result
	is given as a list of spinor values (that is, a list of two-dimensional complex vectors).
	Package: VQM`Dirac1D`.";

QPositionSpaceSpinor::usage = "QPositionSpaceSpinor[t,Kvalues,Kpoints] gives the solution
	of the one-dimensional Dirac equation in position space at time t, given an initial
	spinor Kvalues on a grid Kpoints in momentum space. The result
	is given as a list of spinor values (that is, a list of two-dimensional complex
	vectors). These values give the solution on the set of points described by
	QPositionSpaceGrid[Kpoints]. Package: VQM`Dirac1D`.";

QPositionSpaceSpinorPos::usage = "QPositionSpaceSpinorPos[t,Kvalues,Kpoints] computes
	the positive-energy part of QPositionSpaceSpinor[t,Kvalues,Kpoints]. Package: VQM`Dirac1D`.";

QPositionSpaceSpinorNeg::usage = "QPositionSpaceSpinorNeg[t,Kvalues,Kpoints] computes
	the negative-energy part of QPositionSpaceSpinor[t,Kvalues,Kpoints]. Package: VQM`Dirac1D`.";

QPositionSpaceSpinorRight::usage = "QPositionSpaceSpinorRight[t,Kvalues,Kpoints] computes
	the part with positive velocity of QPositionSpaceSpinorPos[t,Kvalues,Kpoints].
	Package: VQM`Dirac1D`.";

QPositionSpaceSpinorLeft::usage = "QPositionSpaceSpinorLeft[t,Kvalues,Kpoints] computes
	the part with negative velocity of QPositionSpaceSpinor[t,Kvalues,Kpoints].
	Package: VQM`Dirac1D`.";

QPositionSpaceSpinorPosRight::usage = "QPositionSpaceSpinorPosRight[t,Kvalues,Kpoints] computes
	the part with positive energy and positive velocity of QPositionSpaceSpinor[t,Kvalues,Kpoints].
	Package: VQM`Dirac1D`.";

QPositionSpaceSpinorPosLeft::usage = "QPositionSpaceSpinorPosLeft[t,Kvalues,Kpoints] computes
	the part with positive energy and negative velocity of QPositionSpaceSpinor[t,Kvalues,Kpoints].
	Package: VQM`Dirac1D`.";

QPositionSpaceSpinorNegRight::usage = "QPositionSpaceSpinorNegRight[t,Kvalues,Kpoints] computes
	the part with negative energy and positive velocity of QPositionSpaceSpinor[t,Kvalues,Kpoints].
	Package: VQM`Dirac1D`.";

QPositionSpaceSpinorNegLeft::usage = "QPositionSpaceSpinorNegLeft[t,Kvalues,Kpoints] computes
	the part with negative energy and negative velocity of QPositionSpaceSpinor[t,Kvalues,Kpoints].
	Package: VQM`Dirac1D`.";
	
QPositionSpaceGrid::usage = "QPositionSpaceGrid[Kpoints] gives the grid of points in
	position space, where the solution is determined by QPositionSpaceSpinor.
	Package: VQM`Dirac1D`.";

QPositionSpaceInterval::usage = "QPositionSpaceInterval[Kpoints] gives the interval in
	position space, where the solution is determined by QPositionSpaceSpinor.
	Package: VQM`Dirac1D`.";

QPositionSpaceStep::usage = "QPositionSpaceStep[Kpoints] gives the step size of the grid in
	position space, where the solution is determined by QPositionSpaceSpinor.
	Package: VQM`Dirac1D`.";

QComputeMeanPosition::usage = "QComputeMeanPosition[Xvalues,Xpoints] computes
	the mean position of a spinor that is given by a list of spinor-values on points
	Xpoints of the one-dimensional position space.
	Package: VQM`Dirac1D`.";

QInnerProductList::usage = "QInnerProductList[spinorlist1,spinorlist2]";

QSpinorFT::usage = "QSpinorFT[spinorlist,{x1,x2}] computes the component-wise Fourier
	transform of the spinorlist, assuming that spinorlist are C^2-values given on a
	uniform grid of points in the interval {x1,x2}. For example, spinorlist could be
	the numerical approximation of a spinor-wave function. The Fourier-transform is
	computed via a fast Fourier transform method provided by the package FastFourier.m.
	Package: VQM`Dirac1D`.";

QInverseSpinorFT::usage = "QInverseSpinorFT[spinorlist,{k1,k2}] computes
	the component-wise inverse Fourier
	transform of the spinorlist, assuming that spinorlist are C^2-values given on a
	uniform grid of points in the interval {k1,k2}. For example, spinorlist could be
	the numerical approximation of a spinor-wave function in momentum space.
	The inverse Fourier-transform is
	computed via a fast Fourier transform method provided by the package FastFourier.m.
	Package: VQM`Dirac1D`.";

(*-----------------------------------*)
Begin["`Private`"];
(*-----------------------------------*)

(* mysign[k_]:=Which[k>=0, 1, k<0,-1]; *)
(* RM: not a big difference, but neater *)
mysign[k_?NumericQ] := Piecewise[{ {1, k >= 0.}, {-1, k < 0.} }];

QRelativisticEnergy[k_] := Sqrt[$c^2 k^2 + $m^2 $c^4];

QRelativisticVelocity[k_] := $c^2 k/QRelativisticEnergy[k];

QSignedEnergy[k_] := mysign[k] QRelativisticEnergy[k];

QSignedMomentum[en_] := mysign[en] 1/$c Sqrt[en^2 -$m^2 $c^4];

aplus[en_]:= Sqrt[(1 + $m $c^2 / en) /2];

aminus[en_]:= Sqrt[(1 - $m $c^2 / en) /2];

ap[k_]:=aplus[QRelativisticEnergy[k]];

am[k_]:=aminus[QRelativisticEnergy[k]];

QFWMatrix[k_]:= ap[k] QIdentity2 + am[k] mysign[k] I QPauliSigma2;
QInverseFWMatrix[k_] := ap[k] QIdentity2 - am[k] mysign[k] I QPauliSigma2;
	(* == Inverse[QFWMatrix[k]] *)

QBaseSpinorRight[k_] := {aplus[QSignedEnergy[k]], aminus[QSignedEnergy[k]]};
QBaseSpinorLeft[k_]  := mysign[k] {-aminus[QSignedEnergy[k]], aplus[QSignedEnergy[k]]};
QBaseSpinorPos[k_]   := {ap[k], mysign[k] am[k]};
QBaseSpinorNeg[k_]  := {-mysign[k] am[k], ap[k]};

QDiracPlaneWaveRight[k_,x_] := 1/Sqrt[2 Pi] QBaseSpinorRight[k] Exp[I k x];
QDiracPlaneWaveRight[k_,x_,t_] := QDiracPlaneWaveRight[k,x] Exp[- I QSignedEnergy[k] t];

QDiracPlaneWaveLeft[k_,x_] := 1/Sqrt[2 Pi] QBaseSpinorLeft[k] Exp[I k x];
QDiracPlaneWaveLeft[k_,x_,t_] := QDiracPlaneWaveLeft[k,x] Exp[I QSignedEnergy[k] t];
	
QDiracPlaneWavePos[k_,x_] := 1/Sqrt[2 Pi] QBaseSpinorPos[k] Exp[I k x];
QDiracPlaneWavePos[k_,x_,t_] := QDiracPlaneWavePos[k,x] Exp[- I QRelativisticEnergy[k] t];

QDiracPlaneWaveNeg[k_,x_] := 1/Sqrt[2 Pi] QBaseSpinorNeg[k] Exp[I k x];
QDiracPlaneWaveNeg[k_,x_,t_] := QDiracPlaneWaveNeg[k,x] Exp[I QRelativisticEnergy[k] t];
	
OmegaRight[en_,x_] :=
	1/Sqrt[2 Pi] Sqrt[ en/($c^2 QSignedMomentum[en]) ]*
	Exp[ I QSignedMomentum[en] x] {aplus[en], aminus[en]}

OmegaLeft[en_,x_] :=
	mysign[en]/Sqrt[2 Pi] Sqrt[ en/($c^2 QSignedMomentum[en]) ]*
	Exp[I QSignedMomentum[en] x] {-aminus[en], aplus[en]}

OmegaRight[en_,x_,t_] :=
	QDiracPlaneWaveRight[en,x] Exp[-I en t]

OmegaLeft[en_,x_,t_] :=
	QDiracPlaneWaveLeft[en,x] Exp[I en t]

QDiracMatrix1D[k_] := QPauliSigma1 $c k + QPauliSigma3 $m $c^2

QProjectPositiveEnergy[k_] :=
	1/2 (QIdentity2  +  QDiracMatrix1D[k]/QRelativisticEnergy[k])

QProjectNegativeEnergy[k_] :=
	1/2 (QIdentity2  -  QDiracMatrix1D[k]/QRelativisticEnergy[k])

(* Note:
	QProjectPositiveEnergy[k].QBaseSpinorPos[k] == QBaseSpinorPos[k];
	QProjectNegativeEnergy[k].QBaseSpinorPos[k] == 0;
*)

QMomentumSpacePropagator[k_,t_] :=
	{
		{Cos[t*Sqrt[$c^2*(k^2 + $c^2*$m^2)]] -
			(I*$c^2*$m*Sin[t*Sqrt[$c^2*(k^2 + $c^2*$m^2)]])/\
				Sqrt[$c^2*(k^2 + $c^2*$m^2)], 
		(-I*k*$c*Sin[t*Sqrt[$c^2*(k^2 + $c^2*$m^2)]])/\
			Sqrt[$c^2*(k^2 + $c^2*$m^2)]
		}, 
 		{(-I*k*$c*Sin[t*Sqrt[$c^2*(k^2 + $c^2*$m^2)]])/\
 			Sqrt[$c^2*(k^2 + $c^2*$m^2)],
 		Cos[t*Sqrt[$c^2*(k^2 + $c^2*$m^2)]] +
 			(I*$c^2*$m*Sin[t*Sqrt[$c^2*(k^2 + $c^2*$m^2)]])/\
 				Sqrt[$c^2*(k^2 + $c^2*$m^2)]
 		}
 	};

(* QMomentumSpacePropagator is the same as
	Exp[ I  QRelativisticEnergy[k] t] QProjectNegativeEnergy[k] +
	Exp[- I QRelativisticEnergy[k] t] QProjectPositiveEnergy[k]
*)

QDiracEquation1D[spinorfunc_,{x_Symbol,t_Symbol}] := 
Module[{sp = Function[{x,t},spinorfunc]},I D[sp[x,t],t] +
	I $c QPauliSigma1.D[sp[x,t],x] -  $m $c^2 QPauliSigma3.sp[x,t]]

QDiracOperator1D[spinorfunc_,x_Symbol] := 
Module[{sp = Function[{x},spinorfunc]}, -
	I $c QPauliSigma1.D[sp[x],x] + $m $c^2 QPauliSigma3.sp[x]]

(* numerical solution *)

(* QTridiagonalBlockSolve is a variant of TridiagonalSolve in the
standard package LinearAlgebra`Tridiagonal`. *)

QTridiagonalBlockSolve[a_List, b_List, c_List, r_List]:=
    Module[{len=Length[r],
			solution=Array[{0,0}, Length[r]], aux, 
			aux1=Array[{{0,0},{0,0}}, Length[r]],
			a1 = Prepend[a, {{0,0},{0,0}}],
			iter}, 
		aux = Inverse[b[[1]]];
		solution[[1]] = aux.r[[1]] ;
		Do[aux1[[iter]] = aux.c[[iter-1]];
	    	aux = Inverse[b[[iter]]-a1[[iter]].aux1[[iter]]];
	    	solution[[iter]] =  aux.(r[[iter]]-a1[[iter]].solution[[iter-1]]),
	    	{iter, 2, len}]; 			
		Do[solution[[iter]] -= aux1[[iter+1]].solution[[iter+1]],
	    	{iter, len-1, 1, -1}];
		solution
	] /; Length[a] + 1 == Length[b] == Length[c] + 1 == Length[r]

QDiracTimeStep1D[psi_List, pot_List, dt_, dx_] :=
	Module[{v = (2 I/dt)*QIdentity2,
			u = -$c I/(2 dx)*QPauliSigma1,
			dim = Length[psi],
			p = Plus[$m $c^2 QPauliSigma3,#]& /@ pot,
			a,b,c,r},
		a = Table[u,{i,dim-1}]; c = -a;
		c[[1]] = -2 u; a[[dim-1]] = 2 u;
		b = Table[v,{i,dim}] - p;
		b[[1]] += 2 u; b[[dim]] -= 2 u;
		r = Table[-u.psi[[i-1]] + (v + p[[i]]).psi[[i]] + u.psi[[i+1]],
			 	{i,2,dim-1}];
		r = Prepend[r, (v - 2 u + p[[1]]).psi[[1]] + 2 u.psi[[2]] ];
		r = Append[r, (v + 2 u + p[[dim]]).psi[[dim]] - 2 u.psi[[dim-1]] ];
		QTridiagonalBlockSolve[a, b, c, r]
	] /; Length[psi] == Length[pot]

(* numerical solution via fast Fourier transform *)

(*	generate some pre-calculated lists for faster performance.
	Stored as global variables. QInitializeDirac1D must be
	executed before evaluating the spinors in position space *)
(*	It is assumed that the values of a spinor are given on a
	grid of points 'Kpoints' in momentum space *)
   
QInitializeDirac1D[Kpoints_]:=
	Module[{},
		baseSpinorLeftList = QBaseSpinorLeft[#]& /@ Kpoints //N;
		baseSpinorRightList = QBaseSpinorRight[#]& /@ Kpoints //N;
		baseSpinorPosList = QBaseSpinorPos[#]& /@ Kpoints //N;
		baseSpinorNegList = QBaseSpinorNeg[#]& /@ Kpoints //N;
	]

(* auxiliary quantities: *)

(* inner product of two spinors given as lists of elements {z1,z2} *)

QInnerProductList[spinorlist1_,spinorlist2_] :=
	Conjugate[First /@ spinorlist1] * First /@ spinorlist2 +
	Conjugate[Last /@ spinorlist1] * Last /@ spinorlist2  /;
		Dimensions[spinorlist1]==Dimensions[spinorlist2]

(*	compute the inverse Fourier transform of a spinor that is
	given as a list of values Kvalues on the points Kpoints in momentum space *)

QSpinorFT[Xvalues_,Xpoints_] := Transpose[
	{QFourierList[ First /@ Xvalues, {First[Xpoints],Last[Xpoints]}][[1]],
	 QFourierList[ Last  /@ Xvalues, {First[Xpoints],Last[Xpoints]}][[1]]}]

QInverseSpinorFT[Kvalues_,Kpoints_] := Transpose[
	{QInverseFourierList[ First /@ Kvalues, {First[Kpoints],Last[Kpoints]}][[1]],
	 QInverseFourierList[ Last  /@ Kvalues, {First[Kpoints],Last[Kpoints]}][[1]]}]

(*	projections onto positive and negative momenta in momentum space.
	The spinorlist Kvalues are the values of a spinor in momentum space, given on
	a grid of points Kpoints *)

projectPosMom[Kvalues_,Kpoints_] :=
	PadLeft[Extract[Kvalues,Position[N[Kpoints], _?NonNegative]],Length[Kpoints],{{0,0}}]
projectNegMom[Kvalues_,Kpoints_] :=
	PadRight[Extract[Kvalues,Position[N[Kpoints], _?Negative]],Length[Kpoints],{{0,0}}]
	
(*	some auxiliary quantities *)

timeFactorsPos[t_,Kpoints_] := Exp[-I QRelativisticEnergy[#] t]& /@ Kpoints
timeFactorsNeg[t_,Kpoints_] := Exp[ I QRelativisticEnergy[#] t]& /@ Kpoints
timeFactorsSigned[t_,Kpoints_] := Exp[-I QSignedEnergy[#] t]& /@ Kpoints
timeFactorsSignedConj[t_,Kpoints_] := Exp[I QSignedEnergy[#] t]& /@ Kpoints

(*	Spinors in momentum space *)
(*  QInitializeDirac1D[Kpoints] has to be called first!! *)
(*  It is assumed that a spinor is given as Kvalues on points Kpoints in momentum space *)
(*	The following commands compute various parts of this spinor *)

QMomentumSpaceSpinorPos[t_,Kvalues_,Kpoints_] :=
	baseSpinorPosList * QInnerProductList[baseSpinorPosList,Kvalues] *
	timeFactorsPos[t,Kpoints]

QMomentumSpaceSpinorNeg[t_,Kvalues_,Kpoints_] :=
	baseSpinorNegList * QInnerProductList[baseSpinorNegList,Kvalues] *
	timeFactorsNeg[t,Kpoints]

QMomentumSpaceSpinorRight[t_,Kvalues_,Kpoints_] :=
	baseSpinorRightList * QInnerProductList[baseSpinorRightList,Kvalues] *
	timeFactorsSigned[t,Kpoints]

QMomentumSpaceSpinorLeft[t_,Kvalues_,Kpoints_] :=
	baseSpinorLeftList * QInnerProductList[baseSpinorLeftList,Kvalues] *
	timeFactorsSignedConj[t,Kpoints]

QMomentumSpaceSpinorPosRight[t_,Kvalues_,Kpoints_] :=
	projectPosMom[ QMomentumSpaceSpinorPos[t,Kvalues,Kpoints],Kpoints ]

QMomentumSpaceSpinorPosLeft[t_,Kvalues_,Kpoints_] :=
	projectNegMom[ QMomentumSpaceSpinorPos[t,Kvalues,Kpoints],Kpoints ]

QMomentumSpaceSpinorNegRight[t_,Kvalues_,Kpoints_] :=
	projectNegMom[ QMomentumSpaceSpinorNeg[t,Kvalues,Kpoints],Kpoints ]

QMomentumSpaceSpinorNegLeft[t_,Kvalues_,Kpoints_] :=
	projectPosMom[ QMomentumSpaceSpinorNeg[t,Kvalues,Kpoints],Kpoints ]

QMomentumSpaceSpinor[t_,Kvalues_,Kpoints_] :=
	QMomentumSpaceSpinorPos[t,Kvalues,Kpoints] + QMomentumSpaceSpinorNeg[t,Kvalues,Kpoints]

(* spinors in position space *)

QPositionSpaceSpinorPos[t_,Kvalues_,Kpoints_] :=
	QInverseSpinorFT[ QMomentumSpaceSpinorPos[t,Kvalues,Kpoints],Kpoints ];

QPositionSpaceSpinorNeg[t_,Kvalues_,Kpoints_] :=
	QInverseSpinorFT[ QMomentumSpaceSpinorNeg[t,Kvalues,Kpoints],Kpoints ];

QPositionSpaceSpinor[t_,Kvalues_,Kpoints_] :=
	QInverseSpinorFT[ QMomentumSpaceSpinor[t,Kvalues,Kpoints],Kpoints];

QPositionSpaceSpinorRight[t_,Kvalues_,Kpoints_] :=
	QInverseSpinorFT[ QMomentumSpaceSpinorRight[t,Kvalues,Kpoints],Kpoints];

QPositionSpaceSpinorLeft[t_,Kvalues_,Kpoints_] :=
	QInverseSpinorFT[ QMomentumSpaceSpinorLeft[t,Kvalues,Kpoints],Kpoints];

QPositionSpaceSpinorPosRight[t_,Kvalues_,Kpoints_] :=
	QInverseSpinorFT[ QMomentumSpaceSpinorPosRight[t,Kvalues,Kpoints],Kpoints ];

QPositionSpaceSpinorPosLeft[t_,Kvalues_,Kpoints_] :=
	QInverseSpinorFT[ QMomentumSpaceSpinorPosLeft[t,Kvalues,Kpoints],Kpoints ];

QPositionSpaceSpinorNegRight[t_,Kvalues_,Kpoints_] :=
	QInverseSpinorFT[ QMomentumSpaceSpinorNegRight[t,Kvalues,Kpoints],Kpoints ];

QPositionSpaceSpinorNegLeft[t_,Kvalues_,Kpoints_] :=
	QInverseSpinorFT[ QMomentumSpaceSpinorNegLeft[t,Kvalues,Kpoints],Kpoints ];

(* auxiliary quantities in position space *)

QPositionSpaceGrid[Kpoints_]:=
	QFourierGrid[First[Kpoints],Last[Kpoints],Length[Kpoints]]

QPositionSpaceInterval[Kpoints_]:=
	QFourierInterval[First[Kpoints],Last[Kpoints],Length[Kpoints]]

QPositionSpaceStep[Kpoints_]:=
	QFourierStep[First[Kpoints],Last[Kpoints],Length[Kpoints]]

(* Average position of wave packet *)

QComputeMeanPosition[Xvalues_,Xpoints_] :=
	((Abs[First/@ Xvalues]^2 + Abs[Last /@ Xvalues]^2).Xpoints) * (Xpoints[[2]]-Xpoints[[1]]);


(*-----------------------------------*)
End[];      (* end `Private` Context *)
(*-----------------------------------*)

Protect @@ VQM`Dirac1D`Private`Symbols;

$c = 1; $m = 1;

(*-----------------------------------*)
EndPackage[]; (* end package Context *)
(*-----------------------------------*)

If[VQMmsgon,  On[General::"spell" ]];
If[VQMmsgon1, On[General::"spell1"]];
