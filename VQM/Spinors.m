(* ::Package:: *)

(* :Title:   Spinors *)

(* :Name:    VQM`Spinors` *)

(* :Copyright: Copyright 2004 Bernd Thaller *)

(* :Author:  Bernd Thaller,
             Institute of Mathematics,
             University of Graz,
             A-8010 Graz
             bernd.thaller@uni-graz.at
*)

(* :Source:
	Advanced Visual Quantum Mechanics
	Springer-Verlag New York, 2004
	in particular Chapter 3 and 4
*)

(* :Summary:
VQM`Spinors` is a package for 'Visual Quantum Mechanics'.
It defines basic operations with spinors and provides
tools for the visualization of spinors.
*)

(* :Date:    Jul-20-2007 *)
(* :Date:    Apr-25-2018 *)

(* :Package Version:        3.0 *)

(* :Mathematica Version:    11.3.0 *)

(* :Keywords:
    Spinors, Visualization, SU2, Spinor Harmonics
*)

(* :History:
   change opts___ -> opts___?OptionQ
   QSetSpinBasis : ( )  instead of Module[{}, ]
*)

VQMmsgon  = Head[General::"spell"]  =!= $Off;
VQMmsgon1 = Head[General::"spell1"] =!= $Off;

Off[General::spell1,General::spell];

(*-----------------------------------*)
BeginPackage[
	"VQM`Spinors`",
	"VQM`VisualizeVector`",
	"VQM`ColorMaps`"
	];
(*-----------------------------------*)

VQM`Spinors`Private`Symbols=Hold[
	QNorm, QNormalize, $QSpinBasis, QxBasis, QyBasis, QzBasis, QSetSpinBasis, QSpinBasis, QUseBasis,
	QSpinorToComponents, QComponentsToSpinor, QProjectUp, QProjectDown, QProjection,
	QProbabilityUp, QProbabilityDown, QSpinorBarDiagram,
	QPauliSigma1, QPauliSigma2, QPauliSigma3, QPauliSigmaV, QIdentity2,
	QVectorToHermitianMatrix, QVectorToDensityMatrix, QSpinHamiltonian,
	QSpinorUp, QSpinorDown, QConjSpinorUp, QConjSpinorDown,
	QSpinorHarmonicUp, QSpinorHarmonicDown, QSpinorToVector, QVectorToSpinor,
	QVectorLength, QExtractPhase, QHermitianMatrixToVector, QDensityMatrixToVector,
	QSpinorToArrow, QRotationSO3, QRotationSU2,
	QVisualizeSpinor, QVisualizeDensityMatrix];


Unprotect @@ VQM`Spinors`Private`Symbols;
ClearAll @@ VQM`Spinors`Private`Symbols;


QNorm::usage = "QNorm[spinor] determines the norm of a (complex) vector. Package: VQM`Spinors`.";

QNormalize::usage = "QNormalize[spinor] returns a (complex) vector with norm 1. Package: VQM`Spinors`.";

$QSpinBasis::usage = "Gives the standard reference basis in the two-dimensional complex linear space.
	Package: VQM`Spinors`.";

QxBasis::usage = "Basis in the two-dimensional complex linear space consisting of spinors
	that are polarized in the x-direction.
	Package: VQM`Spinors`.";

QyBasis::usage = "Basis in the two-dimensional complex linear space consisting of spinors
	that are polarized in the y-direction.
	Package: VQM`Spinors`.";

QzBasis::usage = "Basis in the two-dimensional complex linear space consisting of spinors
	that are polarized in the z-direction. By default, this is the same as $QSpinBasis.
	Package: VQM`Spinors`.";

QSetSpinBasis::usage = "QSetSpinBasis[{vec1,vec2}] defines QSpinBasis to consist of two complex
	unit vectors in the direction of vec1 and vec2. The complex 2-vectors vec1 and vec2 must
	be orthogonal.
	Package: VQM`Spinors`.";

QSpinBasis::usage = "The basis in C^2. Initially equal to $QSpinBasis. Setting can be changed
	with the command QSetSpinBasis.
	Package: VQM`Spinors`.";

QUseBasis::usage = "Option for QSpinorToComponents etc, specifying the basis in C^2. Default is QSpinBasis.
	Package: VQM`Spinors`.";

QSpinorToComponents::usage = "QSpinorToComponents[spinor] returns the components of the complex vector 'spinor'
	(which is given with respect to the reference basis $QSpinBasis) with respect to
	the basis QSpinBasis (or the basis specified by the QUseBasis-option.
	Hence QSpinorToComponents[spinor].QSpinBasis == spinor.
	Package: VQM`Spinors`.";

QComponentsToSpinor::usage = "QComponentsToSpinor[{c1,c2}] gives the spinor in $QSpinBasis
	that has the components c1 and c2 with respect to QSpinBasis (or the basis specified by the QUseBasis-option).
	QComponentsToSpinor[{c1,c2}] == {c1,c2}.QSpinBasis.
	Package: VQM`Spinors`.";

QProjectUp::usage = "QProjectUp[spinor] projects a spinor into the direction of the first basis
	vector of QSpinBasis (i.e., into the current 'up-direction').
	Package: VQM`Spinors`.";

QProjectDown::usage = "QProjectDown[spinor] projects a spinor into the direction of the second basis
	vector of QSpinBasis (i.e., into the current 'down-direction').
	Package: VQM`Spinors`.";

QProjection::usage := "QProjection[psi] projects onto the one-dimensional eigenspace defined by
	a spinor psi. Hence QProjection[psi] = |psi> <psi|. Here psi = {psi1,psi2{ is
	a spinor consisting of two complex components.
	Package: VQM`Spinors`."

QProbabilityUp::usage = "QProbabilityUp[spinor] gives the probability that the spinor has spin-up with
	respect to the QSpinBasis.
	Package: VQM`Spinors`.";

QProbabilityDown::usage = "QProbabilityUp[spinor] gives the probability that the spinor has spin-down with
	respect to the QSpinBasis.
	Package: VQM`Spinors`.";

QSpinorBarDiagram::usage = "QSpinorBarDiagram[spinor,opts] gives a bar-diagram showing the components of the spinor
	with respect to the current QSpinBasis. The options are Options[Graphics].
	Package: VQM`Spinors`.";

QPauliSigma1::usage = "Defines the Pauli matrix sigma_1.
	Package: VQM`Spinors`.";

QPauliSigma2::usage = "Defines the Pauli matrix sigma_2.
	Package: VQM`Spinors`.";

QPauliSigma3::usage = "Defines the Pauli matrix sigma_3.
	Package: VQM`Spinors`.";

QPauliSigmaV::usage = "Defines a 'vector' whose components are the 
	three Pauli matrices.
	Package: VQM`Spinors`.";

QIdentity2::usage = "Shortcut for the two-by-two identity matrix.
	Package: VQM`Spinors`.";

QVectorToHermitianMatrix::usage = "QVectorToHermitianMatrix[{k0,k1,k2,k3}] = QIdentity2 * k0 + Sum[QPauliSigma_i * k_i, {i,1,3}].
	Converts a four-dimensional real vector into a Hermitian 2x2 matrix.
	Package: VQM`Spinors`.";

QVectorToDensityMatrix::usage = "QVectorToDensityMatrix[{k1,k2,k3}] = 1/2 (QIdentity2+Sum[QPauliSigma_i * k_i, {i,1,3}]).
	Converts a three-dimensional real vector into a Hermitian 2x2 matrix with
	trace 1 (a density matrix). Package: VQM`Spinors`.";

QSpinHamiltonian::usage = "QSpinHamiltonian[{k0,k1,k2,k3}] = QIdentity2 * k0 + Sum[QPauliSigma_i * k_i, {i,1,3}].
	Package: VQM`Spinors`.";

QSpinorUp::usage = "QSpinorUp[{k1,k2,k3}] is the two-dimensional complex vector
	describing 'spin up' in the direction defined by k = {k1,k2,k3}.
	Normalized eigenvector of (QPauliSigmaV . k) belonging to eigenvalue +|k|.
	Package: VQM`Spinors`.";

QSpinorDown::usage = "QSpinorDown[{k1,k2,k3}] is the two-dimensional complex vector
	describing 'spin down' in the direction defined by k = {k1,k2,k3}.
	Normalized eigenvector of (QPauliSigmaV . k) belonging to eigenvalue -|k|.
	Package: VQM`Spinors`.";
	
QConjSpinorUp::usage  = "QConjSpinorUp[{k1,k2,k3}].
	Package: VQM`Spinors`.";

QConjSpinorDown::usage  = "QConjSpinorDown[{k1,k2,k3}].
	Package: VQM`Spinors`.";

QSpinorHarmonicUp::usage = "QSpinorHarmonicUp[j,m,{theta,phi}] gives the spherical harmonic spinor
	with spin parallel to the orbital angular momentum.
	It is an eigenfunction of the angular momentum operators J^2, L^2, J_3, and L.S.
	Package: VQM`Spinors`.";

QSpinorHarmonicDown::usage = "QSpinorHarmonicDown[j,m,{theta,phi}] gives the spherical harmonic spinor
	with spin antiparallel to the orbital angular momentum.
	It is an eigenfunction of the angular momentum operators J^2, L^2, J_3, and L.S.
	Package: VQM`Spinors`.";

QSpinorToVector::usage = "QSpinorToVector[{z1,z2}] defines a real three-dimensional
	vector in the spin-up direction of psi = {z1,z2}. The vector is obtained from the
	expectation value of QPauliSigmaV in the state defined by psi = {z1,z2}.
	The length of the vector is Abs[psi]^2.
	The spinor {z1,z2} must be nonzero.
	Package: VQM`Spinors`.";

QVectorToSpinor::usage = "QVectorToSpinor[{x,y,z}] defines a spinor psi which is spin-up in the
	direction defined by {x,y,z}. The norm of the spinor is related to the norm of the
	vector by |{x,y,z}|=Abs[psi]^2.
	Package: VQM`Spinors`.";

QVectorLength::usage = "QVectorLength is an option for QSpinorToVector and QVectorToSpinor.
	Default is QVectorLength->2 where the vector has the length |spinor|^2 and the spinor has
	|spinor| = Sqrt[|vector|]. QVectorLength->1 gives
	a unit vector resp. spinor, and QVectorLength->3 generates the vector with length |spinor|
	and the spinor with length |vector|.
	Package: VQM`Spinors`.";

QExtractPhase::usage = "QExtractPhase[psi] returns a real number arg for any nonzero psi={z1,z2}.
	The number arg is a phase determined from comparison with QSpinorUp (with respect
	to the direction defined by psi). Since QSpinorUp is defined with a real first component,
	arg is just the argument of z1.
	For a normalized spinor psi, we have psi=Exp[I arg] QSpinorUp[QSpinorToVector[psi]].
	Package: VQM`Spinors`.";
	
QHermitianMatrixToVector::usage = "QHermitianMatrixToVector[{{z1,z2},{z3,z4}}] converts a 
	Hermitian two-by-two matrix into a real, four-dimensional vector, by writing the matrix
	as a linear combination of the identity and the Pauli matrices.
	Package: VQM`Spinors`.";

QDensityMatrixToVector::usage = "QDensityMatrixToVector[{{z1,z2},{z3,z4}}] converts a 
	Hermitian two-by-two matrix into a real, three-dimensional vector, assuming that
	the matrix is a density matrix (i.e., has trace 1).
	Package: VQM`Spinors`.";

QSpinorToArrow::usage = "QSpinorToArrow[pt,spinor] gives a list of lines representing an arrow
	that corresponds to the spinor as defined in QSpinorToVector. The argument pt is optional
	and defaults to {0,0,0}.
	Package: VQM`Spinors`.";

QRotationSO3::usage = "QRotationSO3[3vector] is an orthogonal 3 by 3 matrix that rotates a vector
	around the axis defined by the direction of 3vector through an angle defined by the size
	of 3vector.
	Package: VQM`Spinors`.";

QRotationSU2::usage = "QRotationSU2[3vector] is a unitary 2 by 2 matrix that rotates a spinor
	around the axis defined by the direction of 3vector through an angle defined by the size
	of 3vector.
	Package: VQM`Spinors`.";

QVisualizeSpinor::usage = "QVisualizeSpinor[spinor] converts a spinor into a magnetic needle
	and displays it together with other graphics elements whose appearance is controlled
	by the options QDrawUnitSphere, QDrawAxes, QCoordinateCube, QCoordinateCircles.
	Behavior for QNeedleStyle->True (default): The needle points from the coordinate origin
	to the point QSpinorToVector[spinor]. The upper (lower) half of the needle has a
	color determined from the upper (lower) component of the spinor via QComplexToColor
	(from the package VQM`ComplexPlot`).
	Behavior for QNeedleStyle->False. An arrow is shown instead of the needle.
	Giving the option QHeadColor->QExtractPhase colors the head of the arrow by the phase
	of the first component of the spinor.
	You can give all options from QVisualizeVector (from the package VQM`VisualizeVector`).
	Package: VQM`Spinors`.";

QVisualizeDensityMatrix::usage = "QVisualizeDensityMatrix[matrix] converts a Hermitian matrix
	with trace 1 into an arrow graphics (with the help of QDensityMatrixToVector)
	and displays it together with other graphics elements whose appearance is controlled
	by the options QDrawUnitSphere, QDrawAxes, QCoordinateCube, QCoordinateCircles.
	Package: VQM`Spinors`.";

(*-----------------------------------*)
Begin["`Private`"];
(*-----------------------------------*)

(* constant matrices: *)

QIdentity2 = IdentityMatrix[2];
QPauliSigma1 = {{0,1},{1,0}};
QPauliSigma2 = {{0,-I},{I,0}};
QPauliSigma3 = {{1,0},{0,-1}};
QPauliSigmaV = {QPauliSigma1, QPauliSigma2, QPauliSigma3};

(* eigenvectors of Pauli matrices: *)

QSpinorUp[vec_] := {0,0}/;vec=={0,0,0}
QSpinorDown[vec_] := {0,0}/;vec=={0,0,0}
QSpinorUp[{k1_,k2_,k3_}] := {UnitStep[k3],UnitStep[-k3]}/;(k1==0 && k2==0);
QSpinorDown[{k1_,k2_,k3_}] := {-UnitStep[-k3],UnitStep[k3]}/;(k1==0 && k2==0);
QConjSpinorUp[{k1_,k2_,k3_}] := QSpinorUp[{0,0,k3}]/;(k1==0 && k2==0);
QConjSpinorDown[{k1_,k2_,k3_}] := QSpinorDown[{0,0,k3}]/;(k1==0 && k2==0);

QSpinorUp[k_?VectorQ] := (* this definition has a discontinuity for k on south pole *)
	Module[{absk = Sqrt[k.k]},
		{absk + k[[3]], k[[1]] + I*k[[2]]}/Sqrt[2*absk*(absk + k[[3]])]//Simplify
	];

(* (*This would give a definition where the lower component is continuous *)
QSpinorUp[k_?VectorQ] :=
	Module[{absk = Sqrt[k.k]},
		{k[[1]] - I*k[[2]], absk - k[[3]]}/Sqrt[2*absk*(absk - k[[3]])]//Simplify
	];
*)

QSpinorDown[k_?VectorQ] :=
	Module[{absk = Sqrt[k.k]},
		{-absk + k[[3]], k[[1]] + I*k[[2]]}/Sqrt[2*absk*(absk - k[[3]])]//Simplify
	];

QConjSpinorUp[k_?VectorQ] :=
	Module[{absk = Sqrt[k.k]},
		{absk + k[[3]], k[[1]] - I*k[[2]]}/Sqrt[2*absk*(absk + k[[3]])]//Simplify
	];

QConjSpinorDown[k_?VectorQ] :=
	Module[{absk = Sqrt[k.k]},
		{-absk + k[[3]], k[[1]] - I*k[[2]]}/Sqrt[2*absk*(absk - k[[3]])]//Simplify
	];

QSpinorHarmonicUp[j_,m_,{theta_,phi_}] := 1/Sqrt[2 j]*
	{Sqrt[j+m] SphericalHarmonicY[j-1/2,m-1/2,theta,phi],
	 Sqrt[j-m] SphericalHarmonicY[j-1/2,m+1/2,theta,phi]}

QSpinorHarmonicDown[j_,m_,{theta_,phi_}] := 1/Sqrt[2 j + 2]*
	{Sqrt[j+1-m] SphericalHarmonicY[j+1/2,m-1/2,theta,phi],
	-Sqrt[j+1+m] SphericalHarmonicY[j+1/2,m+1/2,theta,phi]}

(* some useful operations: *)

QNorm[spinor_] :=	Sqrt[Conjugate[spinor].spinor];
QNormalize[spinor_] :=	spinor/QNorm[spinor];

(* The default basis in C^2 is: *)

basespinor1 = {1,0};
basespinor2 = {0,1};
$QSpinBasis = {basespinor1,basespinor2}; (*reference basis, cannot be changed by user*)

QxBasis = {QSpinorUp[{1,0,0}],QSpinorDown[{1,0,0}]};
QyBasis = {QSpinorUp[{0,1,0}],QSpinorDown[{0,1,0}]};
QzBasis = $QSpinBasis;

(* Choose another basis in C^2: *)

QSpinBasis := {basespinor1,basespinor2}; (*initially, QSpinBasis=$QSpinBasis, but can be changed by user *)

QSetSpinBasis[{spinor1_,spinor2_}] :=
	(
		basespinor1 = QNormalize[spinor1];
		basespinor2 = QNormalize[spinor2];
		If[Conjugate[basespinor1].basespinor2 != 0, Message[onb::warning]];
		Return[]
	);

(* Basic operations with respect to a basis: *)

QSpinorToComponents[spinor_,opts___?OptionQ] :=
	Module[{c1,c2, basevectors},
		basevectors = QUseBasis /. Join[{opts},{QUseBasis->QSpinBasis}];
		c1 = Conjugate[basevectors[[1]]].spinor;
		c2 = Conjugate[basevectors[[2]]].spinor;
		Return[{c1,c2}]
	];
		(* components with respect to QSpinBasis *)

QComponentsToSpinor[components_, opts___?OptionQ] := 
	Module[{basevectors},
		basevectors = QUseBasis /. Join[{opts},{QUseBasis->QSpinBasis}];
		components.basevectors
	];
		(* gives the spinor having components {c1,c2} in QSpinBasis *) 

QProjectUp[spinor_,opts___?OptionQ] :=
	Module[{c},
		c = QSpinorToComponents[spinor,opts][[1]];
		If[c==0, Return[{0,0}], Return[QComponentsToSpinor[{c,0}, opts]]
		]
	];

QProjectDown[spinor_,opts___?OptionQ] :=
	Module[{c},
		c = QSpinorToComponents[spinor,opts][[2]];
		If[c==0, Return[{0,0}], Return[QComponentsToSpinor[{0,c}, opts]]
		]
	];

QProjection[{a_,b_}] :=
	{{Abs[a]^2,a Conjugate[b]},{b Conjugate[a], Abs[b]^2}};

QProbabilityUp[spinor_,opts___?OptionQ] :=
	Abs[QSpinorToComponents[QNormalize[spinor],opts][[1]]]^2;
	
QProbabilityDown[spinor_,opts___?OptionQ] :=
	Abs[QSpinorToComponents[QNormalize[spinor],opts][[2]]]^2;

(* operations with spinors and sigma matrices *)

QVectorToHermitianMatrix[k_?VectorQ] :=
	Simplify[(QIdentity2*k[[1]] + QPauliSigma1*k[[2]] +
		QPauliSigma2*k[[3]] + QPauliSigma3*k[[4]])]/;Length[k]==4;

QVectorToHermitianMatrix[k_?VectorQ] := QVectorToHermitianMatrix[{0,k[[1]],k[[2]],k[[3]]}]/;Length[k]==3;

QVectorToDensityMatrix[k_?VectorQ] := 1/2(QVectorToHermitianMatrix[{1,k[[1]],k[[2]],k[[3]]}])/;Length[k]==3;

QSpinHamiltonian[k_?VectorQ] := 2*QVectorToDensityMatrix[k]-QIdentity2/;Length[k]==3;

QHermitianMatrixToVector[P_?MatrixQ] :=
	Module[{p0,p1,p2,p3},
		If[And @@ NumericQ /@ Flatten[P],
			If[Transpose[Conjugate[P]]=!=P,Message[matrix::nonhermitian];Return["Null"]];
    	];
    	p0 = Tr[P]/2;
    	p1 = Tr[P.QPauliSigma1]/2;
    	p2 = Tr[P.QPauliSigma2]/2;
    	p3 = Tr[P.QPauliSigma3]/2;
		Return[{p0,p1,p2,p3}]
	];

QDensityMatrixToVector[rho_?MatrixQ] :=
	Module[{P,p1,p2,p3},
		If[And @@ NumericQ /@ Flatten[rho],
			If[Transpose[Conjugate[rho]]=!=rho,Message[matrix::nonhermitian];Return["Null"]];
			If[Tr[rho] != 1,Message[matrix::nondensity];Return["Null"]];
		];
		P=rho-QIdentity2/2;
    	p1 = Tr[P.QPauliSigma1];
    	p2 = Tr[P.QPauliSigma2];
    	p3 = Tr[P.QPauliSigma3];
		Return[{p1,p2,p3}]
	];

Options[QSpinorToVector] = {QVectorLength->2};

QSpinorToVector[psi_?VectorQ,opts___?OptionQ] :=
	Module[{absf, conj = Conjugate[psi], len},
		len = QVectorLength/.{opts}/.Options[QSpinorToVector];
		If[len==2,absf=1];
		If[len==1,absf=conj.psi];
		If[len==3,absf=Sqrt[conj.psi]];
		If[Re[Chop[absf]]==0.,absf=1];
		v1 = conj.(QPauliSigma1.psi)/absf;
		v2 = conj.(QPauliSigma2.psi)/absf;
		v3 = conj.(QPauliSigma3.psi)/absf;
		Return[Re[{v1, v2, v3}]];
	]/;Length[psi]==2;

QVectorToSpinor[kv_?VectorQ,opts___?OptionQ] :=
	Module[{absf, len, lenkv = Sqrt[kv.kv]},
		len = QVectorLength/.{opts}/.Options[QSpinorToVector];
		If[len==2,absf=Sqrt[lenkv]];
		If[len==3,absf=lenkv];
		If[len==1,absf=1];
		absf QSpinorUp[kv]
	]/;Length[kv]==3;

(*
QExtractPhase[psi_?VectorQ] :=
	Module[{psi1},
		psi1 = QSpinorUp[Chop[QSpinorToVector[psi]]];
		If[     Chop[ psi1[[1]] ] != 0. && Chop[ psi[[1]] ] != 0., Arg[ psi[[1]]/psi1[[1]] ],
			If[ Chop[ psi1[[2]] ] != 0. && Chop[ psi[[2]] ] != 0., Arg[ psi[[2]]/psi1[[2]] ] ]
		]
	]/;Length[psi]==2;
*)
(* As long as QSpinorUp has a real-valued first component,
   we can use the following simpler definition: *)

QExtractPhase[psi_?VectorQ] :=
	Module[{psi1 = psi[[1]]},
		If[ Chop[ psi1 ] != 0, Arg[ psi[[1]] ], 0, 0] 
	]/;Length[psi]==2;

(* Rotation matrices *)

QRotationSO3[vec_?VectorQ] :=
	Module[{alpha = Sqrt[vec.vec], n, R},
		If[Length[vec] != 3, Message[dim::three]; Return[]];
		If[vec=={0,0,0},Return[IdentityMatrix[3]]];
		n = vec/alpha;
		R = Table[
			KroneckerDelta[i,k]*Cos[alpha] +
			n[[i]] n[[k]] (1 - Cos[alpha]) -
			Sum[Signature[{i,k,m}] n[[m]] Sin[alpha],{m,1,3}],
		{i,1,3},{k,1,3}];
		Return[Simplify[R,alpha>0]]
	];

QRotationSU2[vec_?VectorQ] :=
	Module[{alpha = Sqrt[vec.vec], n, U},
		If[Length[vec] != 3, Message[dim::three]; Return[]];
		If[vec=={0,0,0},Return[QIdentity2]];
		n = vec/alpha;
		U = Cos[alpha/2] QIdentity2 - I Sin[alpha/2] (QPauliSigma1*n[[1]]+QPauliSigma2*n[[2]]+QPauliSigma3*n[[3]]);
		Return[Simplify[U,alpha>0]]
	];


QSpinorToArrow[pt_?VectorQ, spinor_?VectorQ, opts___?OptionQ] :=
	Module[{vec, coloropt, color, vshp, opts2},
		coloropt = QHeadColor/.{opts}/.Options[QVectorToArrow];
		vec = QSpinorToVector[spinor, opts];
		opts2=Sequence @@ Select[{opts},!MemberQ[{QHeadColor},First[#]]&];

		If[	coloropt === QExtractPhase,
			   color = Hue[QExtractPhase[spinor]/(2 Pi)];
			   QVectorToArrow[pt,pt+vec, opts2, QHeadColor->Automatic],

			(*else*)QVectorToArrow[pt,pt+vec, opts2,
							QHeadColor ->QComplexToColor[spinor[[1]], FilterRules[Flatten[{opts}], Options @ QComplexToColor]],
							QShaftColor->QComplexToColor[spinor[[2]], FilterRules[Flatten[{opts}], Options @ QComplexToColor]] ]
		]
	]/;(Length[pt]==3 && Length[spinor]==2);

QSpinorToArrow[spinor_?VectorQ, opts___?OptionQ] :=
	QSpinorToArrow[{0,0,0},spinor,opts]/;Length[spinor]==2;

(* Visualization commands *)
Options[QVisualizeSpinor] = Join[Flatten[{QNeedleStyle -> True}],
            Select[Flatten[Options[QVisualizeVector]],
                !MemberQ[{QNeedleStyle},First[#]]&]];

QVisualizeSpinor[psi_?VectorQ, opts___?OptionQ] :=
	Module[{vec, coloropt, color, vshp},
		{ coloropt, vshp} = {QHeadColor, QNeedleStyle}/.{opts}/.Options[QVisualizeSpinor];
		vec = QSpinorToVector[psi, opts];
		If[	coloropt === QExtractPhase
		    , 
		    color = Hue[QExtractPhase[psi]/(2 Pi)];
			QVisualizeVector[vec, opts, QHeadColor->color]
			,
			(*else*) 
			(* Glow would restore the old behaviour ..., but this is not advisable *)
			headColor  = (*Glow @ *)(QComplexToColor[psi[[1]], FilterRules[Flatten[{opts}], Options @ QComplexToColor]]);
			shaftColor = (*Glow @ *)(QComplexToColor[psi[[2]], FilterRules[Flatten[{opts}], Options @ QComplexToColor]]);
			QVisualizeVector[vec, opts,  QNeedleStyle -> True,
						QHeadColor -> headColor,
						QShaftColor-> shaftColor
			]
		]
	];


Options[QVisualizeDensityMatrix] = Options[QVisualizeVector];

QVisualizeDensityMatrix[P_?MatrixQ, opts___?OptionQ] :=
	Module[{vec},
		vec = QDensityMatrixToVector[P];
    	QVisualizeVector[vec,opts]
	];

QSpinorBarDiagram[spinor_,opts___?OptionQ]:=
	Module[{c1,c2,abs1,abs2,clr1,clr2},
		{c1,c2} = QSpinorToComponents[spinor,opts];
		abs1 = Abs[c1]; abs2 = Abs[c2];
		clr1 = If[abs1!=0,Hue[Arg[c1]/(2 Pi)],Hue[0],Hue[0]];
		clr2 = If[abs2!=0,Hue[Arg[c2]/(2 Pi)],Hue[0],Hue[0]];
		Show[
			Graphics[{
				{clr1, Rectangle[{0.05,0},{0.95,abs1}] },
				{clr2, Rectangle[{1.05,0},{1.95,abs2}] }
			}], Flatten[{
			    FilterRules[Flatten[{opts}], Options @ Graphics],
				Frame->True, PlotRange->{{-0.1,2.1},{-0.1,1.1}},
				FrameTicks->{{{0.5,"\!\(c\_1\)"},{1.5,"\!\(c\_2\)"}},Automatic,None,Automatic},
				Epilog->{GrayLevel[0], Line[{{0,0},{2,0}}],Line[{{0,1},{2,1}}]}
				}]
		]
	];


(* Messages: *)

matrix::nonhermitian := "Matrix must be Hermitian.";
matrix::nondensity := "Matrix must have trace 1.";
dim::three := "Argument must be a three-dimensional vector.";
onb::warning := "You specified a basis that is not orthonormal.";

(*-----------------------------------*)
End[];      (* end `Private` Context *)
(*-----------------------------------*)
(*
Protect @@ VQM`Spinors`Private`Symbols;
*)
(*-----------------------------------*)
EndPackage[]; (* end package Context *)
(*-----------------------------------*)

If[VQMmsgon, On[General::"spell"]];
If[VQMmsgon1, On[General::"spell1"]];

(* :Examples:

QVisualizeSpinor[QNormalize[{-1+I,-I/3-1}], QDrawUnitSphere->15, QCoordinateCube->True];

QVisualizeSpinor[QNormalize[{I,-I}], QNeedleStyle->False, QArrowShaft->True,
	Lighting->False, QArrowShape->{10,1/2,1/2,1/2}];
	
*)
