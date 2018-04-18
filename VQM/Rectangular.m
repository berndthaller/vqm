(* :Title:   Rectangular potentials *)

(* :Name:    VQM`Rectangular` *)

(* :Copyright: Copyright 2007 Bernd Thaller *)

(* :Author:  Bernd Thaller,
             Institute of Mathematics,
             University of Graz,
             A-8010 Graz
             bernd.thaller@kfunigraz.ac.at
*)

(* :Source:
	Visual Quantum Mechanics
	Springer-Verlag New York, 2000,
*)

(* :Summary:
This package defines the solutions of the one-dimensional
Schroedinger equation in the presence of rectangular
potentials (steps, wells, barriers). *)

(* :Date:    Jul 20, 2007 *)

(* :Package Version:        0.8 development version *)

(* :Mathematica Version:    6.0 *)

(* :Keywords:
    Rectangular, Potential Well, Potential Step, Potential Barrier
*)

VQMmsgon  = Head[General::"spell"]  =!= $Off;
VQMmsgon1 = Head[General::"spell1"] =!= $Off;

Off[General::spell1,General::spell];


(*-----------------------------------*)
BeginPackage[ "VQM`Rectangular`" ]                
(*-----------------------------------*)


VQM`Rectangular`Private`Symbols = Hold[
QPlaneWaveToRight, QPlaneWaveToLeft,
QReflectionCoefficientJump, QTransmissionCoefficientJump,
QReflectionCoefficientWell, QTransmissionCoefficientWell,
QRadZeroReflWell, QEnZeroReflWell, QPotZeroReflWell,
QTransitionMatrix, QSolutionWellToRight,
QPsiEvenWell, QPsiOddWell,
QDetEvenWell, QDetOddWell,
QCriticalRadiusEvenWell, QCriticalRadiusOddWell,
QCriticalDepthEvenWell, QCriticalDepthOddWell];

Unprotect @@ VQM`Rectangular`Private`Symbols;
ClearAll @@ VQM`Rectangular`Private`Symbols;

 
QPlaneWaveToRight::usage = "QPlaneWaveToRight[...]...Package: VQM`Rectangular`.";

QPlaneWaveToLeft::usage = "QPlaneWaveToLeft[...]... Package: VQM`Rectangular`.";

QReflectionCoefficientJump::usage = "QReflectionCoefficientJump[...]... Package: VQM`Rectangular`.";

QTransmissionCoefficientJump::usage = "QTransmissionCoefficientJump[...]... Package: VQM`Rectangular`.";

QReflectionCoefficientWell::usage = "QReflectionCoefficientWell[...]... Package: VQM`Rectangular`.";

QTransmissionCoefficientWell::usage = "QTransmissionCoefficientWell[...]... Package: VQM`Rectangular`.";

QTransitionMatrix::usage = "QTransitionMatrix[En,V1,V2,s] is the transition matrix as a
	function of the energy for a potential-jump from V[x]=V1 (for x<s) to V[x]=V2 (for x>s).
	Package: VQM`Rectangular`.";

QSolutionWellToRight::usage = "QSolutionWellToRight[...]... Package: VQM`Rectangular`.";

QPsiEvenWell::usage = "QPsiEvenWell[...]... Package: VQM`Rectangular`.";

QPsiOddWell::usage = "QPsiOddWell[...]... Package: VQM`Rectangular`.";

QDetEvenWell::usage = "QDetEvenWell[...]... Package: VQM`Rectangular`.";

QDetOddWell::usage = "QDetOddWell[...]... Package: VQM`Rectangular`.";

QCriticalRadiusEvenWell::usage = "QCriticalRadiusEvenWell[...]... Package: VQM`Rectangular`.";

QCriticalRadiusOddWell::usage = "QCriticalRadiusOddWell[...]... Package: VQM`Rectangular`.";

QCriticalDepthEvenWell::usage = "QCriticalDepthEvenWell[...]... Package: VQM`Rectangular`.";

QCriticalDepthOddWell::usage = "QCriticalDepthOddWell[...]... Package: VQM`Rectangular`.";

(*-----------------------------------*)
Begin["`Private`"]
(*-----------------------------------*)

(* Energy-momentum relations *)
kappa[En_]	:= Sqrt[-2 En]		/; 	En <= 0
k0[En_, V0_]	:= Sqrt[2 (V0+En)]	/;	-En <= V0
k[En_]		:= Sqrt[2 En]		/;	En > 0

(* plane-wave solutions *)
QPlaneWaveToRight[En_,x_] := 1/Sqrt[2 Pi k[En]] Exp[I k[En] x];
QPlaneWaveToLeft[En_,x_] := 1/Sqrt[2 Pi k[En]] Exp[ - I k[En] x];

QTransitionMatrix[En_,V1_,V2_,s_]:=
	1/(2 Sqrt[k[En-V1] k[En-V2]]) {{(k[En-V1]+ k[En-V2])*
          Exp[I (k[En-V1]- k[En-V2])s], -(k[En-V1]- k[En-V2])*
          Exp[-I (k[En-V1]+ k[En-V2])s]},{-(k[En-V1]- k[En-V2])*
          Exp[I (k[En-V1]+ k[En-V2])s],(k[En-V1]+ k[En-V2])*
          Exp[-I (k[En-V1]- k[En-V2])s]}};

(* scattering coefficients for potential jump from V[x]=0 (x<0) to V[x]=V0 (x>0) *)
QReflectionCoefficientJump[En_, V0_] := (k[En]-k[En-V0])/(k[En]+k[En-V0]);
QTransmissionCoefficientJump[En_, V0_] := 2 Sqrt[k[En] k[En-V0]]/(k[En]+k[En-V0]);

(* scattering coefficients for a symmetric rectangular potential of height V0 and radius R.
V0>0 is a barrier, V0<0 is a well *)

QReflectionCoefficientWell[En_,V0_,R_] :=
      ((-1 + E^(4 I R k[En - V0]))*(k[En]^2 - k[En - V0]^2))/
      (E^(2 I R k[En])*(E^(4 I R k[En - V0])*
      (k[En] - k[En - V0])^2 - (k[En] + k[En - V0])^2));
     
QTransmissionCoefficientWell[En_,V0_,R_] := (-4 k[En] k[En - V0])/
     (E^(2 I R (k[En] - k[En - V0]))*(E^(4 I R k[En - V0])*
     (k[En] - k[En - V0])^2 - (k[En] + k[En - V0])^2));

(* Zeros of the reflection coefficient *)

QRadZeroReflWell[n_, En_, V0_] := n*Pi/(2*Sqrt[2]*Sqrt[En - V0]);
QEnZeroReflWell[n_, V0_, R_] := n^2 Pi^2/(8 R^2) + V0;
QPotZeroReflWell[n_, En_, R_] := En - n^2 Pi^2/(8 R^2);


(* Coefficients for the solution inside the barrier: *)

ca[En_,V0_,R_] :=
	(QTransitionMatrix[En,0,V0,-R].{1,QReflectionCoefficientWell[En,V0,R]})[[1]];
cb[En_,V0_,R_] :=
	(QTransitionMatrix[En,0,V0,-R].{1,QReflectionCoefficientWell[En,V0,R]})[[2]];

(* Define the wave function: *)

QSolutionWellToRight[En_,x_,V0_,R_] :=
	Which[
		x <= -R,
			QPlaneWaveToRight[En,x] +
			QReflectionCoefficientWell[En,V0,R] QPlaneWaveToLeft[En,x],
		x > -R && x <= R,
			ca[En,V0,R] QPlaneWaveToRight[En-V0,x] +
			cb[En,V0,R] QPlaneWaveToLeft[En-V0,x],
		x > R,
			QTransmissionCoefficientWell[En,V0,R] QPlaneWaveToRight[En,x]
	]

(* Bound states in a rectangular well *)

QPsiEvenWell[En_, x_, V0_, R0_] := 1/Sqrt[R0 + 1/Sqrt[-2*En]] *
	Which[
		Abs[x] <= R0,
			Cos[Sqrt[2 (En - V0)] x], 
		x > R0,
			Cos[Sqrt[2 (En - V0)] R0] Exp[-Sqrt[-2 En] (x - R0)], 
     	x < -R0,
			Cos[Sqrt[2 (En - V0)] R0] Exp[ Sqrt[-2 En] (x + R0)],
		True, 0]

QPsiOddWell[En_, x_, V0_, R0_] := 1/Sqrt[R0 + 1/Sqrt[-2*En]] *
	Which[
		Abs[x] <= R0, 
    		Sin[Sqrt[2 (En - V0)] x],
		x > R0,
			Sin[Sqrt[2 (En - V0)] R0] Exp[-Sqrt[-2 En] (x - R0)],
		x < -R0,
			-Sin[Sqrt[2 (En - V0)] R0] Exp[Sqrt[-2 En] (x + R0)],
		True, 0]

(* The above are eigenfunctions if En is a zero of: *)

QDetEvenWell[En_, V0_, R0_] := 
	Sqrt[2 (En - V0)] Sin[Sqrt[2 (En - V0)] R0] - 
		Sqrt[-2 En] Cos[Sqrt[2 (En - V0)] R0]

QDetOddWell[En_, V0_, R0_] := 
	Sqrt[2 (En - V0)] Cos[Sqrt[2 (En - V0)] R0] + 
		Sqrt[-2 En] Sin[Sqrt[2 (En - V0)] R0]

(* New bound states appear at these values of R0 and V0: *)
QCriticalRadiusEvenWell[V0_,n_?IntegerQ] := n Pi/(Sqrt[-2 V0])/;(n>=0 && V0<=0)
QCriticalRadiusOddWell[V0_,n_?IntegerQ] := (2 n - 1) Pi/(2 Sqrt[-2 V0])/;(n>=1 && V0<=0)
QCriticalDepthEvenWell[R0_,n_?IntegerQ] := - n^2 Pi^2 / (2 R0^2) /;(n>=0 && R0>=0)
QCriticalDepthOddWell[R0_,n_?IntegerQ] := - (2 n - 1)^2 Pi^2/(8 R0^2)/;(n>=1 && R0>=0)

(*-----------------------------------*)
End[];      (* end `Private` Context *)
(*-----------------------------------*)

Protect@@VQM`Rectangular`Private`Symbols;

(*-----------------------------------*)
EndPackage[]; (* end package Context *)
(*-----------------------------------*)

If[VQMmsgon, On[General::"spell"]];
If[VQMmsgon1, On[General::"spell1"]];
