(* :Title:	Harmonic Oscillator *)

(* :Name:	VQM`Oscillator` *)

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
Defines eigenfunctions and squeezed states for the
one dimensional harmonic oscillator in quantum mechanics.
The time-dependent functions are solutions of the Schroedinger equation
I D[Psi,t] = (1/2m) D[Psi,{x,2}] + (m w^2/2) x^2 Psi.
That is, we use units where hbar is scaled to 1. The value of m is given by
the constant $QOscillatorMass, w is given by $QOscillatorFrequency (defaults 1).
*)

(* :Date:	2006-07-20 *)

(* :Package Version:		2.0 *)

(* :Mathematica Version:	6.0.1 *)

(* :Keywords:
	Free Time Evolution, Gaussian, Wave Packets
*)
    
VQMmsgon  = Head[General::"spell"]  =!= $Off;
VQMmsgon1 = Head[General::"spell1"] =!= $Off;

Off[General::spell1,General::spell];

(*-----------------------------------*)
BeginPackage[ "VQM`Oscillator`"];
(*-----------------------------------*)

ClearAll[$QOscillatorFrequency, $QOscillatorMass];

VQM`Oscillator`Private`Symbols = Hold[
	QOscillatorHamiltonian, QOscillatorEnergy, QOscillatorFunction,
	QOscillatorFunctionT, QOscillatorGaussian, QOscillatorBarDiagram,
	QOscillatorFrequency, QOscillatorMass];

Unprotect @@ VQM`Oscillator`Private`Symbols;
ClearAll @@ VQM`Oscillator`Private`Symbols;

 
QOscillatorHamiltonian::usage = "QOscillatorHamiltonian1D[f[x],x, opts] applies
the Hamiltonian of the one dimensional harmonic oscillator to the function
f[x]. The option QOscillatorFrequency->w (default $QOscillatorFrequency=1)
defines the frequency of the harmonic oscillator.
The option QOscillatorMass->m (default $QOscillatorMass = 1) defines the mass
of the particle (we use units with hbar = 1). Package: VQM`Oscillator`.";

QOscillatorEnergy::usage = "QOscillatorEnergy[n, opts] gives the energy of
the n-th eigenfunction of the harmonic oscillator in one dimension.
The option QOscillatorFrequency->w defines the frequency
of the harmonic oscillator. Package: VQM`Oscillator`.";

QOscillatorFunction::usage = "QOscillatorFunction[n,x,opts] is the n-th
eigenfunction of the harmonic oscillator in one dimension. Possible options
are QOscillatorMass and QOscillatorFrequency. Package: VQM`Oscillator`.";

QOscillatorFunctionT::usage = "QOscillatorFunction[n,x,t,opts] describes the 
time evolution of the n-th
eigenfunction of the harmonic oscillator in one dimension.
This function is just a time-dependent phase
factor times QOscillatorFunction[n,x]. Package: VQM`Oscillator`.";

QOscillatorGaussian::usage = "QOscillatorGaussian[x,t,x0,p0,a,opts] describes
the time evolution of a Gaussian initial function in the field of a harmonic
oscillator potential. The harmonic oscillator is characterized by the options
QOscillatorFrequency->w (default 1) and QOscillatorMass->m (default 1).
The arguments x0, p0, a are also optional (defaults x0=0, p0=1, a=1).
x0 is the average initial position, p0 is the average initial momentum.
a describes the width of the initial position distribution. Package: VQM`Oscillator`.
";

QOscillatorBarDiagram::usage = "QOscillatorBarDiagram[f[x],{x,a,b},n1,n2, opts] plots the
energy representation of a given function in the basis of harmonic
oscillator eigenfunctions. Numerical region for the determination of the expansion
coefficients is the interval (a,b). The graph shows the expansion coefficients c_n
for n between n1 and n2. The interval (a,b) should be large enough so that
all eigenfunctions with quantum numbers larger than n2 are essentially zero
outside that interval. The short version QOscillatorBarDiagram[f[x],{x}] uses the Default values
a=-Infinity, b=Infinity, n1=0, n2=10. Package: VQM`Oscillator`.";

QOscillatorFrequency::usage = "QOscillatorFrequency is an option for QOscillatorHamiltonian,
QOscillatorEnergy, QOscillatorFunction, QOscillatorFunctionT, QOscillatorGaussian,QOscillatorBarDiagram.
QOscillatorFrequency->w sets the frequency of the harmonic oscillator (the coupling
constant of the oscillator potential is m w^2/2). Package: VQM`Oscillator`.";

QOscillatorMass::usage = "QOscillatorMass is an option for QOscillatorHamiltonian,
QOscillatorFunction, QOscillatorFunctionT, QOscillatorGaussian,QOscillatorBarDiagram.
QOscillatorMass->m sets the mass of the harmonic oscillator. Package: VQM`Oscillator`.";

$QOscillatorFrequency::usage = "$QOscillatorFrequency is the default value for the option
QOscillatorFrequency. Can be redefined by the user. Package: VQM`Oscillator`.";

$QOscillatorMass::usage = "$QOscillatorMass is the default value for the option
QOscillatorMass. Can be redefined by the user. Package: VQM`Oscillator`.";

Begin["`Private`"];

OscillatorOptions := {QOscillatorFrequency->$QOscillatorFrequency,
QOscillatorMass->$QOscillatorMass};

Options[QOscillatorHamiltonian] = OscillatorOptions;
Options[QOscillatorEnergy] = OscillatorOptions;
Options[QOscillatorFunction] = OscillatorOptions;
Options[QOscillatorGaussian] = OscillatorOptions;

QOscillatorHamiltonian1D[f_,x_Symbol, opts___?OptionQ] :=
	Module[{m,w},
		{m,w} = {QOscillatorMass,QOscillatorFrequency} /. Join[{opts}, OscillatorOptions];
		func = Function[{x},f];
		-1/(2 m) D[func[x],{x,2}] + m w^2 x^2 func[x]/2
   ];

QOscillatorEnergy[n_Integer,opts___?OptionQ] :=
	Module[{w},
		w = QOscillatorFrequency /. Join[{opts}, OscillatorOptions];
		w (n + 1/2)]

QOscillatorFunction[n_,x_,opts___?OptionQ] :=
	Module[{m,w},
		{m,w} = {QOscillatorMass,QOscillatorFrequency} /. Join[{opts}, OscillatorOptions];
		(m w/Pi)^(1/4) HermiteH[n,x Sqrt[m w]]*Exp[-m w x^2/2]/Sqrt[2^n*n!]
	]

QOscillatorFunction[n_,x_,t_,opts___?OptionQ] :=
	QOscillatorFunction[n,x, opts]*Exp[-I QOscillatorEnergy[n,opts] t]

QOscillatorGaussian[x_,t_,x0_:0,p0_:0,a_:1,opts___?OptionQ] :=
	Module[{m,w},
		{m,w} = {QOscillatorMass,QOscillatorFrequency} /. Join[{opts}, OscillatorOptions];
		(m w)^(1/4) $psi[x*Sqrt[m w], t w, x0*Sqrt[m w], p0/Sqrt[m w], a/(m w)]
	]

QOscillatorBarDiagram[f_, {x_Symbol, x1_: - Infinity, x2_:Infinity}, n1_:0, n2_:10, opts___] := 
	Module[{coeff, coloredLine, lineTable, up, upper, lower, func, w},
		func = Function[{x}, f]; 
		w = QOscillatorFrequency /. Join[{opts}, OscillatorOptions];
		Do[coeff[n] = NIntegrate[QOscillatorFunction[n, x, opts]*
			func[x], {x, x1, x2}], {n, n1, n2}]; 
		coloredLine[n_] := {Thickness[0.01], 
			Hue[Arg[coeff[n]]/(2*Pi)], 
			Line[{{w (n + 1/2), Abs[coeff[n]]},
				{w (n + 1/2), 0.}}]};
		lineTable = 
			Graphics[Table[coloredLine[n], {n, n1, n2}]]; 
		up = Max[Table[Abs[coeff[n]], {n, n1, n2}]]; 
		upper = up + up/10;
		lower = -(up/10); 
		Show[lineTable,
			Frame -> True,
			PlotRange -> {w {n1, n2 + 1}, {lower, upper}}]
  ]

(* auxiliary functions (internal) *)


$psi1[x_,t_,x0_,p0_,a_] :=
	a^(1/4)/(E^((-2*x*(I*p0 + a*x0) + a*(x^2 + x0^2)*Cos[t] +
	I*(p0^2 + x^2 - (2*I)*a*p0*x0)*Sin[t])/
    (2*(Cos[t] + I*a*Sin[t])))*Pi^(1/4)*Sqrt[Cos[t] + I*a*Sin[t]])
    
    
(* adjust for the branch-cut of the square root *)
$sig[t_] := If[OddQ[Floor[(t/Pi-1)/2]],1,-1]

$psi[x_, t_, x0_, p0_, a_] := $sig[t] $psi1[x,t,x0,p0,a]

check[opts___] :=
   Module[{m,w},
      {m,w} = {QOscillatorMass,QOscillatorFrequency} /. Join[{opts}, OscillatorOptions];
      If[m<=0 || w<=0,Message[Oscillator::nonneg]];
      Return[{m,w}]]


Oscillator::nonneg		:= "Mass and oscillator frequency must be positive.";

(*-----------------------------------*)
End[];      (* end `Private` Context *)
(*-----------------------------------*)

Protect @@ VQM`Oscillator`Private`Symbols;

(* not protected symbols: *)
$QOscillatorFrequency=1;
$QOscillatorMass=1;

(*-----------------------------------*)
EndPackage[]; (* end package Context *)
(*-----------------------------------*)

If[VQMmsgon, On[General::"spell"]];
If[VQMmsgon1, On[General::"spell1"]];
