(* ::Package:: *)

(* :Title:   FastFourier *)

(* :Name:    VQM`FastFourier` *)

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
This package computes the one-dimensional Fourier transform
of a complex-valued function f[x] via a fast Fourier transform of a
list representing the function values on a suitable space grid.
The package provides methods
to determine the domains of definition and the discretization in Fourier space
from the discretization in position space and vice versa.
This package can be used, e.g., to compute the Fourier transform of a
numerical solution of the Schroedinger equation.

Version 1.4: Changed RotateLeft in QFourierList etc
*)

(* :Date:	2006-07-20 *)

(* :Package Version:		2.0 *)

(* :Mathematica Version:    6.0.1 *)

(* :Keywords:
    Fourier transform, Quantum mechanics, Wavefunction
*)

VQMmsgon  = Head[General::"spell"]  =!= $Off;
VQMmsgon1 = Head[General::"spell1"] =!= $Off;

Off[General::spell1, General::spell];
    
(*-----------------------------------*)
BeginPackage[
	"VQM`FastFourier`",
	"VQM`ArgColorPlot`"
	];
(*-----------------------------------*)

VQM`FastFourier`Private`Symbols = Hold[
	QFourierList, QInverseFourierList, QFourierListArgColorPlot,
	QInverseFourierListArgColorPlot, QStepSize, QFourierRange,
	QGrid, QLeftBorder, QRightBorder, QSpaceStep, QSpaceInterval,
	QIndexPosition, QFourierGrid, QFourierLeftBorder, QFourierRightBorder,
	QFourierStep, QFourierInterval];

Unprotect @@ VQM`FastFourier`Private`Symbols;
ClearAll @@ VQM`FastFourier`Private`Symbols;
 
QFourierList::usage = "QFourierList[xlist,{left,right}] returns the Fourier transform
	of the function f[x] represented by xlist. The result is given in the form {klist, {kmin, kmax}}.
	xlist is a list of complex numbers that can be interpreted as the discretization of a
	complex-valued function f[x] defined in the interval [left,right].
	xlist should consist of an even number of values. Use f /@ QGrid[a,b,n] to generate
	the values, and left=QLeftBorder[a,b,n], right = QRightBorder[a,b,n].
	klist=QFourierList[xlist,{left,right}][[1]] is a list of complex values in k-space.
	It can be interpreted as the discretization of the Fourier transform of f,
	defined on an appropriate interval in Fourier space (k-space).
	The corresponding interval in k-space is {kmin, kmax} == {-Pi/dx,Pi/dx},
	and the step size in k-space is given by
	dk = QFourierStep[left,right,n] == 2 Pi / (right-left).
	Package: VQM`FastFourier`.";

QInverseFourierList::usage = "QInverseFourierList[klist,{left,right}] returns
	the inverse Fourier transform of the function f[k] represented by klist.
	See QFourierList.
	We have QInverseFourierList[QFourierList[xlist,{left,right}]]=={xlist,{left,right}}.
	The size of the interval in Fourier space is related to the step size in x-space
	by dx = 2 Pi / (right-left).
	The step size dk in k-space determines the size of the x-interval according to
	{xleft, xright} == QInverseFourierList[klist,{left,right}][[2]] == {-Pi/dk, Pi/dk}.
	Package: VQM`FastFourier`.";

QFourierListArgColorPlot::usage = "QFourierListArgColorPlot[xlist,{left,right}] produces
	a QArgColorPlot of the QFourierList obtained from xlist and {left,right}.
	xlist can be interpreted as a discretization of a complex-valued function f[x] on a
	uniform grid of x-values in the interval [left,right].
	The QFourierList is a discretization of the FourierTransform of the function f[x].
	The option QFourierRange->{kmin,kmax} can be used to restrict the plot region in k-space.
	Package: VQM`FastFourier`.";

QInverseFourierListArgColorPlot::usage = "QInverseFourierListArgColorPlot[klist,{left,right}]
	produces a QArgColorPlot of the QInverseFourierList obtained from klist and {left,right}.
	Package: VQM`FastFourier`.";

QStepSize::usage = "QStepSize[list,{left,right}] returns (right-left)/(Length[list]-1),
	that is, the distance between the x-values if list represents values of a function f[x]
	on a uniformly spaced grid of x-values.
	Package: VQM`FastFourier`.";

QFourierRange::usage = "QFourierRange->Automatic is an option for QFourierListArgColorPlot.
	Setting QFourierRange->{kmin,kmax} restricts the k-space interval for the plot to [kmin,kmax].
	The maximal value of the k-space interval is determined from the step-size dx in
	x-space as [-Pi/dx,Pi/dx].
	Package: VQM`FastFourier`.";

QGrid::usage = "QGrid[a,b,n] generates a list of n points within the interval [a,b].
	Defines a space discretization that is useful for a fast Fourier transform.
	The distance between two adjacent points is dx = (b-a)/n = QSpaceStep[a,b,n].
	n has to be an even number.
	Package: VQM`FastFourier`.";

QLeftBorder::usage = "QLeftBorder[a,b,n] is the first element of QGrid[a,b,n].
	Package: VQM`FastFourier`.";

QRightBorder::usage = "QLeftBorder[a,b,n] is the first element of QGrid[a,b,n].
	Package: VQM`FastFourier`.";

QSpaceStep::usage = "QSpaceStep[a,b,n] gives the step size in the list QGrid[a,b,n].
	Package: VQM`FastFourier`.";

QSpaceInterval::usage = "QSpaceInterval[a,b,n] gives the interval
	{QGrid[a,b,n][[1]],QGrid[a,b,n]][[n]]}.
	Package: VQM`FastFourier`.";

QIndexPosition::usage = "QIndexPosition[x,a,b,n] gives the position of x in the
	list QGrid[a,b,n].
	Package: VQM`FastFourier`.";

QFourierGrid::usage = "QFourierGrid[a,b,n] generates a list of n points in Fourier space.
	Defines the domain in Fourier space that corresponds to QGrid[a,b,n] in position space.
	QFourierGrid[a,b,n] are the points where the values of QFourierList are defined.
	Also works for QInverseFourierList.
	Package: VQM`FastFourier`.";

QFourierLeftBorder::usage = "QFourierLeftBorder[a,b,n] determines the left border
	of the interval where the values of QFourierList are defined. This is the first
	element of QFourierGrid[a,b,n].
	Package: VQM`FastFourier`.";

QFourierRightBorder::usage = "QFourierRightBorder[a,b,n] determines the right border
	of the interval where the values of QFourierList are defined. This is the last
	element of QFourierGrid[a,b,n].
	Package: VQM`FastFourier`.";

QFourierStep::usage = "QFourierStep[a,b,n] is the distance between adjacent values
	of QFourierGrid[a,b,n]. This is the step size in Fourier space.
	Package: VQM`FastFourier`.";

QFourierInterval::usage = "QFourierInterval[a,b,n] determines the interval
	in Fourier space in which the values of QFourierList are defined.
	QFourierInterval == {QFourierLeftBorder,QFourierRightBorder}.
	Package: VQM`FastFourier`.";


(*-----------------------------------*)
Begin["`Private`"];
(*-----------------------------------*)

QSpaceStep[a_,b_,n_?EvenQ] := (b-a)/n

QLeftBorder[a_,b_,n_] := Module[{dx = (b-a)/n},
		If[!EvenQ[n],Message[FF::noteven]];
		Floor[a/dx] dx + dx]

QRightBorder[a_,b_,n_] := Module[{dx = (b-a)/n},
		If[!EvenQ[n],Message[FF::noteven]];
		Floor[a/dx] dx + n dx]

QSpaceInterval[a_,b_,n_] := Module[{},
		If[!EvenQ[n],Message[FF::noteven]];
		{QLeftBorder[a,b,n],QRightBorder[a,b,n]}]	

QGrid[a_,b_,n_] :=
	Module[{dx = (b-a)/n},
		If[!EvenQ[n],Message[FF::noteven]];
		Table[Floor[a/dx] dx + k dx, {k, 1, n}]]

QIndexPosition[x_,a_,b_,n_?EvenQ] :=
	Module[{xleft = QLeftBorder[a,b,n], dx = QSpaceStep[a,b,n]},
		If[!EvenQ[n],Message[FF::noteven]];
		Round[(x-xleft)/dx]+1]

QFourierStep[a_,b_,n_] := Module[{},
		If[!EvenQ[n],Message[FF::noteven]];
		2 Pi /(b-a)]

QFourierLeftBorder[a_,b_,n_] := Module[{dx = (b-a)/n},
		If[!EvenQ[n],Message[FF::noteven]];
		(2-n) Pi/(b-a)]

QFourierRightBorder[a_,b_,n_] := Module[{},
		If[!EvenQ[n],Message[FF::noteven]];
		n Pi/(b-a)]

QFourierGrid[a_,b_,n_]:=
	Module[{},
		If[!EvenQ[n],Message[FF::noteven]];
		Range[QFourierLeftBorder[a,b,n], QFourierRightBorder[a,b,n],
				QFourierStep[a,b,n]]
	]

QFourierInterval[a_,b_,n_] := Module[{},
		If[!EvenQ[n],Message[FF::noteven]];
		{QFourierLeftBorder[a,b,n],QFourierRightBorder[a,b,n]}]

FF::noteven = "Number of grid points (last argument) should be a positive
even integer. The result may not have the properties desired for use
together with QFourierList or QInverseFourierList.";

QStepSize[{xlist_,{xleft_,xright_}}] := QStepSize[xlist,{xleft,xright}];

QStepSize[xlist_,{xleft_,xright_}] := (xright-xleft)/(Length[xlist]-1);

QFourierList[{xlist_,{xleft_,xright_}}] := QFourierList[xlist,{xleft,xright}]

QFourierList[xlist_,{xleft_,xright_}]:=
   Module[{n,dx,a,b,klist,kleft,kright,dk},
      n = Length[xlist];
      dx = QStepSize[xlist,{xleft,xright}];
      a=-xleft/dx;
      klist =
      Sqrt[n/(2 Pi)] dx RotateLeft[
         InverseFourier[RotateLeft[xlist,Round[a]-1]],
      Round[n/2+1]];
      kleft = -Pi/dx + 2 Pi/(n dx); kright = Pi/dx;
      {klist,{kleft,kright}}
   ];

QInverseFourierList[{xlist_,{xleft_,xright_}}] := QInverseFourierList[xlist,{xleft,xright}]

QInverseFourierList[xlist_,{xleft_,xright_}]:=
   Module[{n,dx,a,b,klist,kleft,kright,dk},
      n = Length[xlist];
      dx = QStepSize[xlist,{xleft,xright}];
      a=-xleft/dx;
      klist =
      Sqrt[n/(2 Pi)] dx RotateLeft[
         Fourier[RotateLeft[xlist,Round[a]-1]],
      Round[n/2+1]];
      kleft = -Pi/dx + 2 Pi/(n dx); kright = Pi/dx;
      {klist,{kleft,kright}}
   ];

Options[QFourierListArgColorPlot] = Join[Options[QArgColorPlot],{QFourierRange->All}];

QFourierListArgColorPlot[{xlist_,{left_,right_}}, opts___?OptionQ] :=
	QFourierListArgColorPlot[xlist,{left,right}, opts];

QFourierListArgColorPlot[xlist_,{left_,right_}, opts___?OptionQ] :=
	Module[{klist,kleft,kright,dk,ran, kL, kR},
		{klist,{kleft,kright}} = QFourierList[xlist,{left,right}];
		dk = QStepSize[klist,{kleft,kright}];
		ran = QFourierRange/.Flatten[{opts, Options[QFourierListArgColorPlot]}];
		If[ran === All, {kL,kR}={kleft,kright}];
		If[MatchQ[N[ran],{_?NumberQ,_?NumberQ}], {kL,kR}=ran];
		QListArgColorPlot[klist, opts, QHorizontalRange->{{kleft,kright},{kL,kR}}]
];

Options[QInverseFourierListArgColorPlot] = Options[QFourierListArgColorPlot];

QInverseFourierListArgColorPlot[{xlist_,{left_,right_}}, opts___?OptionQ] :=
	QFourierListArgColorPlot[xlist,{left,right}, opts];

QInverseFourierListArgColorPlot[xlist_,{left_,right_}, opts___?OptionQ] :=
	Module[{klist,kleft,kright,dk,indk,partlist, kL, kR},
		{klist,{kleft,kright}} = QInverseFourierList[xlist,{left,right}];
		dk = QStepSize[klist,{kleft,kright}];
		ran = QFourierRange/.Flatten[{opts, Options[QFourierListArgColorPlot]}];
		If[ran === All, {kL,kR}={kleft,kright}];
		If[MatchQ[N[ran],{_?NumberQ,_?NumberQ}], {kL,kR}=ran];
		QListArgColorPlot[klist, opts, QHorizontalRange->{{kleft,kright},{kL,kR}}]
];

(*-----------------------------------*)
End[];      (* end `Private` Context *)
(*-----------------------------------*)

Protect @@ VQM`FastFourier`Private`Symbols;

(*-----------------------------------*)
EndPackage[]; (* end package Context *)
(*-----------------------------------*)

If[VQMmsgon, On[General::"spell"]];
If[VQMmsgon1, On[General::"spell1"]];

(* :Example: *)
(*

(* a complex-valued function (Gaussian): *)

f[x_,x0_,k0_,a_]:= (a/Pi)^(1/4) Exp[-a (x-x0)^2/2 + I k0 (x-x0)];

(* and its Fourier transform: *)

g[k_,x0_,k0_,a_]:= (1/a/Pi)^(1/4) Exp[-(1/a) (k-k0)^2/2 - I x0 k];

a=2; x0=1/2; k0=-2;
QArgColorPlot[f[x,x0,k0,a],{x,-3,4}, PlotRange->All, Axes->{True,False}, Frame->True];
QArgColorPlot[g[k,x0,k0,a],{k,-7,3}, PlotRange->All, Axes->{True,False}, Frame->True];

(* discretization *)

xlist = QGrid[-7,7,512];
flist = f[#,x0,k0,a]& /@ xlist;
xleft = QLeftBorder[-7,7,512]; xright = QRightBorder[-7,7,512];

(* Fourier transform of xlist. Compare with plot of g above: *)

QListArgColorPlot[QFourierList[flist,{xleft,xright}][[1]],
	PlotRange->All, Axes->{True,False}, Frame->True];

QFourierListArgColorPlot[flist,{xleft,xright},
	PlotRange->All, Axes->{True,False}, Frame->True];

QFourierListArgColorPlot[flist,{xleft,xright},
	PlotRange->All, Axes->{True,False}, Frame->True, QFourierRange->{-7,3}];
	(* Note: A larger region in x-space gives smaller step size in k-space *)

(* Inverse transformation: *)

klist = QFourierList[flist,{xleft,xright}];
newxlist=QInverseFourierList[klist][[1]];

Max[Abs[flist-newxlist]] (* identical within numerical accuracy *)

*)
