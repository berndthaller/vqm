(* :Title:	Free Time Evolution *)

(* :Name:	VQM`Free` *)

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
	Defines free Gaussian wave packets
*)

(* :Date:	2007-07-20 *)

(* :Package Version:		2.0 *)

(* :Mathematica Version:	6.0.1 *)

(* :Keywords:
	Free Time Evolution, Gaussian, Wave Packets
*)
    
VQMmsgon  = Head[General::"spell"]  =!= $Off;
VQMmsgon1 = Head[General::"spell1"] =!= $Off;

Off[General::spell1, General::spell];


(*-----------------------------------*)
BeginPackage["VQM`Free`"];              
(*-----------------------------------*)    

ClearAll[$QFreeMass,$QFreeSpaceDimension];

VQM`Free`Private`Symbols = Hold[
	QFreeHamiltonian1D, QFreeGaussian,
	QFreeGaussian1D, QFreeGaussian2D, QFreeGaussian3D,
	QGaussian1D, QGaussian2D, QGaussian3D,
	QFourierGaussian1D, QFreeFourierGaussian1D, QEnergyGaussian1D,
	QFreeFallGaussian1D,
	QFreeMass,QFreeSpaceDimension];

Unprotect @@ VQM`Free`Private`Symbols;
ClearAll @@ VQM`Free`Private`Symbols;

QFreeHamiltonian1D::usage = "QFreeHamiltonian1D[f[x],x, opts] applies
	the one-dimensional free Hamiltonian operator to the function
	f[x]. The option QFreeMass->m (default $QFreeMass=1) defines the mass
	of the particle (we use units with hbar = 1).
	Package: VQM`Free`.";

QFreeGaussian::usage = "QFreeGaussian[x,t,x0,p0,a,opts] is a solution of the
	free Schroedinger equation for a particle with mass m moving in
	n dimensions.
	In one dimension, the Schroedinger equation reads
	I D[psi[x,t],t] == -1/(2 m) D[psi[x,t],{x,2}].
	The mass can be defined by setting the option QFreeMass->m
	(default value is $QFreeMass = 1).
	The space dimension n (the dimension of x) may be specified by
	the option QFreeSpaceDimension->n.
	It can be n = 1, 2, or 3 (default value is $QFreeSpaceDimension=1).
	The arguments x0,p0, and a must have the same dimension.
	The initial wave packet has a width determined by a, an average
	initial position given by x0 and an average momentum given by
	p0. Package: VQM`Free`.";

QFreeGaussian1D::usage = "QFreeGaussian1D[x,t,x0,p0,a,m] is a solution of
	the one dimensional free Schroedinger equation for a particle with mass m.
	The particle has a Gaussian position and momentum
	distribution. The parameter a characterizes the width of the distribution,
	x0 is the average initial position, p0 is the average momentum.
	The arguments x0,p0,a,m may be omitted.
	Default values are x0=0, p0=0, a=1, m=1.
	The initial condition is QFreeGaussian1D[x,0,x0,p0,a,m]==QGaussian1D[x,x0,p0,a].
	Package: VQM`Free`.";

QFreeGaussian2D::usage = "QFreeGaussian2D[x,y,t,x0,p0,a,m] is a solution of
	the two dimensional free Schroedinger equation for a particle with mass m.
	Here x0,p0,a are two dimensional lists.
	The particle has Gaussian position and momentum distributions.
	The parameter a={a1,a2} characterizes the width of the distribution
	in the x and y direction,
	x0={x01,x02} is the average initial position, p0={p01,p02} is the average momentum.
	Only the argument m is optional, default value is m=1.
	The initial condition is QFreeGaussian2D[x,y,0,x0,p0,a,m]==QGaussian2D[x,y,x0,p0,a].
	Package: VQM`Free`.";

QFreeGaussian3D::usage = "QFreeGaussian3D[x,y,z,t,x0,p0,a,m] is a solution of
	the three dimensional free Schroedinger equation for a particle with mass m.
	The arguments x0,p0,a are three dimensional lists.
	The particle has Gaussian position and momentum distributions.
	The parameter a={a1,a2,a3} characterizes the width of the distribution
	in the x, y, and z direction,
	x0={x01,x02,x03} is the average initial position, p0={p01,p02,p03} is the average momentum.
	Only the argument m is optional, default value is m=1.
	The initial condition is QFreeGaussian3D[x,y,z,0,x0,p0,a,m]==QGaussian3D[x,y,z,x0,p0,a].
	Package: VQM`Free`.";

QGaussian1D::usage = "QGaussian1D[x,x0,p0,a] is a normalized Gaussian function in one dimension,
	centered at x0 in position space and at p0 in momentum space.
	The parameter a describes the width in position space.
	Package: VQM`Free`.";

QGaussian2D::usage = "QGaussian2D[x,y,x0,p0,a] is a normalized Gaussian function in two dimensions,
	centered at x0={x01,x02} in position space and at p0={p01,p02} in momentum space.
	The parameter a = {a1,a2} describes the width in position space.
	Package: VQM`Free`.";

QGaussian3D::usage = "QGaussian3D[x,y,z,x0,p0,a] is a normalized Gaussian function in three dimensions,
	centered at x0={x01,x02,x03} in position space and at p0={p01,p02,p03} in momentum space.
	The parameter a = {a1,a2,a3} describes the width in position space.
	Package: VQM`Free`.";

QFreeFourierGaussian1D::usage = "QFreeFourierGaussian1D[p,t,x0,p0,a,m] is the Fourier transform of
	the function QFreeGaussian1D[x,t,x0,p0,a,m].
	The arguments x0,p0,a,m are optional. Default values are x0=0, p0=0, a=1, m=1.
	Package: VQM`Free`.";

QFourierGaussian1D::usage = "QFourierGaussian1D[p,x0,p0,a] is the Fourier transform of
	the function QGaussian1D[x,x0,p0,a,m]. The parameter a describes the width in position
	space, that is, 1/a describes the width in Fourier space.
	Package: VQM`Free`.";

QEnergyGaussian1D::usage = "QEnergyGaussian1D[En,x0,p0,a] is a normalized
	free Gaussian function in the energy representation. The mass of the particle
	is m=1 and the space dimension is 1. x0 is the average initial position
	and p0 is the average momentum of the Gaussian. The parameter a describes
	the width of the Gaussian in position space.
	Package: VQM`Free`.";

QFreeFallGaussian1D::usage = "FreeFallGaussian[x, t, x0, p0, a, c] is a solution of the
	one-dimensional free Schroedinger equation for a particle with mass m=1 moving in a linear
	potential V[x] = c*x. This describes the 'free fall' in quantum mechanics, i.e., a uniformly
	accelerated particle with a Gaussian initial condition. The initial average position is x0,
	the initial average momentum is p0, and the parameter a describes the initial width
	of the Gaussian distribution.
	Package: VQM`Free`.";

QFreeMass::usage = "QFreeMass is an option for QFreeGaussian, QFreeHamiltonian1D. QFreeMass->m sets
	the mass of the particle to m.
	Package: VQM`Free`.";

QFreeSpaceDimension::usage = "QFreeSpaceDimension is an option for QFreeGaussian.
	QFreeSpaceDimension->n (n=1,2,3) sets the space dimension to n. Default is $QFreeSpaceDimension=1.
	Package: VQM`Free`.";

$QFreeMass::usage = "$QFreeMass is the default value for the mass ($QFreeMass = 1) in the package
	Free.m.
	Package: VQM`Free`.";

$QFreeSpaceDimension::usage = "$QFreeSpaceDimension is the default value for the option
	QFreeSpaceDimension. This is the space dimension used
	for the solution QFreeGaussian of the free Schroedinger equation.
	Package: VQM`Free`.";

VQMFree::nonneg := "The mass must be positive.Using QFreeMass = 1.";
VQMFree::baddim := "Dimension of x must be 1,2 or 3. Using QFreeSpaceDimension 1."

(*-----------------------------------*)
Begin["`Private`"];
(*-----------------------------------*)

FreeOptions := {QFreeSpaceDimension->$QFreeSpaceDimension,QFreeMass->$QFreeMass};

Options[FreeHamiltonian] = {QFreeMass->$QFreeMass};

QFreeHamiltonian1D[f_,x_Symbol, opts___?OptionQ] :=
	Module[{m},
		m = QFreeMass /. Join[{opts}, FreeOptions];
		func = Function[{x},f];
		-1/(2 m) D[func[x],{x,2}]
   ];

Options[QFreeGaussian] = FreeOptions;

QFreeGaussian[x_,t_,x0_,p0_,a_,opts___?OptionQ]:=
	Module[{dim,m},
    {dim,m} = {QFreeSpaceDimension,QFreeMass} /. {opts} /. FreeOptions;
	Which[dim==1,
	QFreeGaussian1D[x,t,x0,p0,a,m],
	dim==2,
	QFreeGaussian2D[x[[1]],x[[2]],t,
		{x0[[1]],x0[[2]]},{p0[[1]],p0[[2]]},{a[[1]],a[[2]]},m],
	dim==3,
	QFreeGaussian3D[x[[1]],x[[2]],x[[3]],t,
		{x0[[1]],x0[[2]],x0[[3]]},{p0[[1]],p0[[2]],p0[[3]]},{a[[1]],a[[2]],a[[3]]},m],
	True, 
	Message[VQMFree::baddim]]];

QFreeGaussian1D[x_,t_,x0_:0,p0_:0,a_:1,mass_:1]:=
	Sqrt[Sqrt[a/Pi]/(1+I t a/mass )]*
    Exp[-((a/2) (x-x0)^2-I p0 (x-x0)+I t p0^2/(2 mass))/(1+I t a/mass )];    

QFreeGaussian2D[x1_,x2_,t_,{x01_,x02_},{p01_,p02_},{a1_,a2_},m_:1]:=
	QFreeGaussian1D[x1,t,x01,p01,a1,m]*
	QFreeGaussian1D[x2,t,x02,p02,a2,m];

QFreeGaussian3D[x1_,x2_,x3_,t_,{x01_,x02_,x03},{p01_,p02_,p03_},{a1_,a2_,a3_},m_:1] :=
	QFreeGaussian1D[x1,t,x01,p01,a1,m]*
	QFreeGaussian1D[x2,t,x02,p02,a2,m]*
	QFreeGaussian1D[x3,t,x03,p03,a3,m];

(* normalized Gaussians (initial functions) *)

QGaussian1D[x_,x0_:0,p0_:0,a_:1] := 
	Sqrt[Sqrt[a/Pi]]*Exp[-((a/2) (x-x0)^2-I p0 (x-x0))];

QGaussian2D[x1_,x2_,{x01_,x02_},{p01_,p02_},{a1_,a2_}]:=
	QGaussian1D[x1,x01,p01,a1]*
	QGaussian1D[x2,x02,p02,a2];

QGaussian3D[x1_,x2_,x3_,{x01_,x02_,x03},{p01_,p02_,p03_},{a1_,a2_,a3_}] :=
	QGaussian1D[x1,x01,p01,a1]*
	QGaussian1D[x2,x02,p02,a2]*
	QGaussian1D[x3,x03,p03,a3];

(* Fourier transforms of Gaussian functions *)

QFreeFourierGaussian1D[p_,t_,x0_:0,p0_:0,a_:1,mass_:1] :=
	Sqrt[Sqrt[1/a/Pi]]*
		Exp[-(p-p0)^2/(2 a) - I p x0 - I t p^2/(2 mass)];

QFourierGaussian1D[p_,x0_,p0_,a_] :=
	Sqrt[Sqrt[1/a/Pi]]*
		Exp[-(p-p0)^2/(2 a) - I p x0];

(* Gaussian wave packet in the energy representation *)

k[En_] := Sqrt[2 En] /; En > 0

QEnergyGaussian1D[En_,x0_,p0_,a_] :=
	(1/Sqrt[k[En]])(1/(a Pi))^(1/4)*
		Exp[-(k[En]-p0)^2/(2 a)] Exp[-I x0 k[En]];



(* Uniformly accelerated Gaussian, free fall,
c = field strength in direction of x *)

QFreeFallGaussian1D[x_, t_, x0_, p0_, a_, c_] :=
  QFreeGaussian1D[x,t,x0 - c t^2/2,p0,a] Exp[-I c x t] Exp[1/6 (-I) c^2 t^3];


(*-----------------------------------*)
End[];      (* end `Private` Context *)
(*-----------------------------------*)

Protect @@ VQM`Free`Private`Symbols;

(* not protected symbols *)
$QFreeMass = 1;
$QFreeSpaceDimension = 1;

(*-----------------------------------*)
EndPackage[]; (* end package Context *)
(*-----------------------------------*)

If[VQMmsgon, On[General::"spell"]];
If[VQMmsgon1, On[General::"spell1"]];

