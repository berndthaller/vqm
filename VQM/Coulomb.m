(* :Title:   Coulomb *)

(* :Name:    VQM`Coulomb` *)

(* :Copyright: Copyright 2007 Bernd Thaller *)

(* :Author:  Bernd Thaller,
             Institute of Mathematics,
             University of Graz,
             A-8010 Graz
             bernd.thaller@uni-graz.at
*)

(* :Source:
	Advanced Visual Quantum Mechanics
	Springer-Verlag New York, 2004,
	in particular Chapter 2
*)

(* :Summary:
This package provides definitions for the quantum
mechanical Coulomb system. It gives the solutions
in cartesian as well as polar coordinates in two and
three dimensions. *)

(* :Date:    Jul 20, 2007 *)

(* :Package Version:        2.0 *)

(* :Mathematica Version:    6.0.1 *)

(* :Keywords:
    QCoulombFunction, Coulomb potential, quantum mechanics, Schroedinger equation
*)

VQMmsgon  = Head[General::"spell"]  =!= $Off;
VQMmsgon1 = Head[General::"spell1"] =!= $Off;

Off[General::spell1,General::spell];

(*-----------------------------------*)
BeginPackage[
	"VQM`Coulomb`",
	"VectorAnalysis`"
	];                  (* prevent shadowing *)
(*-----------------------------------*)

ClearAll[$QCoulombSpaceDimension,$QCoulombCoupling];

VQM`Coulomb`Private`Symbols=Hold[
	QPrincipalQuantumNumber,QRadialQuantumNumber,
	QCoulombSpaceDimension,QCoulombCoupling,QCoulombEnergy,
	QCoulombTimePeriod, QRadialCoulombFunction,
	QRadialPositionAmplitude, QCoulombFunction, QEffectiveCoulombPotential];

Unprotect @@ VQM`Coulomb`Private`Symbols; 
ClearAll @@ VQM`Coulomb`Private`Symbols;

QPrincipalQuantumNumber::usage = "QPrincipalQuantumNumber[nrad,ell] gives the principal
	quantum number in terms of the radial quantum number and the angular momentum.
	Package: VQM`Coulomb`.";

QRadialQuantumNumber::usage = "QRadialQuantumNumber[n,ell] is the radial quantum number
	n - ell - 1. It counts the number of radial zeros of the wave function. Package: VQM`Coulomb`.";

QCoulombEnergy::usage = "QCoulombEnergy[n, options] gives the energy of a particle
	in the Coulomb field g/r (here n=nrad+ell+1 is the principal quantum number).
	QCoulombEnergy is an eigenfunction of the differential operator
	-1/2 Delta + g/r.
	The Coulomb coupling constant g can be specified by giving
	the option QCoulombCoupling->g. Default is g=$QCoulombCoupling.
	The constant $QCoulombCoupling is initially set to 1,
	but can be redefined by the user.
	The energy depends on the space dimension.
	Default space dimension is 3;
	you can change this by the option QCoulombSpaceDimension->dim,
	or by redefining the constant $QCoulombSpaceDimension. Package: VQM`Coulomb`.";

QCoulombTimePeriod::usage = "QCoulombTimePeriod[n1,n2,..] gives the time period
	of a superposition of Coulomb eigenstates with principal quantum numbers
	n1,n2, etc. Package: VQM`Coulomb`.";
	
QCoulombFunction::usage = "QCoulombFunction[n,ell,m,{x,y,z},options]
	gives the energy eigenfunction of a particle in the three-dimensional attractive Coulomb field -g/r.
	n is the principal quantum number, ell orbital angular momentum, m is the magnetic quantum number.
	The function is an eigenfunction of (-1/2 Delta - g/r) Psi = QCoulombEnergy Psi.
	The Coulomb coupling constant g can be specified by giving the option QCoulombCoupling->g.
	Default is g=$QCoulombCoupling (attractive for positive values of g).
	The constant $QCoulombCoupling is initially set to 1, but can be redefined by the user.
	The energy depends on the space dimension.
	Default space dimension is 3; you can change this by the option QCoulombSpaceDimension->dim,
	or by redefining the constant $QCoulombSpaceDimension.
	The default coordinate system is Cartesian, as set by the package Calculus`VectorAnalysis`.
	The coordinatesystem can be changed by the command SetCoordinates[Spherical].
	In this case, the function has to be used in the form QCoulombFunction[n,ell,m,{r,theta,phi}].
	QCoulombFunction[n,m,{x,y},options] resp. QCoulombFunction[n,m,{r,phi},options]
	give the eigenfunction for the two-dimensional Coulomb problem in Cartesian resp. Spherical coordinates.
	Package: VQM`Coulomb`.";

QEffectiveCoulombPotential::usage = "QEffectiveCoulombPotential[ell, r] is the effective potential for the radial Coulomb equation.
	Consists of the Coulomb potential and the repulsive angular momentum barrier,
	which depends on the space dimension. Package: VQM`Coulomb`.";

QCoulombHamiltonian::usage = "QCoulombHamiltonian[psi[x,y,z],{x,y,z}] or QCoulombHamiltonian[psi[x,y],{x,y}]
	evaluates the action of the Hamiltonian operator -1/2 Delta - g/r on the wave function psi.
	Package: VQM`Coulomb`.";

QRadialCoulombFunction::usage = "QRadialCoulombFunction[n,ell,r] is the radial part
	of the eigenfunction of the Schroedinger equation in the angular momentum
	subspace described by ell. Here n is the principal quantum number.
	n - ell - 1 (=the radial quantum number) is the number of zeros of the radial Coulomb eigenfunction.
	Package: VQM`Coulomb`.";

QRadialPositionAmplitude::usage = "QRadialPositionAmplitude[n,ell,r] describes the
	amplitude of the radial position density. The square of the radial position amplitude
	gives the probability density for having a position at the distance r from the
	coordinate origin. Here n is the principal quantum number, ell is the orbital angular
	momentum quantum number. Package: VQM`Coulomb`.";

QCoulombSpaceDimension::usage = "QCoulombSpaceDimension->dim is an option used in the package
	VQM`Coulomb`. The formulas for eigenfunctions and energies
	depend on the space dimension. Default value
	is $QCoulombSpaceDimension = 3, but this constant can be redefined by the user. Package: VQM`Coulomb`.";

QCoulombCoupling::usage = "QCoulombCoupling->g is an option used in the package
	VQM`Coulomb`. It describes the strength of the Coulomb potential -g/r
	in the Schroedinger equation. Default value
	is $QCoulombCoupling = 1, but this constant can be redefined by the user. Package: VQM`Coulomb`.";

$QCoulombSpaceDimension::usage = "Constant describing the default value of the space dimension
	in the package VQM`Coulomb`. Initially set to 3. Package: VQM`Coulomb`.";

$QCoulombCoupling::usage = "Constant describing the default value of the strength of the
	Coulomb potential g/r for the package VQM`Coulomb`. Initial value is
	g=1, which describes an attractive Coulomb potential. Package: VQM`Coulomb`.";

Begin["`Private`"];

CoulombOptions := {QCoulombSpaceDimension->$QCoulombSpaceDimension,
QCoulombCoupling->$QCoulombCoupling};

QPrincipalQuantumNumber[nrad_,ell_] := nrad + Abs[ell] + 1;
QRadialQuantumNumber[np_,ell_] := np - Abs[ell] - 1;

QCoulombHamiltonian[f_,{x_Symbol,y_Symbol,z_Symbol}, opts___?OptionQ] :=
   Module[{ell, dim, g, func},
      {ell, dim, g} = check[1, opts];
      If[dim == 2, Message[Coulomb::spacedim];$QCoulombSpaceDimension=3];
      func = Function[{x,y,z},f];
      -1/2(D[func[x,y,z],{x,2}]+D[func[x,y,z],{y,2}]+D[func[x,y,z],{z,2}]) -
      g*func[x,y,z]/Sqrt[x^2+y^2+z^2]
   ];

QCoulombHamiltonian[f_,{x_Symbol,y_Symbol}, opts___?OptionQ] :=
   Module[{ell, dim, g, func},
      {ell, dim, g} = check[1, opts];
      If[dim == 3, Message[Coulomb::spacedim2];$QCoulombSpaceDimension=2];
      func = Function[{x,y},f];
      -1/2(D[func[x,y],{x,2}]+D[func[x,y],{y,2}]) -
      g*func[x,y,z]/Sqrt[x^2+y^2]
   ];

QEffectiveCoulombPotential[ell1_, r_, opts___?OptionQ] :=
   Module[{ell, dim, g},
      {ell, dim, g} = check[ell1, opts];
      -g/r + (1/2) ell(ell+dim-2)/r^2
   ];

QCoulombEnergy[N_, opts___?OptionQ] :=
   Module[{np, dim, g},
      {np, dim, g} = checkN[N, opts];
      If[g<=0,Return[]];
      Which[dim == 2,
               -g^2/(2(np-1/2)^2),
            dim == 3,
               -g^2/(2*np^2)
      ]
   ];

QCoulombTimePeriod[N__Integer,opts___?OptionQ] :=
	Module[{states, g},
		g = QCoulombCoupling/. {opts} /. CoulombOptions;
		states={N};
		2 Pi/g^2 LCM @@ (-g^2/QCoulombEnergy[#,opts]& /@ states)
	];

QRadialCoulombFunction[np_, ell1_, r_, opts___?OptionQ] :=
   Module[{nrad = np - ell1 - 1, ell, dim, g},
      {ell, dim, g} = check[ell1, opts];
      If[g<=0,Return[]];
      Which[dim == 2,
               normalization2D[nrad,Abs[ell],g]*
               Exp[-2*g*r/(1 + 2*Abs[ell] + 2*nrad)]*r^Abs[ell]*
               Hypergeometric1F1[-nrad, 1 + 2*Abs[ell], (4*g*r)/(1 + 2*Abs[ell] + 2*nrad)],
            dim == 3,
               normalization3D[nrad,ell,g]*r^ell*
               E^(-g r/(1 + nrad + ell))*
               Hypergeometric1F1[-nrad, 2 + 2*ell, 2 g r/(1 + nrad + ell)]
      ]
   ];

QRadialPositionAmplitude[np_, ell1_, r_, opts___?OptionQ] :=
   Module[{ell, dim, g},
      {ell, dim, g} = check[ell1, opts];
      If[g<=0,Return[]];
      Which[dim == 2,
               Sqrt[r]*QRadialCoulombFunction[np, Abs[ell], r,
                  QCoulombSpaceDimension->2, QCoulombCoupling->g],
            dim == 3,
               r*QRadialCoulombFunction[np, ell, r,
                  QCoulombSpaceDimension->3, QCoulombCoupling->g]
      ]
   ];
   
QCoulombFunction[np_, ell1_, m_, {c1_,c2_,c3_}, opts___?OptionQ] :=
   Module[{ell, dim, g},
      {ell, dim, g} = check[ell1, opts];
      If[g<=0,Return[]];
      If[dim == 2, Message[Coulomb::spacedim];$QCoulombSpaceDimension=3];
      If[!Integer[m] || Abs[m]>ell,Message[Coulomb::wrongm];Return[]];
      If[CoordinateSystem =!= Cartesian && CoordinateSystem =!= Spherical,
         Message[Coulomb::wrongcoords];Return[]];
      Which[
         CoordinateSystem===Cartesian,
            {ra,th,ph}=CoordinatesFromCartesian[{c1,c2,c3},Spherical];
            QRadialCoulombFunction[np, ell, ra, QCoulombSpaceDimension->3, QCoulombCoupling->g]*
            SphericalHarmonicY[ell,m,th,ph],
         CoordinateSystem===Spherical,
            QRadialCoulombFunction[np, ell, c1, QCoulombSpaceDimension->3, QCoulombCoupling->g]*
            SphericalHarmonicY[ell,m,c2,c3]
      ]
   ];

QCoulombFunction[np_, ell1_, {c1_,c2_}, opts___?OptionQ] :=
   Module[{ell, dim, g},
      {ell, dim, g} = check[ell1, opts];
      If[g<=0,Return[]];
      If[dim == 3, Message[Coulomb::spacedim2];$QCoulombSpaceDimension = 2];
      If[CoordinateSystem =!= Cartesian && CoordinateSystem =!= Spherical,
         Message[Coulomb::wrongcoords];Return[]];
      Which[
         CoordinateSystem===Cartesian,
            {r,th,ph}=CoordinatesFromCartesian[{c1,c2,0},Spherical];
            QRadialCoulombFunction[np, Abs[ell], r, QCoulombSpaceDimension->2, QCoulombCoupling->g]*
            Exp[ I ell ph]/Sqrt[2 Pi],
         CoordinateSystem===Spherical,
            QRadialCoulombFunction[np, Abs[ell], c1, QCoulombSpaceDimension->2, QCoulombCoupling->g]*
            Exp[ I ell c2]/Sqrt[2 Pi]
      ]
   ];
         
(* time dependent functions *)

QRadialCoulombFunction[np_, el_, t_, r_, opts___?OptionQ] :=
	QRadialCoulombFunction[np, el, r, opts] Exp[-I QCoulombEnergy[np, opts] t];

QRadialPositionAmplitude[np_, el_, t_, r_, opts___?OptionQ] :=
	QRadialPositionAmplitude[np, el, r, opts] Exp[-I QCoulombEnergy[np, opts] t];

QCoulombFunction[np_, el_, m_, t_, {c1_,c2_,c3_}, opts___?OptionQ] :=
	QCoulombFunction[np, el, m, {c1,c2,c3}, QCoulombSpaceDimension->3, opts]*
	Exp[-I QCoulombEnergy[np, QCoulombSpaceDimension->3, opts] t];

QCoulombFunction[np_, el_, t_, {c1_,c2_}, opts___?OptionQ] := 
	QCoulombFunction[np, el, {c1,c2}, QCoulombSpaceDimension->2, opts]*
	Exp[-I QCoulombEnergy[np, QCoulombSpaceDimension->2, opts] t];

(* auxiliary functions *)

normalization3D[n_,ell_,g_] := 1/(2 ell + 1)! *
   Sqrt[(n + 2 ell + 1)!/(2 n! (n + ell + 1))]*
   (2 g/(1 + n + ell))^(3/2 + ell);
   
normalization2D[n_,ell_,g_] :=
   ((2 Abs[ell])!)^(-1) Sqrt[(n+2 Abs[ell])!/(4 g n!)] *
      (4*g/(1 + 2*Abs[ell] + 2*n))^(3/2 + Abs[ell]);

checkN[N_,opts___] :=
   Module[{dim1,g,np},
      {dim1,g} = {QCoulombSpaceDimension, QCoulombCoupling} /. {opts} /. CoulombOptions;
      If[g<=0,Message[Coulomb::nonneg]];
      dim = If[dim1 != 2 && dim1 != 3,
               Message[Coulomb::baddim];3,
               dim1];
      np = If[N<1, Message[Coulomb::np];1,N,N];
      Return[{np,dim,g}]]

check[ell1_,opts___] :=
   Module[{dim1,g,ell},
      {dim1,g} = {QCoulombSpaceDimension, QCoulombCoupling} /. {opts} /. CoulombOptions;
      If[g<=0,Message[Coulomb::nonneg]];
      dim = If[dim1 != 2 && dim1 != 3,
               Message[Coulomb::baddim];3,
               dim1];
      ell = If[dim==3 && ell1<0, Message[Coulomb::ellsign];-ell1,ell1,ell1];
      Return[{ell,dim,g}]]

Coulomb::nonneg			:= "Only nonnegative values of the coupling constant
(corresponding to an attractive Coulomb potential)
can lead to eigenfunctions and eigenvalues.";
Coulomb::baddim			:= "Only integer dimensions 2 or 3 are supported, using 3.";
Coulomb::ellsign		:= "Using absolute value of angular momentum quantum number";
Coulomb::np				:= "Principal quantum number must be >=1"
Coulomb::wrongm			:= "In QCoulombFunction[n,ell,m,{x,y,z}],
the quantum number m must be an integer with -ell<=m<=ell."
Coulomb::wrongcoords	:= "The default coordinate system must be either Cartesian or Spherical. Use
SetCoordinates[Cartesian], etc., to set the default coordinate system."
Coulomb::spacedim    	:= "Attention: Default space dimension is set to 2,
but the number of arguments correspond to dimension 3.
Changing default dimension to $QCoulombSpaceDimension=3."
Coulomb::spacedim2   	:= "Attention: Default space dimension is set to 3,
but number of arguments correspond to dimension 2.
Changing default dimension to $QCoulombSpaceDimension=2."

(*-----------------------------------*)
End[];      (* end `Private` Context *)
(*-----------------------------------*)

Protect @@ VQM`Coulomb`Private`Symbols;

(* not protected symbols: *)
$QCoulombSpaceDimension = 3;
$QCoulombCoupling = 1;
$CoulombMass = 1;

(*-----------------------------------*)
EndPackage[]; (* end package Context *)
(*-----------------------------------*)

If[VQMmsgon, On[General::"spell"]];
If[VQMmsgon1, On[General::"spell1"]];

