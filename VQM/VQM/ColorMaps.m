(* ::Package:: *)

(* :Title:   ColorMaps *)

(* :Name:    VQM`ColorMaps` *)

(* :Copyright: Copyright 2007 Bernd Thaller *)

(* :Author:  Bernd Thaller,
             Institute of Mathematics,
             University of Graz,
             A-8010 Graz
             bernd.thaller@uni-graz.at
*)

(* :Source:
	(Advanced) Visual Quantum Mechanics
	Springer-Verlag New York, 2004
	see "Book One" for a description of the color maps
*)

(* :Summary:
    This package defines color maps that associate color values with complex numbers.
    These color maps are used to visualize complex-valued functions.
*)

(* :Date:    Jun 20, 2007 *)

(* :Package Version:        2.0 *)

(* :Mathematica Version:    6.0.1 *)

(* :Keywords:
    Complex Numbers, Color Map
*)

VQMmsgon  = Head[General::"spell"]  =!= $Off;
VQMmsgon1 = Head[General::"spell1"] =!= $Off;

Off[General::spell1,General::spell];

(*-----------------------------------*)
BeginPackage[
    "VQM`ColorMaps`",    (* package Context *)
    "VQM`Spinors`",
    "Utilities`FilterOptions`"
    ];
(*-----------------------------------*)

Clear[$QComplexToColorMap];

$MaxAbsValue=1;

VQM`ColorMaps`Private`Symbols = Hold[
QComplexToColor, QRGBValues, QComplexToRGBValues, QComplexToRGBColor,
QComplexToColorMap, QValueRange, QLightnessRange, QSphereRadius,
QSaturationFromLightness, QBrightnessFromLightness, QHueFromArgument,
QLightnessFromModulus, QProcessColorMapOptions, QVectorToColor, QSpinorToColor];

Unprotect @@ VQM`ColorMaps`Private`Symbols;
ClearAll @@ VQM`ColorMaps`Private`Symbols;

QComplexToColor::usage =
"QComplexToColor[z,opts] associates a color to a complex number z.
The color map is given by the option QComplexToColorMap.
QComplexToColor[z] returns the result in the form
Hue[h,s,b] (in the HSB system). Package: VQM`ColorMaps`.";

QRGBValues::usage =
"QRGBValues[hue,lightness] converts color coordinates from the HLS system
to color coordinates in the RGB system. It associates the
RGB color values {r,g,b} to the given values of hue and lightness,
assuming that the saturation is equal to 1 (maximal saturation
at the given lightness). Package: VQM`ColorMaps`.";

QComplexToRGBValues::usage =
"ComplexToRGB[z] returns the color coordinates {r,g,b} in the
RGB system for a complex number z=x+I y. It uses the default
color map. Package: VQM`ColorMaps`.";

QComplexToRGBColor::usage =
"QComplexToRGBColor[z] returns the color RGBColor[r,g,b] of
the complex number z according to the standard color map. Package: VQM`ColorMaps`.";

QComplexToColorMap::usage =
"QComplexToColorMap is an option for QComplexToColor,
QComplexDensityPlot, and QComplexPlot3D. It specifies the color map
to be used when converting a complex number to a color.
The default setting is
QComplexToColorMap->$QComplexToColorMap, which determines
the color representing a complex number via a stereographic
projection from the complex plane onto the surface of the
color manifold in the HLS system. Package: VQM`ColorMaps`.";

QValueRange::usage =
"QValueRange is an option for QComplexToColor, QComplexDensityPlot,
and QComplexPlot3D. This option is used to determine parameters for the
default color map $QComplexToColorMap. QValueRange->{rmin,rmax},
where rmin < rmax are positive real numbers, causes every z with
Abs[z]<rmin (resp. Abs[z]>rmax) to be colored with minimal (resp.
maximal) lightness. If rmax < rmin, then Abs[z] < rmax has maximal
lightness, Abs[z] > rmin has minimal lightness. The default
Automatic corresponds to QValueRange->{0,Infinity}. The setting
for minimal/maximal lightness is determined by the
QLightnessRange option. Package: VQM`ColorMaps`.";

QLightnessRange::usage =
"QLightnessRange is an option for QComplexToColor, QComplexDensityPlot,
and QComplexPlot3D. QLightnessRange->{lmin,lmax}, where lmin and lmax
are real numbers in the interval [0,1], sets the minimal lightness
to lmin and the maximal lightness to lmax. Default is
QLightnessRange->{0,1}. The setting for this option is used to
determine parameters for the default color map $QComplexToColorMap.
Package: VQM`ColorMaps`.";

QSphereRadius::usage =
"QSphereRadius is an option for QComplexToColor, QComplexDensityPlot,
and QComplexPlot3D. QSphereRadius->R sets the radius of the sphere
used by $QComplexToColorMap to R. $QComplexToColorMap uses a
stereographic projection onto the surface of the color manifold
to determine the color of a complex number. Setting the radius
to R causes complex numbers with Abs[z]=R to be drawn at
lightness 1/2 (which corresponds to maximal brightness and
saturation in the HSB system). Package: VQM`ColorMaps`.";

$QComplexToColorMap::usage =
"$QComplexToColorMap[r,phi,{parameters}] associates a color to
a complex number. Default for QComplexToColor.
The complex number is
given in polar form, r=Abs[z], phi=Arg[z]. The result is returned
as $QComplexToColorMap[r,phi,{}] = Hue[h,s,b], in the HSB color
system. The Hue h of the color is given by phi/2/Pi,
the lightness is determined from r.
The color map can be described as a stereographic projection from
the complex plane onto the surface of the color manifold in the
Hue-Lightness-Saturation system. The optional parameters are
{R,bmin,bmax,lmin,lmax} specifying the radius R of the sphere, the
value bmin and bmax are those values of r which belong to the
minimal and maximal lightness lmin and lmax. Package: VQM`ColorMaps`.";

QSaturationFromLightness::usage = "QSaturationFromLightness[x] is an
auxiliary function provided by the package ColorMaps.
ComputesComputes the value of the saturation in the HSB-system from
the lighness x (0 <= x <= 1) of the color in the HLS-system,
assuming maximal HLS-saturation. Compiled for faster execution.
Package: VQM`ColorMaps`.";

QBrightnessFromLightness::usage = "QBrightnessFromLightness[x] is an
auxiliary function provided by the package ColorMaps.
Computes the value of the brightness in the HSB-system from
the lighness x (0 <= x <= 1) of the color in the HLS-system,
assuming maximal HLS-saturation. Compiled for faster execution.
Package: VQM`ColorMaps`.";

QHueFromArgument::usage = "QHueFromArgument[arg] is a compiled
auxiliary function provided by the
package ColorMaps. Computes the Hue from the argument arg
of a complex number. Package: VQM`ColorMaps`.";

QLightnessFromModulus::usage = "QLightnessFromModulus[parameters] defines a
compiled auxiliary function that depends on certain parameters. 
The parameters are given as a list of the form {R,bmin,bmax,lmin,lmax}.
See the description of $QComplexToColorMap for an explanation.
QLightnessFromModulus[parameters][r] computes the Lightness from the
modulus r of a complex number.
QLightnessFromModulus is provided by the package ColorMaps. Package: VQM`ColorMaps`.";

QProcessColorMapOptions::usage = "QProcessColorMapOptions[opts] is
an auxiliary function provided by the
package ColorMaps. Converts a list of options into a valid list of
parameters for the following color map functions: QLightnessFromModulus,
$QComplexToColorMap. Package: VQM`ColorMaps`.";

QVectorToColor::usage = "QVectorToColor[{u,v,w}] maps a three-dimensional
real vector to a unique color. The color is given as a list of three
real numbers (the RGB values of the color). The saturation depends on the
length of the vector and the hue is defined by the direction.
(red = positive x-direction, etc., the standard color circle in
the xy-plane. White = positive z-direction, black = negative z-direction,
50% gray = zero vector). Package: VQM`ColorMaps`.";

QSpinorToColor::usage = "QSpinorToColor[{u,v}] maps a spinor (i.e., a
C^2-vector) to a unique color. The color is given as a list of three
real numbers (the RGB values of the color).
This is done by first converting the spinor to a
vector (using QSpinorToVector from the package VQM`Spinors`) and
then applying QVectorToColor. Package: VQM`ColorMaps`.";

(*-----------------------------------*)
Begin["`Private`"];
(*-----------------------------------*)

myoptions =
	{QComplexToColorMap -> $QComplexToColorMap,
	QLightnessRange->Automatic,
	QSphereRadius->1., QValueRange->Automatic};


(* --------- QComplexToColor -------- *)

Options[QComplexToColor] = myoptions;

SetAttributes[QComplexToColor,Listable]

QComplexToColor[z_,opts___?OptionQ] :=
    Module[{colormap=QComplexToColorMap/.{opts}/.myoptions},
        If[colormap===Automatic ||
        colormap === $QComplexToColorMap,
        $QComplexToColorMap[Abs[z],Arg[z],QProcessColorMapOptions[opts]],
        colormap[Abs[z],Arg[z]]]
    ]


(* ------ $QComplexToColorMap ------- *)

$QComplexToColorMap[r_,arg_,
        {R_:1.,bmin_:0.,bmax_:1.,lmin_:0.,lmax_:1.}] :=
    Hue[0, 0, If[bmax<0,lmax,lmin]]/;r==0

$QComplexToColorMap[Infinity,arg_,
        {R_:1.,bmin_:0.,bmax_:1.,lmin_:0.,lmax_:1.}] :=
    Hue[
    	QHueFromArgument[arg],
    	QSaturationFromLightness[#],
    	QBrightnessFromLightness[#]
    ]&[If[bmax<0,lmin,lmax]]/;NumericQ[arg]

$QComplexToColorMap[Infinity,arg_,
        {R_:1.,bmin_:0.,bmax_:1.,lmin_:0.,lmax_:1.}] :=
    Hue[0, 0,If[bmax<0,lmin,lmax]]/;!NumericQ[arg]

$QComplexToColorMap[r_,arg_,
        {R_:1.,bmin_:0.,bmax_:1.,lmin_:0.,lmax_:1.}] :=
	Hue[
    	QHueFromArgument[arg], QSaturationFromLightness[#],
    	QBrightnessFromLightness[#]
    ]&[QLightnessFromModulus[{R,bmin,bmax,lmin,lmax}][r]]/;NumericQ[r]

$QComplexToColorMap[r_,arg_,
        {R_:1.,bmin_:0.,bmax_:1.,lmin_:0.,lmax_:1.}] :=
    Hue[0, 0,(lmin+lmax)/2]/;!NumericQ[r]

QLightnessFromModulus[{R_:1.,min_:0.,max_:1.,lmin_:0.,lmax_:1.}]:=
    Compile @@ 
     {{r},PiecewiseExpand[Min[{1,Max[{0,2/Pi ArcTan[r/R] (max-min) + min}]}](
            lmax-lmin)+lmin]};

(*
QLightnessFromModulus[{R_:1.,min_:0.,max_:1.,lmin_:0.,lmax_:1.}]:=
    Min[{1,Max[{0.,2/Pi ArcTan[r/R] (max-min) + min}]}](
            lmax-lmin)+lmin;
*)

QHueFromArgument = Compile[{arg}, If[arg < 0.0, 1.0 + arg/(2. Pi), arg/(2. Pi)]];

QSaturationFromLightness = Compile[{li}, If[li<=.5, 1, 2-2 li]];
QBrightnessFromLightness = Compile[{li}, If[li<=.5, 2 li, 1]];

(* simplified complex to color *)

QRGBValues[h_, l_] :=
	{
		{{2 l, 12 h l, 0}, {1, 2 l - 1 - 12 (l-1) h, 2 l - 1}},
		{{4 (1-3 h) l, 2 l, 0}, {3 + 12 (l-1) h - 2 l, 1, 2 l - 1}},
		{{0, 2 l, 4 (3 h-1) l}, {2 l - 1, 1, 6 l - 5 - 12 (l-1) h}},
		{{0, 4 (2-3 h) l, 2 l}, {2 l - 1, 7 + 12 (l-1) h - 6 l, 1}},
		{{4 (3 h-2) l, 0, 2 l}, {-9 - 12 (l-1) h + 10 l, 2 l - 1, 1}},
		{{2 l, 0, 12 (1-h) l}, {1, 2 l - 1, 11 + 12 (l-1) h - 10 l}}
	}[[1 + Floor[6 h], If[l<=1/2,1,2]]];

QComplexToRGBValues[z_] := QRGBValues[Mod[Arg[z]/(2 Pi),1],2/Pi ArcTan[Abs[z]]];

QComplexToRGBColor[z_] := RGBColor @@ QComplexToRGBValues[z];

(* determine params for $QComplexToColorMap from options *)

QProcessColorMapOptions[opts___] :=
    Module[ {rad,vran,lran,bds},
        {rad,vran,lran} =
            {QSphereRadius,QValueRange,QLightnessRange}
                /.{opts}/.myoptions;
        rad = testrad[rad];
        bds = borders[ testvran[vran]/rad];
        Flatten[{rad,bds,testlran[lran]}]
    ]

borders[{0.,Infinity}] = {0.,1.};
borders[{Infinity,0.}] = {1.,0.};
borders[vran_] :=
    ({#[[1]],#[[1]]+N[Pi]})/(#[[1]]-#[[2]])&[-2 ArcTan[vran]]


(* QSpinorToColor method *)

Arot = N[
	QRotationSO3[{0,1/Sqrt[2],-1/Sqrt[2]} ArcTan[1/Sqrt[2]]].
		QRotationSO3[{-1,0,0} Pi/4]
		];

invabs[u_] := If[u==0,1,1/(2 Abs[u])];
dist[{u_,v_,w_}]:= Min[invabs[u],invabs[v],invabs[w]];
fac[x_]:= 1-1/((3 x)^2 +1);

QVectorToColor[{a_,b_,c_}]:=
  Module[{v,absv,cvec,rotvec,v1},
    If[a==0 && b==0 && c==0, cvec={.5,.5,.5};Return[cvec] ];
    absv = Sqrt[a^2 + b^2 + c^2];
    v1={a,b,c}/absv;
    rotvec= Arot.v1;
    cvec = {.5,.5,.5} + rotvec (dist[rotvec] fac[absv] );
    cvec];

QSpinorToColor[{a_,b_}] := QVectorToColor[QSpinorToVector[{a,b}]];



(* -------- testing options --------- *)

testrad[rad_]:=
    Module[{},
        If[ !NumberQ[N[rad]] || !Positive[N[rad]],
            Message[ComplexPlot::badradius]; Return[1.],
            Return[N[rad]]
        ]
    ]

testvran[vran_]:=
    Module[{},
        If[vran===Automatic || vran=={0,Infinity},
            Return[{0.,Infinity}]];
        If[vran=={Infinity,0},Return[{Infinity,0.}]];
        If[ !VectorQ[vran] || Length[vran] != 2 || 
            Negative[N[vran[[1]]]] ||
            Negative[N[vran[[2]]]] ||
            N[vran[[1]]]==N[vran[[2]]],
            Message[ComplexPlot::badvrange]; {0.,Infinity},
            N[vran] ]
    ]

testlran[lran_]:=
    Module[{},
        If[lran===Automatic,Return[{0.,1.}]];
        If[ !VectorQ[lran]  || Length[lran] != 2 || 
            N[lran[[1]]]<0. || N[lran[[2]]]>1.   ||
            N[lran[[1]]]>1. || N[lran[[2]]]<0. ,
            Message[ComplexPlot::badlrange];{0.,1.},
            N[lran]]
    ]

(* ----------- messages ------------ *)

ColorMaps::badvrange =
"QValueRange should be Automatic or a list of two different
positive numbers (including 0 and Infinity).
Using Automatic instead.";

ColorMaps::badlrange =
"QLightnessRange should be Automatic or a list of two numbers
between 0 and 1. Using Automatic instead.";

ColorMaps::badradius =
"QSphereRadius must be a positive real number.
Using QSphereRadius->1. instead.";

(*-----------------------------------*)
End[];      (* end `Private` Context *)
(*-----------------------------------*)

Protect @@ VQM`ColorMaps`Private`Symbols;

(*-----------------------------------*)
EndPackage[]; (* end package Context *)
(*-----------------------------------*)

If[VQMmsgon, On[General::"spell"]];
If[VQMmsgon1, On[General::"spell1"]];


(* :Examples:

QComplexToColor[3. + 2. I]

QComplexToRGBColor[3. + 2. I]

<<VQM`ComplexPlot`

QColorArrayPlot[Table[QComplexToColor[x+I y],{y,-3,3,.5},{x,-3,3,.5}]]

QColorArrayPlot[
	Table[RGBColor @@ QVectorToColor[{x,y,0}],{y,-Pi,Pi,Pi/5},{x,-Pi,Pi,Pi/5}],
	MeshRange->{{-Pi,Pi},{-Pi,Pi}}]
  
QColorArrayPlot[
	Table[RGBColor @@ QVectorToColor[{x,0,z}],{z,-Pi,Pi,Pi/5},{x,-Pi,Pi,Pi/5}], 
	MeshRange->{{-Pi,Pi},{-Pi,Pi}}]

*)
