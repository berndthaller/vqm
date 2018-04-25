(* ::Package:: *)

(* :Title:   ComplexPlot *)

(* :Name:    VQM`ComplexPlot` *)

(* :Copyright: Copyright 2010 Bernd Thaller *)

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
This package provides commands for visualizing complex-valued
functions by generating two-dimensional QColorDensityGraphics,
ContourGraphics and three dimensional SurfaceGraphics
of complex-valued functions with a color code for complex
numbers.*)

(* :Date:    Nov 12, 2010 *)
(* :Date:    Apr 25, 2018 *)

(* :Package Version:        3.6 *)
(* added option QComplexToValueMap *)
(* QListComplexPlot3D was broken *)
(* added QComplexContourPlot3D *)
(* repaired QListComplexContourPlot *)

(* :Mathematica Version:    11.3.0 *)

(* :Keywords:
    DensityPlot, Complex Function, Wavefunction
*)

(* :Limitations:
    This version needs Mathematica 7.
    'QComplexContourPlot' is just a shorthand for combining
    a QComplexDensityPlot with a ContourPlot of Abs[f].
	QComplexSurfacePlot has no Mesh-option.
*)

(* :Acknowledgement:
    Thanks to Manfred Liebmann and Rolf Mertig for some
    performance-improving suggestions.
*)
    
VQMmsgon  = Head[General::"spell"]  =!= $Off;
VQMmsgon1 = Head[General::"spell1"] =!= $Off;

Off[CompiledFunction::"cfsa"];
Off[CompiledFunction::"cfse"];
Off[CompiledFunction::"cfn"];
Off[CompiledFunction::"cfex"];
Off[General::"ovfl"];
Off[General::"unfl"];
Off[General::"munfl"];

Off[General::spell1,General::spell];

(*-----------------------------------*)
BeginPackage[
    "VQM`ComplexPlot`",			(* package Context *)
	"VQM`ColorMaps`",
    "VQM`ArgColorPlot`"		(* needed for QNiceTicks *)
    ];				
(*-----------------------------------*)

$MaxAbsValue=1;

VQM`ComplexPlot`Private`Symbols=Hold[
QComplexPlot3D,QListComplexPlot3D,
QComplexSurfacePlot,QListComplexSurfacePlot,QComplexDensityPlot,
QSpinorDensityPlot,QListComplexDensityPlot,QComplexContourPlot,
QListComplexContourPlot,QColorArrayPlot,QColorDensityGraphics,
QHighlights,QValueChecking,QScaledValues];

Unprotect @@ VQM`ComplexPlot`Private`Symbols;
ClearAll @@ VQM`ComplexPlot`Private`Symbols;


QComplexPlot3D::usage =
"QComplexPlot3D[f[x,y],{x,xmin,xmax},{y,ymin,ymax},opts]generates a
surface plot of a complex-valued function f of
two variables. The height of the surface is given by the
absolute value, the color is determined by the complex value
of f according to the option QComplexToColorMap.
Package: VQM`ComplexPlot`.";

QListComplexPlot3D::usage =
"QListComplexPlot3D[array,opts] generates a SurfaceGraphics
from a two-dimensional array of complex numbers. The height
of the surface is given by the absolute value, the color is
determined by the complex value of f according to the option
QComplexToColorMap.
Package: VQM`ComplexPlot`.";

QComplexSurfacePlot::usage =
"QComplexSurfacePlot[f[x,y],{x,xmin,xmax},{y,ymin,ymax},opts] is
similar to QListComplexPlot3D, but with a 'real surface look'.
The option QHighlights->On (default) lets the surface appear
glossy. Package: VQM`ComplexPlot`.";

QListComplexSurfacePlot::usage =
"QListComplexSurfacePlot[array,opts] is similar to QListComplexPlot3D,
but with a 'real surface look'.
The option QHighlights->On (default) lets the surface appear
glossy. Package: VQM`ComplexPlot`.";

QComplexDensityPlot::usage =
"QComplexDensityPlot[f[x,y],{x,xmin,xmax},{y,ymin,ymax},opts]
generates a QColorDensityGraphics of a complex-valued function f
of two variables x and y. It is similar to DensityPlot. The
complex value f[x,y] is mapped one-to-one to a color. The color
map is given by the option QComplexToColorMap. The default
$QComplexToColorMap determines the hue of the color from Arg[z],
and the lightness from Abs[z]. Package: VQM`ComplexPlot`.";

QSpinorDensityPlot::usage =
"QSpinorDensityPlot[{f,g}, {x,xmin,xmax}, {y,ymin,ymax}, opts]
visualizes a C^2-valued function in two dimensions by interlacing
colored density plots of the components f and g. Each 'pixel' thus
has an upper part with a color corresponding to the complex value
of the upper component f, and a lower part with a color corresponding
to the lower component g. Package: VQM`ComplexPlot`.";

QListComplexDensityPlot::usage =
"QListComplexDensityPlot[array,opts] gives a QColorDensityGraphics
of a two-dimensional array of complex numbers. It is similar to
ListDensityPlot. Each complex number is mapped one-to-one to a
color. The color map is determinded by the option
QComplexToColorMap. The default $QComplexToColorMap determines
the hue of the color from Arg[z], and the lightness from
Abs[z]. Package: VQM`ComplexPlot`.";

QComplexContourPlot::usage =
"QComplexContourPlot[f[x,y],{x,xmin,xmax},{y,ymin,ymax},opts]
visualizes a complex-valued function f of two variables
x and y. QComplexContourPlot combines a QColorDensityGraphics
with a ContourGraphics of the absolute value of f. Package: VQM`ComplexPlot`.";

QListComplexContourPlot::usage =
"QListComplexContourPlot[array,opts] generates a
QColorDensityGraphics of a two-dimensional array of complex
numbers and combines it with a ContourGraphics of Abs[array].
Package: VQM`ComplexPlot`.";

QComplexContourPlot3D::usage =
"QComplexContourPlot3D[f[x,y],{x,xmin,xmax},{y,ymin,ymax},opts]
visualizes a complex-valued function f of two variables
x and y. QComplexContourPlot3D combines QComplexPlot3D
with a ContourGraphics of the absolute value of f. The
contourlines appear as level lines on the phase-colored
three dimensional surface representing the absolute value
of the function. Package: VQM`ComplexPlot`.";

QColorArrayPlot::usage =
"QColorArrayPlot[list, opts] makes a 2D raster plot with colors
given by list. Here list is a two-dimensional array of color
directives. Package: VQM`ComplexPlot`.";

QColorDensityGraphics::usage =
"QColorDensityGraphics[absarray,colorarray,{opts}] is a
representation of a two-dimensional plot of an array of complex
numbers. It can be converted to SurfaceGraphics, ContourGraphics,
DensityGraphics and Graphics objects. Package: VQM`ComplexPlot`.";

QValueChecking::usage =
"QValueChecking is an option for QComplexDensityPlot and
QComplexContourPlot. It has no effect on the Plot3D commands.
The complex values 0 and ComplexInfinity have no well defined
color, in particular when QLightnessRange is {lmin,lmax} instead
of {0,1}. The default QValueChecking->On plots 0 as
GrayLevel[lmin], ComplexInfinity as GrayLevel[lmax] and
Indeterminate values in an intermediate gray level.
For slightly better performance in the case that there are
no exceptional points, use the setting QValueChecking->Off.
Package: VQM`ComplexPlot`.";

QScaledValues::usage = "QScaledValues is an option for QComplexDensityPlot
and QComplexContourPlot. It has no effect on the Plot3D commands.
If QScaledValues is set to True, then all values are scaled
so that the maximal absolute value is $MaxAbsValue (=1 by default).
Package: VQM`ComplexPlot`.";

QHighlights::usage =
"QHighlights is an option for QComplexSurfacePlot and QListComplexSurfacePlot.
QHighlights->On (default) lets the surface appear with Specularity[White, 15]. 
QHighlightes -> Specularity[ color, n ] can be used different specularities.
Package: VQM`ComplexPlot`.";

$QComplexToValueMap::usage = "$QComplexToValueMap=Abs.";

QComplexToValueMap::usage = "QComplexToValueMap for QComplexPlot3D and QListComplexPlot3D.
Specifies a transformation to be apply to the function in the argument of these plot
commands. Default is QComplexToValueMap -> Abs.
QComplexPlot3D[w[x,y],{x,x1,x2},{y,y1,y2},QComplexToValueMap -> Re] plots the surface
Re[func] colored with the usual color map for the complex values w[x,y].
QComplexPlot3D[w[x,y],{x,x1,x2},{y,y1,y2}, QComplexToValueMap->Re, QComplexToColorMap->(Hue[#1 Sin[#2]]&)],
Plots the real part of the function with a color derived from the imaginary part.";


(*-----------------------------------*)
Begin["`Private`"];
(*-----------------------------------*)

myoptions =
    {QComplexToColorMap->Automatic,
    QComplexToValueMap->Automatic,
        QLightnessRange->Automatic,
        QSphereRadius->1., QValueRange->Automatic, QScaledValues->False,
        QValueChecking->On, QHighlights->On};


$defaultPlotPoints = 25;

(* since option settings, e.g. for ListDensityPlot, changed considerably, 
   define a utitility function for how to migrate here:
*)

$adjust = {
  (Mesh -> None) -> (Mesh -> True),
  (PlotPoints -> Automatic) -> (PlotPoints -> $defaultPlotPoints),
  (PlotRange->{Full,Full,Automatic}) -> (PlotRange -> Automatic)
};

adjustOptions[ListDensityPlot] = Options[ListDensityPlot] /. $adjust;
adjustOptions[DensityPlot] = Options[DensityPlot] /. $adjust;
(*
Notice: Since Graphics is unhappy with the changed default PlotRange -> {Full,Full,Automatic} 
setting from DensityPlot the following will *not* work:
adjustOptions=Options;
*)


(* --------- QComplexPlot3D --------- *)

Options[QComplexPlot3D] =
    Join[myoptions,Options[Plot3D]];

$QComplexToValueMap = Abs;

SetAttributes[QComplexPlot3D,HoldAll]

QComplexPlot3D[func_,{x_Symbol,xmin_,xmax_},
        {y_Symbol,ymin_,ymax_}, opts___?OptionQ] := 
Module[{colfun, fli,cls1,cls,
        colormap=QComplexToColorMap/.Join[{opts}, myoptions],
        valuemap=QComplexToValueMap/.Join[{opts}, myoptions]},
       If[colormap === Automatic || colormap === $QComplexToColorMap,
          fli = QLightnessFromModulus[QProcessColorMapOptions[opts]];
          cls1 = Hue[QHueFromArgument[#2], QSaturationFromLightness[#1], QBrightnessFromLightness[#1]]&;
          cls = cls1[fli[#1],#2]&,
       (*else*) cls = colormap];
       If[valuemap === Automatic || valuemap === $QComplexToValueMap, valuemap=$QComplexToValueMap];
colfun  = Function @@ { {x,y,z}, cls[Abs[func],Arg[func]] };
Plot3D[valuemap[func],
              {x,xmin,xmax}, {y,ymin,ymax},
              ColorFunction ->  colfun,
              ColorFunctionScaling -> False,
              Evaluate[FilterRules[Flatten[{opts}],Options[Plot3D]]]
     ]
]/;And @@ NumericQ /@ {xmin,xmax,ymin,ymax};


(* --------- QComplexPlot3D --------- *)

Options[QComplexPlot3D] =
    Join[myoptions,Options[Plot3D]];

SetAttributes[QComplexPlot3D,HoldAll]

QComplexContourPlot3D[func_,{x_Symbol,xmin_,xmax_},
        {y_Symbol,ymin_,ymax_}, opts___?OptionQ] :=
	Module[{gr,cgr,absf,line3d,contStyle,linegraph, f=Function[{x,y},func],
        valuemap=QComplexToValueMap/.Join[{opts}, myoptions]},
       	If[valuemap === Automatic || valuemap === $QComplexToValueMap, valuemap=$QComplexToValueMap];
		gr = ContourPlot[valuemap[func],{x,xmin,xmax},{y,ymin,ymax},
				ContourShading->False,
				Evaluate[FilterRules[Flatten[{opts}], Options@ContourPlot]],
				PlotRange->All];
		cgr = Last/@ Graphics[gr][[1]];
		absf[n_] := valuemap[f[cgr[[n,1,3]][[1]],cgr[[n,1,3]][[2]]]];
		line3d[n_] := Map[Append[#,absf[n]+0.01]&, cgr[[n]],{2}];
		contStyle = ContourStyle/.Join[{opts},{ContourStyle->{GrayLevel[0.5], Thickness[0.0025]}}];
		contStyle=If[ListQ[contStyle],Sequence@@contStyle,contStyle];
		linegraph = {contStyle,Table[line3d[n],{n,1,Length[cgr]}]};
		Show[
			{QComplexPlot3D[func,
				{x, xmin, xmax}, {y, ymin, ymax},
				opts,
				Mesh -> False,
				QLightnessRange -> {0,1}],
			Graphics3D[linegraph]}
		]
]/;And @@ NumericQ /@ {xmin,xmax,ymin,ymax};

(* ------- QListComplexPlot3D ------- *)

Options[QListComplexPlot3D] =
    Join[myoptions, Options[ListPlot3D]];
    
QListComplexPlot3D[array_, opts___?OptionQ] := 
Module[{ab1=$QComplexToValueMap[array],absarr,plotarr,
        argarr=N[Arg[Drop[#,-1]& /@ Drop[array,-1]]],
        abss,ligarr,hues,stns,brts,colors,
        checkopts={QComplexToColorMap,QValueChecking}/.Join[{opts}, myoptions],
        scaledvl=QScaledValues/.Join[{opts}, myoptions],
        valuemap=QComplexToValueMap/.Join[{opts}, myoptions]},
       If[valuemap === Automatic || valuemap === $QComplexToValueMap,
             plotarr=ab1,
             (*else*) plotarr=Map[valuemap,array] ];
       If[scaledvl === True,absarr = $MaxAbsValue ab1/Max[ab1],absarr = ab1];
       Remove[ab1];
       abss=Drop[#,-1]& /@ Drop[absarr,-1];
       If[checkopts[[1]] === Automatic ||
          checkopts[[1]] === $QComplexToColorMap,
             ligarr = Map[QLightnessFromModulus[params], absarr, {2}];
             hues   = Map[QHueFromArgument, argarr, {2}];
             Remove[argarr];
	         stns   = Map[QSaturationFromLightness, ligarr, {2}];
	         brts   = Map[QBrightnessFromLightness, ligarr, {2}];
              Remove[ligarr];
	         colors = MapThread[Hue, {hues,stns,brts}, 2];
              Remove[hues,stns,brts];
	         If[checkopts[[2]]===On || checkopts[[2]] === True,
	         colors = checkvalues[arr,colors,params],
	         Remove[arr]],
          (*else*) colors =
            MapThread[checkopts[[1]][#1,#2]&,{abss,argarr},2];
            Remove[argarr]];
            ListPlot3D[ plotarr, colors,
            Evaluate[FilterRules[Flatten[{opts}], Options @ ListPlot3D]]
       ]
]/; MatrixQ[array]

(* ------- QComplexSurfacePlot ------- *)

Options[QComplexSurfacePlot] =
    Join[myoptions, {Compiled->True}, Options[Plot3D]];

SetAttributes[QComplexSurfacePlot,HoldAll]

QComplexSurfacePlot[func_,{x_Symbol,xmin_,xmax_},
        {y_Symbol,ymin_,ymax_}, opts___?OptionQ] := 
    Module[{meshchange, pl, comp, dx, f=Function[{x,y},func],
            bmax=N[{xmax,ymax}],bmin=N[{xmin,ymin}], gc, array},
        {pl,comp} = {PlotPoints,Compiled}/.
                     Join[{opts}, Options[QComplexSurfacePlot], adjustOptions[DensityPlot]];
        pl = testplotpoints[pl];
        dx  = (bmax-bmin)/(pl-1.);
        If[comp == True,
           gc = fcomp[f,bmin,dx];
           array =
              Table[gc[j,i],{j,0,pl[[2]]-1},{i,0,pl[[1]]-1}],
           array = Table[f[x,y],
                {y, N[ymin], ymax, dx[[2]]},
                {x, N[xmin], xmax, dx[[1]]}] ];
        QListComplexSurfacePlot[array,opts,
        	Ticks->{QNiceTicks[xmin,xmax,dx[[1]],6],
        		QNiceTicks[ymin,ymax,dx[[2]],6],Automatic}]
    ]/; And @@ NumericQ /@ {xmin,xmax,ymin,ymax}

(* ------- QListComplexSurfacePlot ------- *)

Options[QListComplexSurfacePlot] =
    Join[myoptions, Options[ListPlot3D]]/.(Mesh->Automatic) :> (Mesh->False);

$defaultSpecularity = Specularity[White,15];

QListComplexSurfacePlot[array_?MatrixQ, opts___?OptionQ] := 
Module[{interp, colorfun,
        clgs,hues,stns,brts,colors,
        hlite = QHighlights /.Join[{opts}, Options[QListComplexSurfacePlot]]},
       If[hlite === On, hlite = $defaultSpecularity];
       If[Head[hlite] =!= Specularity, hlite = Sequence[]];
       clgs = With[{qops=QProcessColorMapOptions[opts]},
					Function[z, QSaturationFromLightness[QLightnessFromModulus[qops][z]]]];

       brts = With[{qops=QProcessColorMapOptions[opts]},
					Function[z, QBrightnessFromLightness[QLightnessFromModulus[qops][z]]]];

(* RM: create a complex-valued interpolation function *)
       interp   =Interpolation[Flatten[Table[{i,j,array[[i,j]]},{i,Length[array]},{j,Length[array[[1]]]}],1]];
       colorfun[x_,y_,z_] := Block[{f=interp[x,y],ab},
               ab = Abs[f];
               Directive[Hue[Mod[QHueFromArgument[Arg[f]],1], clgs@ab, brts@ab],
                         hlite 
                        ]         ];
meshchange[ops___][z__] := If[FreeQ[{ops,z}, Mesh], 
                              Hold[Sequence][z, Mesh -> False], 
                              Hold[Sequence][z]
				             ] // ReleaseHold;

          ListPlot3D[Abs[array],ColorFunction -> colorfun,
            ColorFunctionScaling -> False,
            Evaluate[(meshchange[opts]@FilterRules[Flatten[{opts}], Options@ListPlot3D])],
            BoxRatios->{1,1,.4}
       ]
];


(* ------- QComplexDensityPlot ------- *)

Options[QComplexDensityPlot] =
    Sort @ Join[myoptions, adjustOptions[DensityPlot],
    {Compiled -> (*RM2018 True*)False}];

SetAttributes[QComplexDensityPlot,HoldAll]

QComplexDensityPlot[func_,
        {x_Symbol,xmin_,xmax_}, {y_Symbol,ymin_,ymax_},
        opts___?OptionQ] :=
    Module[{pl, comp, ops, dx, f=Function[{x,y},func],
            bmax=N[{xmax,ymax}],bmin=N[{xmin,ymin}], gc, array},  
        {pl,comp} = {PlotPoints,Compiled}/.
                     Join[{opts}, Options[QComplexDensityPlot],adjustOptions[DensityPlot]];     
 (*RM: the default setting of PlotPoints is not "Automatic" ... *) 
       If[pl === Automatic, pl = $defaultPlotPoints];
        pl = testplotpoints[pl];
        ops = JoinOptions[DataRange->{{xmin,xmax}, {ymin,ymax}},opts];

        dx  = (bmax-bmin)/(pl-1.);
        If[comp == True,
           gc = fcomp[f,bmin,dx];
           array =
              Table[gc[j,i],{j,0,pl[[2]]-1},{i,0,pl[[1]]-1}],
       array = Table[f[x,y]
(*RM:
Commenting this out 2018 ... 
TimeConstrained[ TimeConstrained[f[x,y], 1] /. 
                 {$Aborted:>10.^10000, Overflow[]:>10.^10000, Underflow[]:> -10.^10000}, 1
      ]/.$Aborted:>-10.^10000
*)
                 ,
                {y, N@(ymin+dx[[2]]/112.), ymax, dx[[2]]},
                {x, N@xmin, xmax, dx[[1]]}] ];
 QListComplexDensityPlot[array,ops]
    ]/; And @@ NumericQ /@ {xmin,xmax,ymin,ymax}

Options[QSpinorDensityPlot] = {Compiled -> True};
QSpinorDensityPlot[{func1_,func2_},{x_Symbol,xmin_,xmax_},{y_Symbol,
		ymin_,ymax_},
		opts___?OptionQ] :=
	Module[{pl,comp,ops,dx,f1=Function[{x,y},func1],f2=Function[{x,y},func2],
          bmax=N[{xmax,ymax}],bmin=N[{xmin,ymin}],gc1, gc2, array},
        {pl,comp} = {PlotPoints,Compiled}/.
        Join[{opts}, Options[QSpinorDensityPlot], adjustOptions[DensityPlot]];
        pl = testplotpoints[pl];
        ops = JoinOptions[DataRange->{{xmin,xmax},{ymin,ymax}},opts];
        dx  = (bmax-bmin)/(pl-1.);
        If[comp == True,
           gc1 = fcomp[f1,bmin,dx];gc2 = fcomp[f2,bmin,dx];
        array=Flatten[ 
            Table[{ Table[gc2[j,i],{i,0,pl[[1]]-1}],
                	Table[gc1[j,i],{i,0,pl[[1]]-1}] },
			{j,0,pl[[2]]-1}],1],
        array=Flatten[ 
			Table[{ Table[f2[x,y],{x,N[xmin],xmax,dx[[1]]}],
					Table[f1[x,y],{x,N[xmin],xmax,dx[[1]]}]},
			{y,N[ymin], ymax, dx[[2]]}],1]
		];
		QListComplexDensityPlot[array, ops, SpinorMesh->True]
	]/; And @@ NumericQ /@ {xmin,xmax,ymin,ymax}

fcomp[f_, {xmin_, ymin_}, {dx_, dy_}] := 
	Compile @@ {{{j,_Integer},{i,_Integer}},
	f[(xmin+i*dx),(ymin+j*dy)]}


(* ---- QListComplexDensityPlot ---- *)

Options[QListComplexDensityPlot] =
    Join[myoptions, adjustOptions[ListDensityPlot]];
    
QListComplexDensityPlot[array_?MatrixQ, opts___?OptionQ] :=
(
    Module[{arr = N[array],ab1,absarr,argarr,ligarr,hues,stns,brts,
            colors,params=QProcessColorMapOptions[opts],
            checkopts= {QComplexToColorMap,QValueChecking}/.Join[{opts}, myoptions],
            scaledvl=QScaledValues/.Join[{opts}, myoptions]
           },
       ab1 = Abs[arr];
       If[scaledvl === True, absarr = $MaxAbsValue ab1/Max[ab1], absarr = ab1];
       absarr =  Developer`ToPackedArray[ absarr /.Overflow[]:>10.^10000/.Indeterminate:>10.^10000];
       argarr = N[Arg[arr] /. Overflow[] :> 10.^10000];
       If[!MatchQ[argarr, {{__?NumberQ}..}],  argarr = Map[Replace[#, a_ /; !NumberQ[a] -> 0.]&, argarr, {2}]];
       If[!(FreeQ[ argarr, Indeterminate] ||  Max[argarr] < 10.^10000), 
          argarr = argarr/.{r_Real:>10.^10000 /;r>10.^10000, r_Real:>-10.^10000 /;r<-10.^10000, Indeterminate:>10.^10000};
         ];
       argarr = Developer`ToPackedArray[argarr // N];
       If[checkopts[[1]] === Automatic ||
          checkopts[[1]] === $QComplexToColorMap,
             ligarr = Map[QLightnessFromModulus[params], absarr, {2}];
             hues   =  Developer`ToPackedArray @ N[ Mod[Map[QHueFromArgument, argarr, {2}],1] ];
Remove[argarr];
	         stns   = Map[QSaturationFromLightness, ligarr, {2}];(*/.
                          {r_Real :> 0. /; r < 0. , r_Real :> 1 /; r > 1. } *)
	         brts   = Map[QBrightnessFromLightness, ligarr, {2}]; (*/.
                          {r_Real :> 0. /; r < 0. , r_Real :> 1 /; r > 1. } *)
Remove[ligarr];

	         colors = MapThread[List, {hues,stns,brts}, 2];
Remove[hues,stns,brts];
	        
	         If[checkopts[[2]]===On || checkopts[[2]] === True,
	         colors = checkvalues[arr,colors,params]; Remove[arr],
Remove[arr]
                   ],
          (*else*) colors =
            MapThread[checkopts[[1]][#1,#2]&,{absarr,argarr},2];
Remove[argarr];
         ];
If[!FreeQ[colors, Indeterminate], colors = colors /. Indeterminate:>0.];

If[Max[absarr] > 10.^10000 ,
   absarr = absarr/.r_Real:>10.^10000 /;r>10.^10000;
  ];

      QColorDensityGraphics[ absarr, colors,  
   {JoinOptions[{opts},Options[QColorDensityGraphics]]}]
    ]
);


SetAttributes[checkvalues, HoldAll];
checkvalues[arr_,colors_,params_]:=
   Module[{colors1=colors,bmax=params[[3]],
           lmin=params[[4]],lmax=params[[5]], hue=Head[colors[[1,1]]]},
      {li,la}=If[bmax<=0,{lmax,lmin},{lmin,lmax}] // N;
      colors1=ReplacePart[colors1, hue[0.,0.,li],
                      Position[arr,z_/;z==0]];
      colors1=ReplacePart[colors1, hue[0.,0.,.5 (lmin+lmax)],
                      Position[arr,z_/;z==Indeterminate]];
      ReplacePart[colors1,hue[0.,0.,la],
                      Position[arr,z_/;z==ComplexInfinity]]

   ];

(* ------ QColorDensityGraphics ------ *)

Options[QColorDensityGraphics] =
    Join[adjustOptions[ListDensityPlot],{ SpinorMesh->False}];

QColorDensityGraphics[abs_,colors_,opts___?OptionQ] :=
    Module[{ops=FilterRules[ {opts}, Options[QColorArrayPlot] ]},
        QColorArrayPlot[colors,ops]
    ]

QColorDensityGraphics[SurfaceGraphics[abs_,colors_,opts___],opts2___] :=
    QColorDensityGraphics[abs,colors,
    FilterRules[{JoinOptions[{opts2},opts]}, Options@QColorDensityGraphics ]]

QColorDensityGraphics/:
SurfaceGraphics[QColorDensityGraphics[abs_,colors_,opts___],opts2___] :=
    Module[{colorarray},
        If[Dimensions[abs]==Dimensions[colors],
            colorarray = Drop[#,-1]& /@ Drop[colors,-1],
            colorarray = colors];
        SurfaceGraphics[abs, colorarray,
         FilterRules[{JoinOptions[{opts2},opts]}, Options@SurfaceGraphics]
        ]
    ]

QColorDensityGraphics/:
DensityGraphics[QColorDensityGraphics[abs_,colors_,opts___],opts2___] :=
    DensityGraphics[abs,
         FilterRules[{JoinOptions[{opts2},opts]}, Options@DensityGraphics]
        ]

QColorDensityGraphics/:
ContourGraphics[QColorDensityGraphics[abs_,colors_,opts___],opts2___] :=
    ContourGraphics[abs,
        FilterRules[{JoinOptions[{opts2},opts]}, Options@ContourGraphics]
        ]

QColorDensityGraphics/:
Graphics[QColorDensityGraphics[abs_,colors_,opts___]] :=
    ColorArrayGraphics[colors,opts]


(* -------- QColorArrayPlot --------- *)
$huescale=256;

Options[QColorArrayPlot] = Options[QColorDensityGraphics];

QColorArrayPlot[colors_,opts___?OptionQ] :=
    ColorArrayGraphics[colors,opts]

ColorArrayGraphics[colorsin_,  opts___?OptionQ ]:=
	ColorArrayGraphics[colorsin, Hue, opts];

ColorArrayGraphics[colorsin_, colorfunhead_Symbol, opts___?OptionQ ]:=
    Module[{newcolorfunctionflag,colfunhead=colorfunhead,datran,msh,mshsty,
            dims=Dimensions[colorsin],colors=colorsin, fac,
            ops=JoinOptions[Flatten[{opts}],Options[QColorArrayPlot]]},
        {datran,msh,mshsty,mshtype}={DataRange,Mesh,MeshStyle,SpinorMesh}/.{ops};
        If[mshtype === True, fac = 2, fac = 1];
        If[ datran === Automatic,
            datran = { {0,dims[[2]]}, {0,dims[[1]]} } ];
        ops=FilterRules[{ops}, Options@Graphics];
        If[ mshsty === Automatic, mshsty = GrayLevel[0] ];

(* this is the default *)
If[!MatchQ[colors, {{{__Real}..}..}], 
   colfunhead = Head[colors[[1,1]]] ;
   If[colfunhead === Hue, newcolorfunctionflag = True; 
    colors = colors/.Hue[a_,b___] -> {Mod[a,1],b},
      newcolorfunctionflag = False
     ];
   , 
   (* else *)
   If[Max[colors]>1 || Min[colors]<-1,
      colors = colors/. {a_Real,b__}:> {Mod[a,1],b}
     ];
   newcolorfunctionflag = True
  ];

If[!FreeQ[colors, Overflow[]], colors = colors /. Overflow[] :> 1];

        If[msh === True || msh === Automatic,
           (* RM: This is kind of a hack, but things look really strange otherwise.
                   When Antialiasing is implemented in a better way: fix this 
           *)
        Function[ras, Style[Graphics[{ras,
                                      Flatten[ List @@ meshlines[datran,{dims[[2]],dims[[1]]/fac},mshsty]]
                                     },  ops
                                    ], Antialiasing -> False
                           ] 
                ],
        (* no meshlines : *)
        Function[ras, Graphics[ras, ops]]
          ] @
(* this gives Raster[  ] *)
If[TrueQ[newcolorfunctionflag], 
   colors = Round[$huescale colors];
   Raster[colors, datran // Transpose, {1,$huescale}, ColorFunction -> colfunhead]
   ,
   (*else *)
   (* keep the old setup for compatibility *)
   If[colfunhead === Hue,
      Raster[colors /. Hue[a_,b___] :> N[{Mod[a,1],b}],datran//Transpose,ColorFunction -> Hue ],
      Raster[colors /. colfunhead[a_,b___] :> N[{a,b}],datran//Transpose,ColorFunction->colfunhead]
     ]
  ]
   ]


(* ------- QComplexContourPlot ------- *)

Options[QComplexContourPlot]=
    Join[myoptions, {Antialiasing -> False},
    Options[ContourPlot]]/. (ColorFunction -> _) :>
    (ColorFunction -> (Hue[Mod[#,1],.1,1,0]&) );

SetAttributes[QComplexContourPlot,HoldAll];

QComplexContourPlot[func_,
        {x_Symbol,xmin_,xmax_}, {y_Symbol,ymin_,ymax_},
        opts___?OptionQ]:=
Module[{gr1,gr2,antial},
    If[TrueQ[Antialiasing /. Flatten[{opts}] /. Options[QComplexContourPlot]],
       antial = Identity,
       antial = Style[#, Antialiasing -> False]&
      ];
    gr1=QComplexDensityPlot[func,{x,xmin,xmax},{y,ymin,ymax},
        Mesh->False, PlotRange->All, opts];
    gr2=ContourPlot @@ {Abs[func], {x,xmin,xmax}, {y,ymin,ymax},
        ContourShading->False,
        FilterRules[Flatten[{opts}], Options[ContourPlot]] };
   antial @ Show[gr1,gr2,PlotRange -> All, FilterRules[{opts}, Options@Graphics]]
]/; And @@ NumericQ /@ {xmin,xmax,ymin,ymax}


(* ----- QListComplexContourPlot ----- *)

Options[QListComplexContourPlot]=
    Join[myoptions,Options[ListContourPlot]];

QListComplexContourPlot[array_,opts___?OptionQ]:=
	Module[{gr1,gr2,dims=Dimensions[array],datran,epi,ep,
		ops=JoinOptions[Flatten[{opts}],Options[QListComplexContourPlot]]},
        {datran,epi}={DataRange,Epilog}/.{ops};
        If[ datran === Automatic,
            datran = { {0,dims[[2]]}, {0,dims[[1]]} } ];
    gr2=ListContourPlot[Abs[array],ContourShading->False,
        DataRange->datran,
        FilterRules[Flatten[{opts}], Options[ListContourPlot]]
        ][[1]];
        ep={epi,gr2};
   QListComplexDensityPlot[array,
        Mesh->False,PlotRange->All,Epilog->ep,opts]
    ]/; MatrixQ[N[array]]

(* -------- testing options --------- *)

testplotpoints[pl_] :=
    Module[{},
    	If[ IntegerQ[pl] && Positive[pl], {pl,pl}, pl]
    	] 


(* ------ dealing with options ------ *)

(* joining two lists of rules with list1 having precedence
over list2: *)

(* RM: this will e.g. allow QListComplexDensityPlot to still accept MeshRange *)
$fixfor6 = {System`MeshRange -> DataRange};

JoinOptions[list1_,list2___]:=
    Module[{namelist=First /@ Flatten[{list1}/. $fixfor6 ]},
        Sequence @@ 
        Sort[Join[Flatten[{list1}/. $fixfor6 ],
            Select[Flatten[{list2}/. $fixfor6 ],
                !MemberQ[namelist,First[#]]&]]
        ]
    ] 


(* ------ generate mesh lines ------ *)
Clear[meshlines];
meshlines[m_,dms_,style_] :=
(* RM: this *should* work, but the FE discards Style ... *)    
(*Style[*)Graphics[
        Flatten[{{style},
        Table[Line[{{i,m[[2,1]]},{i,m[[2,2]]}}],
            {i,m[[1,1]],m[[1,2]],(m[[1,2]]-m[[1,1]])/dms[[1]]}],
        Table[Line[{{m[[1,1]],i},{m[[1,2]],i}}],
            {i,m[[2,1]],m[[2,2]],(m[[2,2]]-m[[2,1]])/dms[[2]]}]
        }]
    ](*, Antialiasing -> False]*)


(* ----------- messages ------------ *)

ComplexPlot::badvrange =
"QValueRange should be Automatic or a list of two different
positive numbers (including 0 and Infinity).
Using Automatic instead.";

ComplexPlot::badlrange =
"QLightnessRange should be Automatic or a list of two numbers
between 0 and 1. Using Automatic instead.";

ComplexPlot::badradius =
"QSphereRadius must be a positive real number.
Using QSphereRadius->1. instead.";

(*-----------------------------------*)
End[];      (* end `Private` Context *)
(*-----------------------------------*)

Protect @@ VQM`ComplexPlot`Private`Symbols;

(*-----------------------------------*)
EndPackage[]; (* end package Context *)
(*-----------------------------------*)

If[VQMmsgon, On[General::"spell"]];
If[VQMmsgon1, On[General::"spell1"]];


(* :Examples:

(* Turn off error messages *)
msgonindet = Head[General::"indet"]=!=$Off;
Off[Arg::"indet" ];

QComplexDensityPlot[x+I y, {x,-4,4}, {y,-4,4}]

gr1=QComplexDensityPlot[Sin[x + I y],
        {x,-Pi,Pi}, {y,-1,2.5},
        Mesh->False, PlotPoints->30]

(* There is a probelm here in 6 ... : gr1 is already Graphics[Raster[]]
Probably need another function ... 
Show[SurfaceGraphics[gr1],
    AspectRatio->Automatic,Axes->True,Mesh->True]
*)

QComplexContourPlot[x + I y, {x,-4,4}, {y,-4,4},
    QValueRange->{2,5},
    PlotRange->{0,3},Contours->5,
    ContourStyle->GrayLevel[0.5]]

QComplexContourPlot[Zeta[x + I y],
    {x,-.7,2.5}, {y,-2,42},
    PlotPoints->{10,50}, PlotRange->{-1.5,8},
    Contours->8, AspectRatio->5]

myfunc[{y_,m_,c_}] := Hue[-2 y,1-m,c];
collist=Table[myfunc[Mod[{x,x+y,x-y},1]],
        {y,0,1,1/20.}, {x,0,1,1/20.}];
QColorArrayPlot[collist,Mesh->False]

lis = Table[N[Tan[x + I y]],
        {y,-1.,1,.1},{x,-N[Pi],Pi,N[Pi]/10}];
QListComplexDensityPlot[lis,
    DataRange ->{{-Pi,Pi},{-1,1}}]

QComplexPlot3D[2 Exp[-x^2 - y^2 - 3 I x] Sin[x + I y],
    {x,-1,1},{y,-1,1}]

f[x_,y_]:=
   Which[Abs[x+I y] > 1.1, (*Indeterminate*) 0 ,
		Abs[x+I y] > 0.8, 0,
		Abs[x+I y] > 0.5, ComplexInfinity,
		True, DirectedInfinity[x+I y]];

Unprotect[Arg]; Arg[ComplexInfinity]=0;

QComplexDensityPlot[f[x,y],{x,-1,1},{y,-1,1},Compiled->False,
    QLightnessRange->{.2,.8}] 

QComplexDensityPlot[f[x,y],{x,-1,1},{y,-1,1},QValueChecking->Off,
    Compiled->False,QLightnessRange->{.2,.8}]

arr = Table[Exp[I 3 (x+y) - x^2-y^2] +
            Exp[-I 2 (x+y) - (x-.5)^2-(y-.5)^2],
               {y,-3,3,.1},{x,-3,3,.1}];
QListComplexSurfacePlot[arr,
    PlotRange->All, QSphereRadius->0.5,
    QLightnessRange->{0.2,1.}];               

QColorArrayPlot[
	Table[QComplexToRGBColor[x + I y], {y,-2.,2,.2}, {x,-2.,2,.2}],
        Mesh->False]

QSpinorDensityPlot[{Exp[I x],Exp[I y]}, {x,-Pi,Pi}, {y,-Pi,Pi}, 
	PlotPoints->20]

If[msgonindet, On[Arg::"indet"]]; (* turn on if it was on before *)

BUGPrint["C4a "];MP;
*)
