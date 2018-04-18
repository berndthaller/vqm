(* :Title: VisualizeFunction1D *)

(* :Name: VQM`VisualizeFunction1D` *)

(* :Copyright: Copyright 2005 Bernd Thaller *)

(* :Author:  Bernd Thaller,
             Institute of Mathematics,
             University of Graz,
             A-8010 Graz
             bernd.thaller@uni-graz.at
*)

(* :Source:
	Visual Quantum Mechanics
	Springer-Verlag New York, 2000
*)

(* :Summary:
VQM`VisualizeFunction1D` is a package for 'Visual Quantum Mechanics',
mainly supporting notebooks in Chapter 1.
It provides auxiliary commands for visualizing complex-valued functions
of one variable, without using the phase-coloring method of
the package ArgColorPlot.
*)

(* :Date:	2005-07-24 *)

(* :Package Version: 1.1
Changed for Mathematica 7 2009-06-10
*)

(* :Mathematica Version: 5.1 *)

(* :History:
    1.0 for Visual Quantum Mechanics, 2nd ed.
*)

VQMmsgon  = Head[General::"spell"]  =!= $Off;
VQMmsgon1 = Head[General::"spell1"] =!= $Off;

Off[General::spell1, General::spell];
    

(*-----------------------------------*)
BeginPackage["VQM`VisualizeFunction1D`"];
(*-----------------------------------*)

VQM`VisualizeFunction1D`Private`Symbols=Hold[
	QShowComplexPoint, QShowComplexPointPolar];

Unprotect@@VQM`VisualizeFunction1D`Private`Symbols;
ClearAll@@VQM`VisualizeFunction1D`Private`Symbols;
    
QShowComplexPoint::usage = "QComplexPoint[z] generates a visualization of a
	single complex number z in the complex plane.
	The number z is represented by a black dot with coordinates (Re(z), Im(z)).
	A red line symbolizes the real part of z, a yellow-green line symbolizes
	the imaginary part. The visible region of the complex plane may be
	specified by PlotRange->{{xmin,xmax},{ymin,ymax}}
	Package: VQM`VisualizeFunction1D`.";

QShowComplexPointPolar::usage = "QComplexPointPolar[z] generates a visualization
	of a single complex number z in the complex plane, together with some
	graphics elements symbolizing the polar coordinates of z.
	The number z is represented by a black dot with coordinates (Re(z), Im(z)).
	The modulus or absolute value |z| is shown by the length of a gray line.
	The phase or argument arg(z) is the angle between the gray line and
	the positive x-axis. Here arg(z) equals the length of a blue circular arc.
	It is assumed that arg(z) is in the interval 
	Package: VQM`VisualizeFunction1D`.";

QPlotRe::usage = "QPlotRe[f[x],{x,xmin,xmax},opts] plots only the real part of a
    complex-valued function as a red curve.
	Package: VQM`VisualizeFunction1D`.";

QPlotIm::usage = "QPlotIm[f[x],{x,xmin,xmax},opts] plots only the imaginary part
	of a complex-valued function. The function graph is yellow-green.
	Package: VQM`VisualizeFunction1D`.";

QPlotReIm::usage = "QPlotReIm[f[x],{x,xmin,xmax},opts] combines plots of the
	real part and of the imaginary part of the complex-valued function f into
	a single graph. The real part is shown in red, the imaginary part is
	yellow-green.
	Package: VQM`VisualizeFunction1D`.";

QPlotAbs::usage = "QPlotAbs[f[x],{x,xmin,xmax},opts] shows the absolute value
	of a complex-valued function f.
	Package: VQM`VisualizeFunction1D`.";

QPlotArg::usage = "QPlotArg[f[x],{x,xmin,xmax},opts] shows the argument (phase)
	of a comples-valued function f.
	Package: VQM`VisualizeFunction1D`.";

QPlotAbsArg::usage = "QPlotAbsArg[f[x],{x,xmin,xmax},opts] combines plots of
	absolute value (gray line) and argument (blue line) of a complex-valued
	function f.
	Package: VQM`VisualizeFunction1D`.";

QComplexFunctionGraph::usage = "QComplexFunctionGraph[f[x],{x,xmin,xmax},opts]
	shows the graph of a complex-valued function of a real variable x in 3D,
	with projections of the absolute value, the real part, and the imaginary
	part. The options QBoxSize and QProjectionAt set the dimensions of the
	bounding box and the positions of the projections.";

QBoxSize::usage = "QBoxSize->{{x1,x2},{y1,y2},{z1,z2}} is an option for 
	QComplexFunctionGraph, specifying the size of the boundary box.
	Its values are used for positioning coordinate axes, etc.
	Default value: Automatic."

QProjectionAt::usage = "QProjectionAt->{a,b,c} is an option for 
	QComplexFunctionGraph, specifying the location of the graphs
	of the real and imaginary parts, and the absolute value.
	Automatic projects the graphs onto the boundary box."



(*-----------------------------------*)
Begin["`Private`"];
(*-----------------------------------*)

(* Plot complex number z with Re and Im: *)

QShowComplexPoint[z_,opts___?OptionQ] :=
  Module[
    {zpoint={Re[z],Im[z]},
      realpartcolor = RGBColor[1,0,0], (* color of the line showing real part *)
      imagpartcolor=RGBColor[.5,1,0],  (* color of the line showing imag part *)
      dotcolor=GrayLevel[0],           (* color of the point z - actually a small disk *)
      relativeThickn = 0.015},		   (* determines thickness of line and diameter of point *)
    Show[
      Graphics[{
          {imagpartcolor,Thickness[relativeThickn],
            Line[{{zpoint[[1]],0},zpoint}]},{realpartcolor,
            Thickness[relativeThickn],
            Line[{{0,zpoint[[2]]},zpoint}]},{dotcolor,
            PointSize[2*relativeThickn],Point[zpoint]} 
          }], (*the diameter of the point is twice the diameter of the line*)
      opts,
      AspectRatio->Automatic,
      Axes->True
      ]
    ]

Options[QShowComplexPoint] = Options[Graphics];

QShowComplexPointPolar[z_,opts___?OptionQ]:=
  Module[
    {zpoint={Re[z],Im[z]},
      relativeThickn = 0.015,
      arglinecolor = RGBColor[0,0,1],
      abslinecolor = GrayLevel[0.4],
      dotcolor = GrayLevel[0],
      angles = If[Arg[z]>0,{0,Arg[z]},{Arg[z],0},{0,0}]
      },
    Show[
      Graphics[{
          {abslinecolor,Thickness[relativeThickn/3],
            Circle[{0,0},1]}, (*unit circle*)
          {abslinecolor,
            Thickness[relativeThickn],Line[{{0,0},zpoint}]},
          {abslinecolor,Thickness[relativeThickn/3],
            Line[{{0,0}, zpoint/Abs[z]}]},{RGBColor[0,0,1],Thickness[.015],
            Circle[{0,0},1,angles]},{dotcolor,PointSize[2*relativeThickn],
            Point[zpoint]}
          }],
      opts,
      AspectRatio->Automatic,
      Axes->True]
    ]

Options[QShowComplexPointPolar] = Options[Graphics];

QPlotRe[func_,{x_Symbol,a_,b_},opts___?OptionQ] :=
  Module[{f=Function[{x},func]},
    Plot[Re[f[x]],{x,a,b},
      Evaluate[Sequence @@ FilterRules[Flatten[{opts}], Options[Plot]]],
      PlotStyle->{Thickness[0.015],Hue[0]}]
    ]

Options[QPlotRe] = Options[Plot];

QPlotIm[func_,{x_Symbol,a_,b_},opts___?OptionQ] :=
  Module[{f=Function[{x},func]},
    Plot[Im[f[x]],{x,a,b},
      Evaluate[Sequence @@ FilterRules[Flatten[{opts}], Options[Plot]]],
      PlotStyle->{Thickness[0.015],Hue[0.25]}]
    ]

Options[QPlotIm] = Options[Plot];

QPlotReIm[func_,{x_Symbol,a_,b_},opts___?OptionQ] :=
  Module[{f=Function[{x},func],gr1,gr2},
  	gr1 = QPlotRe[f[x], {x, a, b}, opts];
	gr2 = QPlotIm[f[x], {x, a, b}, opts];
    Show[{gr1, gr2},
 		Evaluate[Sequence @@ FilterRules[Flatten[{opts}], Options[Graphics]]]
 		]
    ]

Options[QPlotReIm] = Options[Plot];

QPlotAbs[func_,{x_Symbol,a_,b_},opts___?OptionQ] :=
  Module[{f=Function[{x},func]},
    Plot[Abs[f[x]],{x,a,b},
      Evaluate[Sequence @@ FilterRules[Flatten[{opts}], Options[Plot]]],
      PlotStyle->{Thickness[0.015],GrayLevel[0.4]}]
    ]

Options[QPlotAbs] = Options[Plot];

QPlotArg[func_,{x_Symbol,a_,b_},opts___?OptionQ] :=
  Module[{f=Function[{x},func]},
    Plot[Arg[f[x]],{x,a,b},
      Evaluate[Sequence @@ FilterRules[Flatten[{opts}], Options[Plot]]],
      PlotStyle->{Thickness[0.015],Hue[0.66]}]
    ]

Options[QPlotArg] = Options[Plot];

QPlotAbsArg[func_,{x_Symbol,a_,b_},opts___?OptionQ] :=
  Module[{f=Function[{x},func],gr1,gr2},
  	gr1=QPlotAbs[f[x],{x,a,b},opts];
    gr2=QPlotArg[f[x],{x,a,b},opts];
    Show[{gr1,gr2},
      Evaluate[Sequence @@ FilterRules[Flatten[{opts}], Options[Graphics]]]
      ]
    ]

Options[QPlotAbsArg] = Union[Options[Plot],Options[Graphics]];

axes[x1_, x2_, y1_, y2_, z1_, z2_] :=
  
  Graphics3D[{Dashing[{0.02, 0.02}], Thickness[0.0075], GrayLevel[0],
    Line[{{x1, 0, 0}, {x2, 0, 0}}],
    Line[{{0, y1, 0}, {0, y2, 0}}],
    Line[{{0, 0, z1}, {0, 0, z2}}]}];


prAxes[x1_, x2_, y1_, y2_, z1_, z2_, a_, b_, c_] :=
  Module[{xp},
   Graphics3D[{Dashing[{0.02, 0.02}], Thickness[0.0075], 
     GrayLevel[0.75],
     (*x-axis projected on surfaces y=b and z=c*)
     
     Line[{{x1, b , 0}, {x2, b , 0}}],
     Line[{{x1, 0, c }, {x2, 0, c }}],
     (*y-axis projected on surfaces x=a and z=c*)
     
     Line[{{0, y1, c }, {0, y2, c }}],
     Line[{{a, y1, 0}, {a, y2, 0}}],
     (*z-axis projected on surfaces y=b and x=a*)
     
     Line[{{0, b , z1}, {0, b , z2}}],
     Line[{{a, 0, z1}, {a, 0, z2}}]
     }]
   ];


myoptions =
	{QBoxSize -> Automatic,
	QProjectionAt->Automatic};

QComplexFunctionGraph[func_, {x_, le_, ri_}, opts___?OptionQ] := 
	Module[{f = Function[{x}, func], curve, imcurve, recurve, prcurve, 
			thickn = 0.015, box, pr, x1, x2, y1, y2, z1, z2
			}, 
		{box, pr} = {QBoxSize, QProjectionAt} /. {opts} /. myoptions;
		If[(box === Automatic) || (Dimensions[box] != {3, 2}),
			box = {{le, ri}, {-2, 2}, {-2, 2}}];
		{{x1, x2}, {y1, y2}, {z1, z2}} = box;
		If[(pr === Automatic) || (Dimensions[pr] != 3),
			pr = {x1,y2,z1}];
		{a,b,c} = pr;
  curve = 
   ParametricPlot3D[Evaluate[{t, Re[f[t]], Im[f[t]]}], {t, le, ri}, 
    PlotPoints -> 400, PlotStyle -> {Thickness[thickn], GrayLevel[0]}];
  imcurve = 
   ParametricPlot3D[Evaluate[{t, b, Im[f[t]]}], {t, le, ri}, 
    PlotPoints -> 400, 
    PlotStyle -> {Thickness[thickn], Hue[0.25, 0.5, 1]}];
  recurve = 
   ParametricPlot3D[Evaluate[{t, Re[f[t]], c}], {t, le, ri}, 
    PlotPoints -> 400, 
    PlotStyle -> {Thickness[thickn], Hue[0, 0.5, 1]}];
  prcurve = 
   ParametricPlot3D[Evaluate[{a, Re[f[t]], Im[f[t]]}], {t, le, ri}, 
    PlotPoints -> 400, 
    PlotStyle -> {Thickness[thickn], GrayLevel[0.75]}];
  Show[{curve, imcurve, recurve, prcurve, 
    axes[x1, x2, y1, y2, z1, z2], 
    prAxes[x1, x2, y1, y2, z1, z2, a, b, c]},
    Evaluate[Sequence @@ FilterRules[Flatten[{opts}], Options[Graphics3D]]], 
   BoxStyle -> {GrayLevel[0], Thickness[0.0025]}, PlotRange -> All, 
   Axes -> True, AxesStyle -> Thickness[0.0075], 
   AxesLabel -> {"x", "Re", "Im"}, LabelStyle -> 14, 
   AxesEdge -> {{-1, -1}, {1, -1}, {1, 1}}, 
   Ticks -> {Automatic, {-1, 0, 1}, {-1, 0, 1}}]]


(*-----------------------------------*)
End[];      (* end `Private` Context *)
(*-----------------------------------*)

Protect@@VQM`VisualizeFunction1D`Private`Symbols;

(*-----------------------------------*)
EndPackage[]; (* end package Context *)
(*-----------------------------------*)

If[VQMmsgon, On[General::"spell"]];
If[VQMmsgon1, On[General::"spell1"]];

(* Examples:
	QShowComplexPoint[3+2 I, PlotRange->{{-2,4},{-1.5,3}}]
	QShowComplexPointPolar[3+2 I, PlotRange->{{-2,4},{-1.5,3}}]
	QPlotReIm[Exp[2 I x],{x,-Pi,Pi}]
	QPlotAbsArg[Exp[2 I x], {x,-2 Pi,2 Pi}, Axes->True, PlotRange->{-5, 5}]
	QComplexFunctionGraph[Exp[I*3*x - x^2/2], {x, -Pi, Pi}]
	QComplexFunctionGraph[Exp[I*3*x - x^2/2], {x, -Pi, Pi},
         QBoxSize -> {{-Pi, Pi}, {-1.5, 1.5}, {-1.5, 1.5}},  QProjectionAt -> {0, -1.5, 1.5}]
*)
