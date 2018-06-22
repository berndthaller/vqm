(* :Title: VisualizeFourierSum1D *)

(* :Name: VQM`VisualizeFourierSum1D` *)

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
VQM`VisualizeFourierSum1D` is a package for 'Visual Quantum Mechanics',
mainly supporting notebooks in Chapter 1.
It provides auxiliary commands for visualizing complex-valued functions
of one variable, without using the phase-coloring method of
the package ArgColorPlot.
*)

(* :Date:	2007-01-02 *)
(* :Date:	2018-06-22 *)

(* :Package Version: 1.0 *)

(* :Mathematica Version: 11.3 *)

(* :History:
    1.0 for Visual Quantum Mechanics, 2nd ed.
*)

VQMmsgon  = Head[General::"spell"]  =!= $Off;
VQMmsgon1 = Head[General::"spell1"] =!= $Off;

Off[General::spell1, General::spell];
    

(*-----------------------------------*)
BeginPackage["VQM`VisualizeFourierSum1D`",
	{"VQM`ArgColorPlot`"}
];
(*-----------------------------------*)

VQM`VisualizeFourierSum1D`Private`Symbols=Hold[
	QWaveNumberPlot, QFourierSumPlot, QWaveStackPlot, QFourierSumWaveNumberPlot,
	QWaveStackWaveNumberPlot, QFourierSumWaveStackPlot, QStackWavesGraph,
	QVisualizeSuperposition, QXRegion, QKRegion, QLineSize, QEpilog, QBorders,
	QXPlotLabel, QKPlotLabel, QSPlotLabel, QAmplitudeFactor, QHalfPeriod,
	QSeparationLines];

Unprotect@@VQM`VisualizeFourierSum1D`Private`Symbols;
ClearAll@@VQM`VisualizeFourierSum1D`Private`Symbols;
    
QWaveNumberPlot::usage = "QWaveNumberPlot[{{k1,A1},...,{kn,An}}, options]
	plots vertical lines of length N[L]*Abs[Aj] at positions kj. Here N[L] is a
    normalization constant set by the option QAmplitudeFactor. The color of the
    lines is determined by the complex phases of the coefficients Aj.
    The thickness of the lines is changed by QLineSize->d, where d
    is the thickness of the lines relative to the coordinates of the plot
    (default is 0.2). The plot region is determined automatically, but can also
    be specified by the option QKRegion->{{kleft,kright},{klower,kupper}}.
    You may add additional graphics elements by using the option QEpilog. 
	Package: VQM`VisualizeFourierSum1D`.";

QFourierSumPlot::usage = "QFourierSumPlot[{{k1,A1},...,{kn,An}}, options]
	plots a superposition (linearcombination) of the plane waves
	Aj Pi/(L Sqrt[2 Pi]) Exp[I kj x]. Here it is assumed that all plane waves,
	and hence the superposition, have the periodicity interval [-L,L].
	We attempt to determine L from the kj. It can also be specified by the
	option QHalfPeriod->L. The plot region can be given as
	QXRegion->{{xleft,xright},{xlower,xupper}}. The periodicity interval
	is indicated by vertical gray lines in the plot. This can be supressed
	by QBorders->False. Usually, the normalization constant of the plane wave
	can be changed by setting the option QAmplitudeFactor. The setting
	QAmplitudeFactor->1 simply plots the sum of Aj*Exp[I kj x].
	Default is QHalfPeriod->L, QAmplitudeFactor->Pi/(L Sqrt[2 Pi]). 
	Package: VQM`VisualizeFourierSum1D`.";

QWaveStackPlot::usage = "QWaveStackPlot[{{k1,A1},...,{kn,An}}, options] 
	plots plane waves with wave numbers kj and (complex) amplitudes N[L] An,
	one on top of the other. The factor N[L] is by default Pi/(L Sqrt[2 Pi])
	but it can also be specified with the option QAmplitudeFactor.
	Here[-L,L] is the interval of periodicity. L is determined from the
	smallest nonzero Abs[kj], but it can also be set with the option
	QHalfPeriod->L. We try to guess the plot region, but it can also be
	specified by QXRegion->{{xleft,xright},{xlower,xupper}}.
	With the default QBorders->True, the interval of periodicity is outlined
	by vertical gray lines, QBorders->False supresses those lines.
	The plane waves are separated by thin black lines, QSeparationLines->False
	supresses these lines. 
	Package: VQM`VisualizeFourierSum1D`.";

QFourierSumWaveNumberPlot::usage = "QFourierSumWaveNumberPlot[{{k1, A1}, ..., {kn, An}}, options]
	combines QWaveNumberPlot and QFourierSumPlot into a GraphicsArray.
	The graphics can be labeled separately; use the option QXPlotLabel
	for the Fourier sum and QKPlotLabel for the plot of wave number lines.
	Simply use the option PlotLabel for the combined graphics.
	Note that if you label only one of the graphs in the array,the alignment
	might be poor. Give the options QKRegion, QXRegion to manipulate the
	plot ranges of the graphs separately. QEpilog applies only to the
	wave number plot, so you can use this option without disturbing
	the plot of the Fourier sum. Use Epilog to plot something over the
	Fourier sum without disturbing the graph showing the wave number lines. 
	Package: VQM`VisualizeFourierSum1D`.";

QWaveStackWaveNumberPlot::usage = "QWaveStackWaveNumberPlot[{{k1,A1},...,{kn,An}}, options]
	combines QWaveNumberPlot and QWaveStackPlot into a GraphicsArray.
	The graphics can be labeled separately; use the option QSPlotLabel
	for the plot of stacked plane waves and QKPlotLabel for wave number plot.
	Simply use PlotLabel for the combined graph. Note that if you label
	only one of the graphs in the array, the alignment would be poor.
	Give the options QKRegion, QXRegion to manipulate the
	plot ranges of the graphs separately. QEpilog applies only to the
	wave number plot, so you can use this option without disturbing
	the plot of the plane waves.  
	Package: VQM`VisualizeFourierSum1D`.";

QFourierSumWaveStackPlot::usage = "QFourierSumWaveStackPlot[{{k1,A1},...,{kn,An}}, options]
	combines QWaveStackPlot and QFourierSumPlot into a GraphicsArray.
	The graphics can be labeled separately; use the option QXPlotLabel
	for the Fourier sum and QSPlotLabel for the plot of stacked plane waves.
	Simply use PlotLabel for the combined graph. Note that if you label
	only one of the graphs in the array, the alignment would be poor.
	QXRegion->{{xleft,xright},{xlower,xupper}} gives the plot range
	for both parts of the graph. 
	Package: VQM`VisualizeFourierSum1D`.";

QStackWavesGraph::usage = "QStackWavesGraph[{{k1,A1},...,{kn,An}},  options]
	is similar to QFourierSumWaveStackPlot, but slightly different in
	appearance, because it combines QWaveStackPlot and the QFourierSumPlot within
	a single frame. 
	Package: VQM`VisualizeFourierSum1D`.";

QVisualizeSuperposition::usage = "QVisualizeSuperposition[{{k1,A1},...,{kn,An}},  options]
	combines QWaveNumberPlot,QWaveStackPlot,and QFourierSumPlot into a
	single GraphicsArray with neatly arranged parts. The parts can be labeled
	separately, use the option QXPlotLabel for the Fourier sum, QSPlotLabel
	for the plot of stacked plane waves, and QKPlotLabel for the wave number plot.
	Simply use the option PlotLabel for the combined graphics.
	Note that if you label only one part of the array, the alignment would be poor.
	Describe the horizontal and vertical plot ranges with the options QXRegion and
	QKRegion. 
	Package: VQM`VisualizeFourierSum1D`.";

QXRegion::usage = "QXRegion->{{xleft,xright},{xlower,xupper}} allows to specify
	the plot range for QWaveStackPlot, QFourierSumPlot, and related plot commands.
	Package: VQM`VisualizeFourierSum1D`.";
	
QKRegion::usage = "QKRegion->{{kleft,kright},{klower,kupper}} allows to specify
	the plot range for QWaveNumberPlot, and related plot commands.
	Package: VQM`VisualizeFourierSum1D`.";
	
QLineSize::usage = "QLineSize->d is an option for QWaveNumberPlot. It allows
	to specify the thickness of the wave number lines.
	Package: VQM`VisualizeFourierSum1D`.";
	
QEpilog::usage = "QEpilog->graphics is an option for QWaveNumberPlot. It places
	additional graphics elements like the option Epilog in ordinary Plot
	commands. In combined commands like QFourierSumWaveNumberPlot one can use
	QEpilog to place graphics elements in the wave number plot, and Epilog
	to place graphics elements in the Fourier sum plot. 
	Package: VQM`VisualizeFourierSum1D`.";
	
QBorders::usage = "QBorders->True is an option for QWaveStackPlot and QFourierSumPlot.
	It places vertical lines at x=+L and x=-L in order to mark the borders of the
	interval of periodicity. L can be changed by the option QHalfPeriod. 
	Package: VQM`VisualizeFourierSum1D`.";

QXPlotLabel::usage = "QXPlotLabel->text is an option for QVisualizeSuperposition,
	QFourierSumWaveStackPlot, and QFourierSumWaveNumberPlot. It places a plot label
	on the part of the plot that displays the Fourier sum in x-space. 
	Package: VQM`VisualizeFourierSum1D`.";

QKPlotLabel::usage = "QKPlotLabel->text is an option for QVisualizeSuperposition,
	and QFourierSumWaveNumberPlot. It places a plot label
	on the part of the plot that displays the wave number lines in k-space. 
	Package: VQM`VisualizeFourierSum1D`.";

QSPlotLabel::usage = "QSPlotLabel->text is an option for QVisualizeSuperposition,
	and QFourierSumWaveStackPlot. It places a plot label
	on the part of the plot that displays the stack of plane waves in x-space. 
	Package: VQM`VisualizeFourierSum1D`.";

QAmplitudeFactor::usage = "QAmplitudeFactor->N sets the norming factor for the
	plane wave in QWaveStackPlot, QFourierSumPlot, and related plot commands.
	By default, the norming factor depends on the interval of periodicity [-L,L]
	and is given by Pi/(L Sqrt[2 Pi]). With QAmplitudeFactor->1 in
	QFourierSumWaveNumberPlot, and QVisualizeSuperposition, etc, the length
	of the wave number lines and the amplitude of the corresponding plane
	waves are equal. 
	Package: VQM`VisualizeFourierSum1D`.";
	
QHalfPeriod::usage = "QHalfPeriod->L is an option for QFourierSumPlot, QWaveStackPlot,
	and related plot commands. It specifies the interval of periodicity [-L,L] for
	all functions. L is used to determine the plot region for QFourierSumPlot and
	QWaveStackPlot, and the interval [-L,L] is indicated by vertical lines unless
	the option QBorders is set to False. QFourierSumPlot and QWaveStackPlot
	try to determine L automatically from the provided wave numbers as L = Pi/Abs[k1],
	where k1 is the nonzero wave number closest to k=0, but this might fail
	(for example, if the wave number k1 is missing in the superposition).
	It is not checked whether the choice made for L is consistent
	with the wave numbers provided to the plot commands. 
	Package: VQM`VisualizeFourierSum1D`.";
	
QSeparationLines::usage = "QSeparationLines->True is an option for QWaveStackPlot.
	This command plots several plane waves on top of each other, separated by
	horizontal lines. Setting this option ot False supresses the separation
	lines, which looks better if there is a large number of plane waves. 
	Package: VQM`VisualizeFourierSum1D`.";

(*-----------------------------------*)
Begin["`Private`"];
(*-----------------------------------*)

(* -- Auxiliary -- *)

RemoveOptions[list1_, list2___] :=
  Sequence @@ Sort[Select[Flatten[{list2}],  !MemberQ[list1, First[#1]] & ]];

JoinOptions[list1_, list2___] := 
  Module[{namelist = First /@ Flatten[{list1}]}, 
   Sequence @@ Sort[Join[Flatten[{list1}], Select[Flatten[{list2}], 
        !MemberQ[namelist, First[#1]] & ]]]];

(* -- QWaveNumberPlot -- *)

colorLine[{mom_, coeff_}, linethickness_] :=
  {Hue[Mod[Arg[coeff]/2/Pi, 1]],
    Polygon[{{mom - linethickness/2, Abs[coeff]}, {mom - linethickness/2, 0.},
      {mom + linethickness/2, 0.}, {mom + linethickness/2, Abs[coeff]}}]}; 

lineTable[listOfPlaneWaves_, linethickness_] := (colorLine[#1, linethickness] & ) /@ listOfPlaneWaves;

Options[QWaveNumberPlot] = RemoveOptions[{PlotRange, PlotStyle, QShiftPlot, QBottomLine, Epilog}, 
   Join[Options[QArgColorPlot], {QLineSize -> 0.2, QKRegion -> Automatic, QEpilog -> {}}]];

QWaveNumberPlot[listOfPlaneWaves_, opts___] := 
   Module[{linethickness = QLineSize /. {opts} /. Options[QWaveNumberPlot], 
         plotreg = QKRegion /. {opts} /. Options[QWaveNumberPlot], 
         kmin = Min[First /@ listOfPlaneWaves],
         kmax = Max[First /@ listOfPlaneWaves], 
         Amax = Max[Abs /@ Last /@ listOfPlaneWaves], 
         myopts = FilterRules[Flatten @ {opts}, Options @ QWaveNumberPlot],
         kleft, kright, klower, kupper, msgon, epiloggr, gr}, 
      If[plotreg == Automatic,
      plotreg = {{kmin - 3*linethickness, kmax + 3*linethickness}, {-0.1*Amax, 1.1*Amax}}]; 
      {{kleft, kright}, {klower, kupper}} = plotreg; 
      msgon = Head[Arg::indet] =!= $Off; Off[Arg::indet]; 
      epiloggr = Join[lineTable[listOfPlaneWaves, linethickness], QEpilog /. {opts} /. Options[QWaveNumberPlot]]; 
      gr = QArgColorPlot[0, {k, kleft, kright},
         Compiled -> False, PlotRange -> {klower, kupper}, myopts, Epilog -> epiloggr, 
         AxesLabel -> {"k", None}];
      If[msgon, On[Arg::indet]];
      gr];
      
(* Note: QWaveNumberPlot uses QArgColorPlot only because want to combine the result
in a graphics array with other QArgColorPlots. The plot of the function 0 just provides
a plot region that aligns well with other QArgColorPlots. *)

(* -- QFourierSumPlot -- *)

Options[QFourierSumPlot] =
  Join[Options[QArgColorPlot],
    {QXRegion -> Automatic, QHalfPeriod -> Automatic, QAmplitudeFactor -> Automatic, QBorders -> True}];

QFourierSumPlot[listOfPlaneWaves_, opts___] := 
   Module[{superposition, xleft, xright, xlower, xupper, maxval, epiloggr={},
         L = QHalfPeriod /. {opts} /. Options[QFourierSumPlot],
         fac = QAmplitudeFactor /. {opts} /. Options[QFourierSumPlot],
         plotreg = QXRegion /. {opts} /. Options[QFourierSumPlot], 
         borders = QBorders /. {opts} /. Options[QFourierSumPlot]}, 
      If[L === Automatic, L = Pi/Min[Select[Abs /@ First /@ listOfPlaneWaves, (# != 0)& ]]];
      If[L==0, L=Pi];
      If[fac === Automatic, fac = Sqrt[Pi/2]/L];
      If[plotreg === Automatic, plotreg = {{-1.05 L, 1.05 L}, Automatic}];
      {xleft,xright}=plotreg[[1]];
      maxval = fac*(Plus @@ Abs /@ Last /@ listOfPlaneWaves);
      If[plotreg[[2]]===Automatic, {xlower,xupper}={-0.1 maxval, 1.1 maxval}, {xlower,xupper}=plotreg[[2]]];
      If[borders===True, epiloggr = {GrayLevel[0.5], Thickness[0.01], 
          Line[{{-L, xlower}, {-L, xupper}}], Line[{{L, xlower}, {L, xupper}}]}, epiloggr={}];
      superposition[x_] = fac*Exp[I*First /@ listOfPlaneWaves*x] . Last /@ listOfPlaneWaves;
      QArgColorPlot[superposition[x], {x, xleft, xright}, opts, PlotRange -> Evaluate[plotreg[[2]]], 
        AxesLabel -> {"x", None}, Epilog -> epiloggr]]; 
        
(* -- QWaveStackPlot -- *)

Options[QWaveStackPlot] = Join[Options[QArgColorPlot], 
    {QXRegion -> Automatic, QHalfPeriod -> Automatic, 
     QAmplitudeFactor -> Automatic, QBorders -> True, 
     QSeparationLines -> True}];
     
QWaveStackPlot[listOfPlaneWaves_, opts___] := 
  Module[{plwav, height, maxval, gr,
        L = QHalfPeriod /. {opts} /. Options[QWaveStackPlot],
        fac = QAmplitudeFactor /. {opts} /. Options[QWaveStackPlot],
        plotreg = QXRegion /. {opts} /. Options[QWaveStackPlot],
        borders = QBorders /. {opts} /. Options[QWaveStackPlot],
        sep = QSeparationLines /. {opts} /. Options[QWaveStackPlot]}, 
    If[L === Automatic,
       L = Pi/Min[Select[Abs /@ First /@ listOfPlaneWaves, #1 != 0 & ]]];
    If[L == 0, L = Pi];
    If[fac === Automatic, fac = Sqrt[Pi/2]/L];
     plwav[m_, x_] := 
     fac*Exp[I*listOfPlaneWaves[[m,1]]*x]*listOfPlaneWaves[[m,2]]; 
    height[m_] := fac*Plus @@ Abs[N[Last /@ Take[listOfPlaneWaves, m]]]; 
    maxval = height[Length[listOfPlaneWaves]];
    If[plotreg === Automatic, 
     plotreg = {{-1.05*L, 1.05*L}, {-0.1*maxval, 1.1*maxval}}, Null]; 
    {{xleft, xright}, {xlower, xupper}} = plotreg; 
    If[borders === True,
      epiloggr = {GrayLevel[0.5], Thickness[0.01], 
        Line[{{-L, xlower}, {-L, xupper}}],
        Line[{{L, xlower}, {L, xupper}}]}, 
      epiloggr = {}];
    msgon = Head[Arg::indet] =!= $Off;
    Off[Arg::indet]; 
    If[sep === True, 
      gr =
        Show[
        {Table[QArgColorPlot[plwav[m, x], {x, xleft, xright}, 
           QShiftPlot -> Evaluate[height[m - 1]],
           Evaluate[FilterRules[Flatten[{opts}], Options @ QArgColorPlot]], 
           AxesLabel -> {"x", None}, PlotRange -> {xlower, xupper}], 
        {m, 1, Length[listOfPlaneWaves]}
        ],
        Graphics[epiloggr]}, 
        Evaluate[FilterRules[Flatten[{opts}], Options @ QArgColorPlot]]
        ]
      , 
      gr =
        Show[
        Flatten[{Table[Graphics[{EdgeForm[],
          QArgColorPlot[plwav[m, x], {x, xleft, xright},
            QShiftPlot -> Evaluate[height[m - 1]]
            ][[1,1]]}], 
           {m, 1, Length[listOfPlaneWaves]}],
           Graphics[epiloggr]}], 
        Evaluate[FilterRules[Flatten[{opts}], Options @ Graphics]],
        Axes -> True, 
        AxesFront -> True,
        AxesLabel -> {"x", None}, 
        PlotRange -> {xlower, xupper}
        ];
    ];
    If[msgon, On[Arg::indet]];
    gr
  ]

(* -- QWaveStackWaveNumberPlot -- *)

Options[QWaveStackWaveNumberPlot] =
	RemoveOptions[{GridLines, Epilog}, Join[Options[Graphics], {QSPlotLabel -> None, QKPlotLabel -> None}]];
	
QWaveStackWaveNumberPlot[listOfPlaneWaves_, opts___] :=
   Module[{myopts = RemoveOptions[{PlotLabel}, opts], 
         combinedopts = FilterRules[Flatten[{opts}], Join[Options @ Graphics, Options @ QWaveStackWaveNumberPlot]], 
         slabel = QSPlotLabel /. {opts} /. Options[QWaveStackWaveNumberPlot], 
         klabel = QKPlotLabel /. {opts} /. Options[QWaveStackWaveNumberPlot]}, 
      Show[
         Graphics[
            {Rectangle[{0, 0}, {1, 1/2},
                QWaveNumberPlot[listOfPlaneWaves, myopts,
                   PlotLabel -> klabel, Frame -> True, 
                   Axes -> {True, None}, AspectRatio -> 1/2.2
                ]
             ],
             Rectangle[{0, 1/2}, {1, 1}, 
                QWaveStackPlot[listOfPlaneWaves, myopts,
                   PlotLabel -> slabel, Frame -> True,
                   Axes -> {True, None}, AspectRatio -> 1/2.2
                ]
             ]}
         ], combinedopts, AspectRatio -> 1
      ]
   ]

(* -- QFourierSumWaveNumberPlot -- *)

Options[QFourierSumWaveNumberPlot] =
    RemoveOptions[{GridLines, Epilog}, Join[Options[Graphics], {QKPlotLabel -> None, QXPlotLabel -> None}]];

QFourierSumWaveNumberPlot[listOfPlaneWaves_, opts___] :=
   Module[{myopts = RemoveOptions[{PlotLabel}, opts], 
         combinedopts = FilterRules[Flatten[{opts}], Join[Options[Graphics], Options[QFourierSumWaveNumberPlot]] ],
         xlabel = QXPlotLabel /. {opts} /. Options[QFourierSumWaveNumberPlot], 
         klabel = QKPlotLabel /. {opts} /. Options[QFourierSumWaveNumberPlot]}, 
      Show[
         Graphics[
            {Rectangle[{0, 0}, {1, 1/2},
                QWaveNumberPlot[listOfPlaneWaves, myopts,
                   PlotLabel -> klabel, Frame -> True, 
                   Axes -> {True, None}, AspectRatio -> 1/2.2
                ]
             ],
             Rectangle[{0, 1/2}, {1, 1}, 
                QFourierSumPlot[listOfPlaneWaves, myopts,
                   PlotLabel -> xlabel, Frame -> True,
                   Axes -> {True, None}, AspectRatio -> 1/2.2
                ]
             ]}
         ], combinedopts, AspectRatio -> 1
      ]
   ]


(* -- QFourierSumWaveStackPlot -- *)

Options[QFourierSumWaveStackPlot] = RemoveOptions[{GridLines, Epilog}, 
    Join[Options[Graphics], {QXPlotLabel -> None, QSPlotLabel -> None}]];

QFourierSumWaveStackPlot[listOfPlaneWaves_, opts___] :=
   Module[{myopts = RemoveOptions[{PlotLabel}, opts], 
         combinedopts = FilterRules[ Flatten[{opts}], Join[Options[Graphics], Options[QFourierSumWaveStackPlot]] ], 
         xlabel = QXPlotLabel /. {opts} /. Options[QFourierSumWaveStackPlot], 
         slabel = QSPlotLabel /. {opts} /. Options[QFourierSumWaveStackPlot]}, 
      Show[
         Graphics[
            {Rectangle[{0, 0}, {1, 1/2},
                QWaveStackPlot[listOfPlaneWaves, myopts,
                   PlotLabel -> slabel, Frame -> True, 
                   Axes -> {True, None}, AspectRatio -> 1/2.2
                ]
             ],
             Rectangle[{0, 1/2}, {1, 1}, 
                QFourierSumPlot[listOfPlaneWaves, myopts,
                   PlotLabel -> xlabel, Frame -> True,
                   Axes -> {True, None}, AspectRatio -> 1/2.2
                ]
             ]}
         ], combinedopts, AspectRatio -> 1
      ]
   ]


(* -- QStackWavesGraph -- *)

Options[QStackWavesGraph] = {JoinOptions[Options[QFourierSumPlot], 
     Options[QWaveStackPlot]]};

QStackWavesGraph[listOfPlaneWaves_, opts___] := 
  Module[{ posY, yticks, tickstep, 
    L = QHalfPeriod /. {opts} /. Options[QStackWavesGraph], 
    fac = QAmplitudeFactor /. {opts} /. Options[QStackWavesGraph], 
    plotreg = QXRegion /. {opts} /. Options[QStackWavesGraph], maxval, whitespaceY
    }, 
   If[L === Automatic,
      L = Pi/Min[Select[Abs /@ First /@ listOfPlaneWaves, #1 != 0 & ]]];
   If[L == 0, L = Pi];
   If[fac === Automatic, fac = Sqrt[Pi/2]/L];
   maxval = fac*Plus @@ Abs /@ N /@ Last /@ listOfPlaneWaves; 
   whitespaceY = maxval/10;
   posY = whitespaceY + maxval; 
   tickstep = Which[maxval > 3,  Round[maxval/3.],
                    maxval > 0.6,  0.5, 
                    maxval > 0.25, 0.2,
                    maxval > 0.1, 0.1];
   If[plotreg === Automatic, 
      plotreg = {{-1.05*L, 1.05*L},
         {-whitespaceY, posY + maxval + whitespaceY/2}}];
   {{xleft, xright}, {xlower, xupper}} = plotreg;
   yticks = Table[{posY + k, k}, {k, 0, Ceiling[maxval], tickstep}];
   Show[
      {
       QFourierSumPlot[listOfPlaneWaves,
          QHalfPeriod -> L, 
          QXRegion -> {{xleft, xright}, {xlower, xupper}}, 
          Evaluate[FilterRules[Flatten[{opts}], Options[QFourierSumPlot]]],
          QShiftPlot -> posY
       ],
       QWaveStackPlot[listOfPlaneWaves, 
          QHalfPeriod -> L,
          QXRegion -> {{xleft, xright}, {xlower, xupper}}, 
          Evaluate[FilterRules[Flatten[{opts}], Options@QWaveStackPlot]] 
       ]
      },
      Evaluate[Flatten[{opts}], Options[Graphics]],
      Frame -> True, FrameTicks -> {Automatic, yticks, None, None}, 
      Epilog -> Line[{{xleft, posY}, {xright, posY}}],
      AspectRatio -> 1
      (*
      TODO: find out if something like AxesFront is still needed
      , 
      AxesFront -> True
      *)
   ]
]
  
(* -- QVisualizeSuperposition -- *)

Options[QVisualizeSuperposition] =
    RemoveOptions[{GridLines, Epilog}, 
      Join[Options[Graphics],{QXPlotLabel -> None, QKPlotLabel -> None,
          QSPlotLabel -> None, QXRegion -> Automatic,
          QKRegion -> Automatic, QHalfPeriod -> Automatic,
          QAmplitudeFactor -> Automatic, QBorders -> True, 
          QSeparationLines -> True, QLineSize -> 0.2, QEpilog -> {}}]];

QVisualizeSuperposition[listOfPlaneWaves_, opts___] := 
  Module[{myopts = RemoveOptions[{PlotLabel}, opts], 
    combinedopts = FilterRules[Flatten[{opts}], Join[Options @ Graphics, Options @ QVisualizeSuperposition]],
    xlabel = QXPlotLabel /. {opts} /. Options[QVisualizeSuperposition], 
    slabel = QSPlotLabel /. {opts} /. Options[QVisualizeSuperposition], 
    klabel = QKPlotLabel /. {opts} /. Options[QVisualizeSuperposition]}, 
   Show[
     Graphics[
       {Rectangle[{0, 0}, {1, 1/3},
         QWaveNumberPlot[ listOfPlaneWaves,
           myopts, PlotLabel -> klabel, Frame -> True, 
           Axes -> {True, None}, AspectRatio -> 1/3.2
           ]
        ],
        Rectangle[{0, 1/3}, {1, 2/3},
         QWaveStackPlot[ listOfPlaneWaves,
           myopts, PlotLabel -> slabel, Frame -> True, 
           Axes -> {True, None}, AspectRatio -> 1/3.2
           ]
        ],
        Rectangle[{0, 2/3}, {1, 1},
         QFourierSumPlot[ listOfPlaneWaves,
           myopts, PlotLabel -> xlabel, Frame -> True, 
           Axes -> {True, None}, AspectRatio -> 1/3.2
         ]
        ]}
     ],
     combinedopts, AspectRatio -> 1
   ]
  ]




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

k[n_, L_] := n Pi/L;

exampleListOfPlaneWaves = {{-1, 4}, {0, 5 + 5*I}, {1, 2 - 2*I}, {2, -1}, 
    {3, 2*I}, {5, -5*I}, {4, 0}}; 

sinListOfPlaneWaves = Array[{#1, Sin[#1]} & , 5]; 

GaussListOfPlaneWaves[L_] := Table[{k[n, L], Exp[-k[n, L]^2/2]}, 
   {n, -Ceiling[L], Ceiling[L]}]

GaussListOfPlaneWavesShift[L_] := 
  Table[{k[n, L], Exp[-k[n, L]^2/2 + I*k[n, L]]}, 
   {n, -Ceiling[L], Ceiling[L]}]

gr = Plot[Exp[-k^2/2], {k, -Pi, Pi}, PlotStyle -> GrayLevel[0.5], 
     DisplayFunction -> Identity][[1]];

QWaveNumberPlot[exampleListOfPlaneWaves]
QFourierSumPlot[exampleListOfPlaneWaves]

QWaveNumberPlot[GaussListOfPlaneWavesShift[3*Pi], 
  QKRegion -> {{-3.2, 3.2}, {-0.3, 1.2}}, QLineSize -> 0.4, Frame -> True, 
  AxesLabel -> None, PlotStyle -> Thickness[0.15], PlotRange -> {-5, 5}]

QWaveNumberPlot[GaussListOfPlaneWaves[3*Pi], 
  QKRegion -> {{-3.2, 3.2}, {-0.3, 1.2}}, QLineSize -> 0.4, Frame -> True, 
  AxesLabel -> None, PlotStyle -> Thickness[0.15], PlotRange -> {-5, 5}]

QWaveNumberPlot[GaussListOfPlaneWaves[3*Pi], 
  QKRegion -> {{-3.2, 3.2}, {-0.3, 1.2}}, QLineSize -> 0.1, Frame -> True, 
  QEpilog -> gr]

QFourierSumPlot[GaussListOfPlaneWaves[3 Pi], PlotRange->All]

QFourierSumPlot[GaussListOfPlaneWavesShift[2*Pi], 
  QXRegion -> {{-12, 12}, {-0.2, 1.2}}, Frame -> True]

QFourierSumPlot[GaussListOfPlaneWavesShift[5], QHalfPeriod -> 10, 
  PlotRange -> All, Frame -> True, Axes -> None, QBorders -> None]

QFourierSumPlot[GaussListOfPlaneWaves[Pi], QAmplitudeFactor -> 1, 
  PlotRange -> All]

QFourierSumPlot[GaussListOfPlaneWaves[Pi], QAmplitudeFactor -> 1, 
  QHalfPeriod -> Pi, QXRegion -> {{-3*(Pi/2), 3*(Pi/2)}, {-0.5, 6.4}}, 
  Frame -> True]

QFourierSumPlot[{{0, 1}}, QXRegion -> {{-4, 4}, {-0.2, 0.7}}]

QFourierSumPlot[{{2, 1}}, QHalfPeriod -> Pi, QXRegion -> {{-4, 4}, {-0.2, 0.9}}]

QFourierSumPlot[{{2, 1}}, QXRegion -> {{-4, 4}, {-0.2, 0.9}}]

QWaveStackPlot[exampleListOfPlaneWaves, QBorders -> False]

QWaveStackPlot[{{0, 1}}, Frame -> True, QAmplitudeFactor -> 1]

QWaveStackPlot[GaussListOfPlaneWaves[2*Pi], Frame -> True, 
  QAmplitudeFactor -> 1, QXRegion -> {{-4*Pi, 4*Pi}, {-0.5, 5.6}}]

QWaveStackPlot[GaussListOfPlaneWaves[2*Pi], Frame -> True, 
  QAmplitudeFactor -> 1, QXRegion -> {{-4*Pi, 4*Pi}, {-0.5, 5.6}}, 
  QSeparationLines -> False, Axes -> False]

QWaveStackWaveNumberPlot[sinListOfPlaneWaves]

QWaveStackWaveNumberPlot[GaussListOfPlaneWaves[Pi], QLineSize -> 0.4, 
  QAmplitudeFactor -> 1, QXRegion -> 
   {{-7, 7}, {-0.1*Sqrt[2*Pi], 1.1*Sqrt[2*Pi]}}, 
  QKRegion -> {{-7, 7}, {-0.1, 1.1}}, FrameTicks -> 
   {Automatic, {0, 1, 2}, Automatic, Automatic}, QEpilog -> gr, 
  PlotLabel -> "QWaveStackWaveNumberPlot", 
  QSPlotLabel -> "pile of plane waves", QKPlotLabel -> "wave numbers/amplitudes"]

QFourierSumWaveNumberPlot[GaussListOfPlaneWaves[Pi], 
  QXPlotLabel -> "position space"]

QFourierSumWaveNumberPlot[GaussListOfPlaneWaves[Pi], QLineSize -> 0.4, 
  QAmplitudeFactor -> 1,
  QXRegion -> {{-7, 7}, {-0.1*Sqrt[2*Pi], 1.1*Sqrt[2*Pi]}}, 
  QKRegion -> {{-7, 7}, {-0.1, 1.1}},
  FrameTicks -> {Automatic, {0, 1, 2}, Automatic, Automatic}, QEpilog -> gr, 
  PlotLabel -> StringForm["L = ``", PaddedForm[N[Pi], {6, 5}]], 
  QXPlotLabel -> "position space", QKPlotLabel -> "momentum space"]

QFourierSumWaveNumberPlot[{{0, 1}}, QHalfPeriod -> Sqrt[2*Pi], 
  QXRegion -> {{-7, 7}, {-0.2, 1.5}}, QKRegion -> {{-7, 7}, {-0.2, 1.5}}, 
  FrameTicks -> {Automatic, {0, 1}, Automatic, Automatic}]

QFourierSumWaveNumberPlot[{{0, 1}}, QHalfPeriod -> Sqrt[2*Pi]/2, 
  QXRegion -> {{-7, 7}, {-0.2, 1.5}}, QKRegion -> {{-7, 7}, {-0.2, 1.5}}, 
  FrameTicks -> {Automatic, {0, 1}, Automatic, Automatic}, AxesLabel -> None]

QFourierSumWaveNumberPlot[{{0, 1}}, QHalfPeriod -> Sqrt[2*Pi]/2, 
  QXRegion -> {{-7, 7}, {-0.2, 1.5}}, QKRegion -> {{-7, 7}, {-0.2, 1.5}}, 
  FrameTicks -> {Automatic, {0, 1}, Automatic, Automatic}, AxesLabel -> None, 
  QEpilog -> gr]

QFourierSumWaveNumberPlot[{{0, 1}}, QHalfPeriod -> Sqrt[2*Pi]/2, 
  QXRegion -> {{-7, 7}, {-0.2, 1.5}}, QKRegion -> {{-7, 7}, {-0.2, 1.5}}, 
  FrameTicks -> {Automatic, {0, 1}, Automatic, Automatic}, AxesLabel -> None, 
  Epilog -> gr]

QStackWavesGraph[sinListOfPlaneWaves, QHalfPeriod -> 2*Pi]

QStackWavesGraph[GaussListOfPlaneWaves[3*Pi], QSeparationLines -> False]

QStackWavesGraph[exampleListOfPlaneWaves, QSeparationLines -> False, 
  QBorders -> False]

QFourierSumWaveStackPlot[GaussListOfPlaneWaves[Pi], QSeparationLines -> False]

QFourierSumWaveStackPlot[GaussListOfPlaneWaves[Pi], QAmplitudeFactor -> 1, 
  QXRegion -> {{-7, 7}, {-0.1*Sqrt[2*Pi], 1.1*Sqrt[2*Pi]}}, 
  FrameTicks -> {Automatic, {0, 1, 2}, Automatic, Automatic}, 
  PlotLabel -> StringForm["L = ``", PaddedForm[N[Pi], {6, 5}]], 
  QXPlotLabel -> "position space", QSPlotLabel -> "plane waves"]

QVisualizeSuperposition[GaussListOfPlaneWaves[Pi], QLineSize -> 0.4, 
  QAmplitudeFactor -> 1, QXRegion -> 
   {{-7, 7}, {-0.1*Sqrt[2*Pi], 1.1*Sqrt[2*Pi]}}, 
  QKRegion -> {{-7, 7}, {-0.1, 1.1}}, FrameTicks -> 
   {Automatic, {0, 1, 2}, Automatic, Automatic}, QEpilog -> gr, 
  PlotLabel -> StringForm["L = ``", PaddedForm[N[Pi], {6, 5}]], 
  QXPlotLabel -> "position space", QKPlotLabel -> "wave number space", 
  QSPlotLabel -> "stacked waves"]

QVisualizeSuperposition[GaussListOfPlaneWaves[2*Pi], QEpilog -> gr, 
  QKRegion -> {{-7, 7}, {-0.2, 1.2}}, QXRegion -> {{-7, 7}, {-0.2, 1.2}}, 
  QSeparationLines -> None]

QVisualizeSuperposition[GaussListOfPlaneWaves[2*(Pi/3)], QEpilog -> gr, 
  QKRegion -> {{-7, 7}, {-0.2, 1.2}}, QXRegion -> {{-7, 7}, {-0.2, 1.2}}, 
  QSeparationLines -> None, FrameTicks -> None, Ticks -> None, 
  AxesLabel -> None]
 
*)

