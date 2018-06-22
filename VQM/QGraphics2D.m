(* :Title:	Quantum Kernel Plot *)

(* :Name:	VQM`QGraphics2D` *)

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
*)

(* :Summary: 
This package provides some plot commands for use together
with the QuantumKernel application.
*)

(* :Date:	2006-07-20 *)

(* :Package Version:		2.0 *)

(* :Mathematica Version:	6.0.1 *)

(* :Keywords:
	Plot, QuantumKernel, ComplexPlot
*)
   
(*-----------------------------------*)
BeginPackage[
    "VQM`QGraphics2D`",
    {"VQM`QuantumKernel`",
    "VQM`ComplexPlot`",
    "VQM`ColorMaps`"}]  ;              
(*-----------------------------------*)
    

VQM`QGraphics2D`Private`Symbols = Hold[
	QPrepareOptions, QExtractPart, QGetSpinorAndDensityPlot,
	QMakeTable, QZeroTable, QParameters, QGetSpinorAndDensityPlotTwo,
	QGetAndSpinorToColorPlot, QGetAndSpinorToColorPlotTwo
];

Unprotect @@ VQM`QGraphics2D`Private`Symbols;
ClearAll @@ VQM`QGraphics2D`Private`Symbols;
    
QPrepareOptions::usage = "QPrepareOptions[{dx,numLeft,numRight,plotLeft,plotRight,skipFac}]
	turns a list of parameters into a list of Options for QGraphics2D commands.
	Auxiliary function. See description of QParameters. Package: VQM`QGraphics2D`.";

QExtractPart::usage = "QExtractPart[qparams]
	generates part specifications. It returns lists of indices suitable
	for use with the Part function. Here qparams = {dx,numLeft,numRight,plotLeft,plotRight,skipFac}.
	QExtractPart is useful for extracting a part of a large two-dimensional numerical array
	of numbers. Assume that the array describes function values on a fine grid of
	space points in the region defined by numLeft and numRight and grid constant dx. But you want to
	plot only the values in region defined by plotLeft and plotRight. If you want to plot at
	a lower resolution, you may want to keep only every n-th
	value in the x-direction and only every m-th value in the y-direction. Then choose
	skipfac = {n,m}. The smaller array of numbers containing only the values to be plotted
	is then obtained by Part[array, QExtractPart[qparams]]. See also the description of QParameters.
	Package: VQM`QGraphics2D`.";

QGetAndDensityPlot::usage = "QGetAndDensityPlot[psi,T,QParameters->qparams,opts] is a utility function for
	visualizing a numerically determined function psi. It is assumed that psi is a function object
	defined in QuantumKernel (see VQM`QuantumKernel`). Via MathLink, QGetAndDensityPlot gets the numerical array
	of complex numbers representing psi from QuantumKernel. Then it extracts the part of psi that is needed
	for the visualization, as specified by qparams = {dx,numLeft,numRight,plotLeft,plotRight,skipFac}. See
	the description of QParameters. The option QParameters must be given.
	Finally, QGetAndDensityPlot produces a density plot of Abs[psi] with PlotLabel t=T.
	Package: VQM`QGraphics2D`.";

QGetAndComplexDensityPlot::usage = "QGetAndComplexDensityPlot[psi,T,QParameters->qparams,opts] is a utility function for
	visualizing a numerically determined function psi. It is assumed that psi is a function object
	defined in QuantumKernel (see VQM`QuantumKernel`). Via MathLink, QGetAndComplexDensityPlot gets the numerical array
	of complex numbers representing psi from QuantumKernel. Then it extracts the part of psi that is needed
	for the visualization, as specified by qparams = {dx,numLeft,numRight,plotLeft,plotRight,skipFac}. See
	the description of QParameters. The option QParameters must be given.
	Finally, QGetAndComplexDensityPlot produces a colored density plot of psi with PlotLabel t=T.
	Package: VQM`QGraphics2D`.";

QGetSpinorAndDensityPlot::usage = "QGetSpinorAndDensityPlot[psi,T,QParameters->qparams,opts] is a utility function for
	visualizing a numerically determined spinor psi.
	It is assumed that QuantumKernel computes a spinor psi
	(given by a two-dimensional array of 4 real numbers representing real and imaginary parts
	of upper and lower components).
	QGetSpinorAndDensityPlot extracts (via MathLink) the array psi of spinors from QuantumKernel
	and visualizes the absolute value with a density plot (grayscale image).
	psi is the name of the function object defined in QuantumKernel.
	T is the time variable for the PlotLabel.
	The required option QParameters specifies the parameters for extracting from psi
	the values that are needed for the visualization.
	See also the description of QParameters.  Package: VQM`QGraphics2D`.";

QMakeTable::usage = "QMakeTable[f,QParameters->qparams] turns a function f into an array of numerical values
	by computing the values of f on a two-dimensional grid of points as specified by qparams.
	The option QParameters must be given. Here qparams = {dx, numleft, numright, ...} describes the
	numerical region and the spacing of grid points. See also the description of QParameters.
	The numerical array can then be passed to QuantumKernel. Package: VQM`QGraphics2D`.";

QZeroTable::usage = "QZeroTable[QParameters->qparams] generates a table of zero values matching the two-dimensional
	grid defined by qparams = {dx, numleft, numright, ...}. See also the description of QParameters.
	Package: VQM`QGraphics2D`.";

QParameters::usage = "QParameters is an option for the plot commands defined in the
	package QGraphics2D.m. The form is
	QParameters->{dx,numleft,numright,plotleft,plotright,skipfac}.
	Here dx is the step size in space (optionally dx={dx1,dx2}).
	numleft = {xl,yl} and numright = {xr,yr} are the borders
	of the region where the numerical computation is done.
	plotleft = {pxl,pyl} and plotright = {pxr,pyr} are the borders of the plot region, which is
	usually smaller than the numerical region. skipfac describes which data points
	actually become pixels in the plot. It indicates which
	values are plotted and which are dropped in order to make the final graphics smaller.
	For example, skipfac = {3,2} means that only every third value in the x-direction
	and only every second value in the y-direction is actually plotted.
	Package: VQM`QGraphics2D`.";

QGetSpinorAndDensityPlotTwo::usage = "QGetSpinorAndDensityPlotTwo[psiUp,psiDown,T,QParameters->qparams,opts]
	is a utility function for visualizing a numerically determined spinor.
	It is assumed that QuantumKernel computes two complex functions psiUp, psiDown,
	representing upper and lower components of the spinor.
	QGetSpinorAndDensityPlotTwo extracts (via MathLink) these arrays from QuantumKernel
	combines them into a spinor and visualizes the absolute value with a density plot (grayscale image).
	psiUp, psiDown are the names of the wave function objects defined in QuantumKernel.
	T is the time variable for the PlotLabel.
	The required option QParameters specifies the parameters for extracting from psi the values that are
	needed for the visualization. See also the description of QParameters.  Package: VQM`QGraphics2D`.";

QGetAndSpinorToColorPlot::usage = "QGetAndSpinorToColorPlot[psi,T,QParameters->qparams,opts]
	is a utility function for visualizing a numerically determined spinor psi.
	It is assumed that QuantumKernel computes a spinor psi
	(given by a two-dimensional array of 4 real numbers representing real and imaginary parts
	of upper and lower components).
	QGetAndSpinorToColorPlot extracts (via MathLink) the array spinor from QuantumKernel
	and visualizes it by associating a color to the local spin direction.
	psi is the name of the function object defined in QuantumKernel.
	T is the time variable for the PlotLabel.
	The required option QParameters specifies the parameters for extracting from psi
	the values that are needed for the visualization.
	See also the description of QParameters.  Package: VQM`QGraphics2D`.";

QGetAndSpinorToColorPlotTwo::usage = "QGetAndSpinorToColorPlotTwo[psiUp,psiDown,T,QParameters->qparams,opts]
	is a utility function for visualizing a numerically determined spinor.
	It is assumed that QuantumKernel computes two complex functions psiUp, psiDown,
	representing upper and lower components of the spinor.
	QGetSpinorAndDensityPlotTwo extracts (via MathLink) these arrays from QuantumKernel
	combines them into a spinor and visualizes it by associating a color to the local spin-direction (color array plot).
	psiUp, psiDown are the names of the wave function objects defined in QuantumKernel.
	T is the time variable for the PlotLabel.
	The required option QParameters specifies the parameters for extracting from psi the values that are
	needed for the visualization. See also the description of QParameters.  Package: VQM`QGraphics2D`.";


(*-----------------------------------*)
Begin["`Private`"];
(*-----------------------------------*)
 

Options[QPlotCommands] = {QParameters->None};


QPrepareOptions[qparams_]:=
	Module[{dx,numLeft,numRight,plotLeft,plotRight,skipFac},
		{dx,numLeft,numRight,plotLeft,plotRight,skipFac}=qparams;
		Sequence@@{Mesh->False, Background->GrayLevel[0],
			DataRange->({plotLeft,plotRight}//Transpose),
			AspectRatio->(plotRight[[2]]-plotLeft[[2]])/(plotRight[[1]]-
                plotLeft[[1]])}
  ];


QExtractPart[{dx_, numLeft_, numRight_, plotLeft_, plotRight_, skipFac_}] := 
	Module[{a = Floor[(plotLeft - numLeft)/dx + {1, 1}], 
		b = (plotRight - numLeft)/dx, skipF},
		skipF=If[VectorQ[skipFac]===False,{skipFac,skipFac},skipFac,skipFac];
		Sequence @@ {Table[i, {i, a[[2]], b[[2]], skipF[[2]]}], 
			Table[i, {i, a[[1]], b[[1]], skipF[[1]]}]}
  ];


QMakeTable[f_, opts___Rule]:= 
	Module[{qparams,dx,numLeft,numRight},
		qparams = QParameters/.{opts}/.Options[QPlotCommands];
		{dx,numLeft,numRight} = Take[qparams,3];
		dx=If[VectorQ[dx]===False,{dx,dx},dx,dx];
		Table[f[x,y],
			{y,numLeft[[2]]+dx[[2]]/2,numRight[[2]]-dx[[2]]/2,dx[[2]]},
    		{x,numLeft[[1]]+dx[[1]]/2,numRight[[1]]-dx[[1]]/2,dx[[1]]}
    	]
    ];


QZeroTable[opts___Rule] :=
	Module[{qparams,dx,numLeft,numRight},
	qparams = QParameters/.{opts}/.Options[QPlotCommands];
	{dx,numLeft,numRight} = Take[qparams,3];
    Table[0.,{y,numLeft[[2]]+dx/2,numRight[[2]]-dx/2,dx},
    	{x,numLeft[[1]]+dx/2,numRight[[1]]-dx/2,dx}]
    ];


QGetAndDensityPlot[psi_, T_, opts___] :=
	Module[{psi1 = QGetArray[psi], psi3, qparams},
		psi3 = psi1[[1]]^2 + psi1[[2]]^2;
		qparams = QParameters/.{opts}/.Options[QPlotCommands];
		If[qparams=!=None, psi3 = Part[psi3,QExtractPart[qparams]]];
		Remove[psi1];
		ListDensityPlot[psi3, FilterRules[Flatten[{opts}], Options @ ListDensityPlot],
			QPrepareOptions[qparams],
			PlotLabel->StringForm["t =``", PaddedForm[T,{3,2}]]
		]
	];


QGetAndComplexDensityPlot[psi_, T_, opts___] :=
	Module[{psi1 = QGetArray[psi], psi3, qparams},
		psi3 = psi1[[1]] + I*psi1[[2]];
		qparams = QParameters/.{opts}/.Options[QPlotCommands];
		If[qparams=!=None, psi3 = Part[psi3,QExtractPart[qparams]]];
		Remove[psi1];
		QListComplexDensityPlot[psi3, FilterRules[Flatten[{opts}], Options @ ListDensityPlot],
			QPrepareOptions[qparams],
			PlotLabel->StringForm["t =``", PaddedForm[T,{3,2}]]
		]
	];


QGetSpinorAndDensityPlot[psi_, T_, opts___] :=
	Module[{psi1 = QGetArray[psi], psi3, qparams},
		psi3 = psi1[[1]]^2 + psi1[[2]]^2 + psi1[[3]]^2 + psi1[[4]]^2;
		qparams = QParameters/.{opts}/.Options[QPlotCommands];
		If[qparams=!=None, psi3 = Part[psi3,QExtractPart[qparams]]];
		Remove[psi1]; (* because of memory management efficiently *)
		ListDensityPlot[psi3, FilterRules[Flatten[{opts}], Options @ ListDensityPlot],
			QPrepareOptions[qparams],
			PlotLabel->StringForm["t =``", PaddedForm[T,{3,2}]]
		]
	];


QGetSpinorAndDensityPlotTwo[psiUp_, psiDown_, T_, opts___] :=
	Module[{psi1 = QGetArray[psiUp],
			psi2 = QGetArray[psiDown], psi3, qparams},
		psi3 = psi1[[1]]^2 + psi1[[2]]^2 + psi2[[1]]^2 + psi2[[2]]^2;
		qparams = QParameters/.{opts}/.Options[QPlotCommands];
		If[qparams=!=None, psi3 = Part[psi3,QExtractPart[qparams]]];
		Remove[psi1,psi2];
		ListDensityPlot[psi3, FilterRules[Flatten[{opts}], Options @ ListDensityPlot],
			QPrepareOptions[qparams],
			PlotLabel->StringForm["t =``", PaddedForm[T,{3,2}]]
		]
	];


QGetAndSpinorToColorPlot[psi_, T_, opts___] :=
	Module[{psi1 = QGetArray[psi], psi3, qparams}, 
    	qparams = QParameters /. {opts} /. Options[QPlotCommands];
		If[qparams === None, Return[]]; 
		psi3 = Transpose[
			{Transpose[(psi1[[1]] + I*psi1[[2]])[[QExtractPart[qparams]]]], 
        	 Transpose[(psi1[[3]] + I*psi1[[4]])[[QExtractPart[qparams]]]]}, {3, 2, 1}];
        Remove[psi1]; 
		QColorArrayPlot[Apply[RGBColor, Chop[Map[QSpinorToColor, psi3, {2}]], {2}], 
			FilterRules[Flatten[{opts}], Options @ QColorArrayPlot],
			QPrepareOptions[qparams], 
			PlotLabel -> StringForm["t =``", PaddedForm[T, {3, 2}]]
		]
	]; 


QGetAndSpinorToColorPlotTwo[psiUp_, psiDown_, T_, opts___] :=
	Module[{psi1 = QGetArray[psiUp],
			psi2 = QGetArray[psiDown],
			psi3, qparams},
		qparams = QParameters/.{opts}/.Options[QPlotCommands];
		If[qparams === None, Return[]];
		psi3 = Transpose[{Transpose[Part[psi1[[1]]+I psi1[[2]],QExtractPart[qparams]]],
    					  Transpose[Part[psi2[[1]]+I psi2[[2]],QExtractPart[qparams]]]},{3,2,1}];
		Remove[psi1,psi2];
		QColorArrayPlot[Apply[RGBColor,Chop[Map[QSpinorToColor,psi3,{2}]],{2}],
			FilterRules[Flatten[{opts}], Options @ QColorArrayPlot],
			QPrepareOptions[qparams],
			PlotLabel->StringForm["t =``", PaddedForm[T,{3,2}]]
		]
	];


(*-----------------------------------*)
End[];      (* end `Private` Context *)
(*-----------------------------------*)

Protect@@VQM`QGraphics2D`Private`Symbols;

(*-----------------------------------*)
EndPackage[]; (* end package Context *)
(*-----------------------------------*)

