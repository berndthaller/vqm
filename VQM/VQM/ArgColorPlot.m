(* ::Package:: *)

(* :Title:	Arg Color Plot *)

(* :Name:	VQM`ArgColorPlot` *)

(* :Copyright: Copyright 2007 by Bernd Thaller *)

(* :Author:	Bernd Thaller,
			Institute of Mathematics,
			University of Graz,
			A-8010 Graz
			bernd.thaller@uni-graz.at
*)

(* :Source:
	(Advanced) Visual Quantum Mechanics
	Springer-Verlag New York, 2004
	see "Book One" for a description of the color map
*)

(* :Summary: 
	Plot the absolute value Abs[f[x]] of a complex-valued
	function f depending on a real variable x and
	fill the space between the plotted function
	and the x-axis with a color corresponding to the argument
	Arg[f[x]]. The saturation and brightness
	of the colors can be set using the options QSaturation
	and QBrightness.
*)

(* :Date:	2007-07-20 *)

(* :Package Version:		2.0.0 *)

(* :Mathematica Version:	6.0.1 *)

(* :Keywords:
	Plot, Area under the curve, ArgColors, Wavefunction
*)
    
VQMmsgon  = Head[General::"spell"]  =!= $Off;
VQMmsgon1 = Head[General::"spell1"] =!= $Off;

Off[General::spell1, General::spell];
    

(*-----------------------------------*)
BeginPackage[
    "VQM`ArgColorPlot`",
    "Utilities`FilterOptions`"
];
(*-----------------------------------*)


VQM`ArgColorPlot`Private`Symbols = Hold[
    QArgColorPlot,
    QListArgColorPlot, QCombinedPlot,
    QListCombinedPlot, QSpinorPlot, QListSpinorPlot, 
    QSpinorCombinedPlot, QListSpinorCombinedPlot,
    QNiceTicks, QSaturation, QBrightness, QBottomLine,
    QShiftPlot, QHorizontalRange, QPlotDown, QSquared, QCurveStyle];

Unprotect @@ VQM`ArgColorPlot`Private`Symbols;
ClearAll @@ VQM`ArgColorPlot`Private`Symbols;

(* 
this is a local flag to whether to use ListLinePlot or not,
could be exported if needed here.
Default is True.
$UseListLinePlot::usage="test (RM)";
*)

QArgColorPlot::usage = "QArgColorPlot[f[x],{x,x0,x1},opts] is used
like the usual Plot command. It gives a two-dimensional plot
of a complex-valued function f of a single real variable x
in the range {x0,x1}.
The plot shows the curve Abs[f] with area between the curve
and the x-axis colored by Hue[Arg[f[x]]/(2 Pi)].
The default options of Plot are changed to Axes->{True,None},
Frame->True. Package: VQM`ArgColorPlot`";

QListArgColorPlot::usage = "QListArgColorPlot[f,{x,x0,x1},opts]
plots a Abs[f], where f is a list of complex numbers. The 
points of the list Abs[f] are joined by a line.
The area between the curve and the x-axis is colored at each
point by Hue[Arg[f]/(2 Pi)]. Package: VQM`ArgColorPlot`";

QCombinedPlot::usage = "QCombinedPlot[{f[x],g[x]},{x,x0,x1},opts]
works like QArgColorPlot with respect to f. The curve g is drawn
in front of the QArgColorPlot of f. Package: VQM`ArgColorPlot`";

QListCombinedPlot::usage = "QListCombinedPlot[{list,f[x]},{x,x0,x1},opts]
works like QListArgColorPlot with respect to list.
It is assumed that list represents the discretized values of a function
defined on the interval [x0,x1]. The color list plot
is then combined with an ordinary plot of f on the same scale and
with the Ticks automatically adjusted. Package: VQM`ArgColorPlot`";

QSpinorPlot::usage = "QSpinorPlot[{func1,func2},{x,x0,x1},opts] provides a
method to visualize C^2-valued functions (for example, spinor wavefunctions
in quantum mechanics). The QSpinorPlot combines a QArgColorPlot of func1
with a QArgColorPlot of func2 (upside down, with less saturation)
Both curves are plotted with the option QSquared->True (that is, a plot
of the curve Abs[func]^2 is filled with a color describing the phase).
In the background, a filled plot of Abs[func1]^2 + Abs[func2]^2
displays the corresponding density. Package: VQM`ArgColorPlot`";

QListSpinorPlot::usage = "QListSpinorPlot[list,opts] visualizes
a spinor-valued list of complex numbers. Each element of list is
a C^-vector, that is, list = {{z11,z12},{z21,z22},...}.
Alternatively, list = {list1,list2} with two lists of complex numbers,
list1 giving the upper component of the spinor-valued wave function, and
list2 giving the lower component. The lower component is plotted
upside down with less saturation. See also the description of QSpinorPlot.
Package: VQM`ArgColorPlot`";

QSpinorCombinedPlot::usage =
"QSpinorCombinedPlot[{func1,func2},{x,x0,x1},opts] combines
a QSpinorPlot of func1 with an ordinary Plot of a real-valued function func2.
See the description of QCombinedPlot and
of QSpinorPlot. Package: VQM`ArgColorPlot`";

QListSpinorCombinedPlot::usage =
"QListSpinorCombinedPlot[{list,f[x]},{x,x0,x1},opts] combines
a QListSpinorPlot of list1 with an ordinary Plot  of a real-valued function f.
See the description of QListCombinedPlot and of QListSpinorPlot.
Package: VQM`ArgColorPlot`";

QNiceTicks::usage = "QNiceTicks[xmin,xmax,dx] provides a list of nice positions
for use in the Ticks or FrameTicks option in a ListPlot, where it is assumed
that the list of values ranges between xmin and xmax in steps dx. Package: VQM`ArgColorPlot`";

QSaturation::usage = "Option for QArgColorPlot. QSaturation->s causes the colors
in the plot to appear at saturation s. Default value is 1. Package: VQM`ArgColorPlot`";

QBrightness::usage = "Option for QArgColorPlot. QBrightness->b causes the colors
to appear with brightness b. Default value is 1. Package: VQM`ArgColorPlot`";

QBottomLine::usage = "Option for QArgColorPlot.
QBottomLine->number fills the region between the graph of
the function and the horizontal line at y=number. Package: VQM`ArgColorPlot`";

QShiftPlot::usage = "Option for QArgColorPlot. QShiftPlot->number performs a vertical
shift of the whole plot by number. Package: VQM`ArgColorPlot`";

QHorizontalRange::usage = "Option for QListArgColorPlot that describes the range of
coordinates on the horizontal axis. QHorizontalRange->{{x1,x2},{y1,y2}} means that
the list represents values given on a regular grid of points in
the interval [x1,x2], but only the values in the interval [y1,y2] are to be plotted.
The interval [y1,y2] is the visible coordinate range
on the horizontal axis. Ticks for axes or frame are defined automatically.
If [y1,y2] is larger than [x1,x2], then the list is padded with zeros.
QHorizontalRange->{x1,x2} can be used as a shortcut for QHorizontalRange->{{x1,x2},{x1,x2}}.
Package: VQM`ArgColorPlot`";

QPlotDown::usage = "Option for QArgColorPlot. If set to True, the plot is 'upside down',
that is, the graph is drawn in the negative y-direction Package: VQM`ArgColorPlot`";

QSquared::usage = "Option for QArgColorPlot. If set to True, then the graph shows
the absolute square of the function, otherwise the graph shows the
absolute value. Package: VQM`ArgColorPlot`";

QCurveStyle::usage = "Option for QCombinedPlot. Controls how the curve representing
the real function is drawn. Works like PlotStyle, which affects only
the complex part. Package: VQM`ArgColorPlot`";


(*-----------------------------------*)
Begin["`Private`"];
(*-----------------------------------*)


If[!ValueQ[$UseListLinePlot],
$UseListLinePlot=True;
];

Options[QArgColorPlot] := Sort @ {JoinOptions[
        {
        PlotStyle -> None
        },
        Options[ If[$UseListLinePlot, ListLinePlot, Plot] ],
		{
        Opacity -> 1,         
        Compiled->True, 
		QSaturation->1, QBrightness->1, QBottomLine->0, QShiftPlot->0, QPlotDown->False, 
		QSquared->False
        } ] };

Options[QListArgColorPlot] =  Sort @  Join[Options[QArgColorPlot],{QHorizontalRange->All}];

SetAttributes[QArgColorPlot,HoldAll];

QArgColorPlot[func_,{x_Symbol,xmin_,xmax_},opts___?OptionQ] :=
	Module[{comp,plot1,opts1,xvars,yvals},
        comp = Compiled/.Join[{opts}, Options[QArgColorPlot]];
        If[comp === True,
          abf = Compile[{x},Abs[func],{{func,_Complex}}];
          fnc = Compile[{x},func,{{func,_Complex}}], (*else*)
          abf = Abs[func];
          fnc = Function[{x},func];
        ];
        opts1   = FilterOptions[Plot, opts];
		plot1   = Plot @@ { abf, {x,xmin,xmax}, 
                            {opts1} /. ( PlotStyle -> _ ) :> Sequence[]};
		xvars   = First /@ Level[plot1[[1]], {5}];
        yvals   = fnc /@ xvars; 

      QListArgColorPlot[{xvars,yvals}//Transpose, opts]
	]/;NumericQ[xmin] && NumericQ[xmax]


QListArgColorPlot[list_List,opts___?OptionQ] :=
	Module[{ran, auxopts},
		ran = QHorizontalRange/.Join[{opts}, Options[QListArgColorPlot]];
		If[MatchQ[N[ran],{_?NumberQ,_?NumberQ}],ran = {ran,ran}];
		If[MatchQ[N[ran],{{_?NumberQ,_?NumberQ},{_?NumberQ,_?NumberQ}}],
			xvars = Range[ran[[1,1]],ran[[1,2]],(ran[[1,2]]-ran[[1,1]])/(Length[list]-1)];
			auxopts = JoinOptions[QHorizontalRange->{ran[[2,1]],ran[[2,2]]},opts];,
			xvars	= Range[Length[list]];
			auxopts = opts
		];
		QListArgColorPlot[{xvars,list}//Transpose, auxopts]
	]/;VectorQ[list]

QListArgColorPlot[list:{ {_, _}.. }, opts___?OptionQ] :=
	Module[{opac,interp, colorfun, res,  sat,bri,style,dir,squ,bl,nbl,sh,shift,ran,
			xvars = list[[All,1]],
			yvals = list[[All,-1]],
			hues, auxlist=list,
			optsWithDefaults
           },
         optsWithDefaults = JoinOptions[{opts}, {Axes->True}, Options[QlistArgColorPlot]];

		If[Not[And @@ (NumberQ[#]& /@ Flatten[N[list]])],
			Message[LACP::listform]; Return[] ];
		    {opac, sat,bri,style,dir,squ,bl,sh,ran} =
			{Opacity, QSaturation,QBrightness,PlotStyle,QPlotDown,QSquared,
			 QBottomLine,QShiftPlot,QHorizontalRange}/.{opts}/.Options[QListArgColorPlot];
        nbl    = If[NumberQ[N[bl]],N[bl],0.,0.];
        shift  = If[NumberQ[N[sh]],N[sh],0.,0.];
		If[MatchQ[N[ran],{{_?NumberQ,_?NumberQ},{_?NumberQ,_?NumberQ}}],
			ran = ran[[2]] ];
		If[MatchQ[N[ran],{_?NumberQ,_?NumberQ}],
		    auxlist = TrimList[list,{ran[[1]],ran[[2]]}];
(* RM: This is necessary for successful Interpolation later *)
            auxlist = Union[Chop @ auxlist];  
			xvars    = First /@ auxlist;
            yvals	= Last /@ auxlist; 
         ];
         
         yvals	= If[squ === True, Abs[yvals]^2, Abs[yvals], Abs[yvals]];
         yvals	= If[dir === True, -1*yvals, yvals, yvals];

If[$UseListLinePlot,
         (* create a complex-valued interpolation function *)
         interp   = Interpolation[auxlist];
         If[opac == 1, opac=Sequence[]];
         colorfun = Function[{x,y},(Hue[#,sat,bri,opac]& @ (Mod[Arg[interp[x]]/(2Pi), 1]))];
         res = fillit[xvars,colorfun,yvals,style,nbl,shift, optsWithDefaults];
    (* else *)
         ,
         hues = Hue[Mod[Arg[#]/(2 Pi),1],sat,bri, opac]& /@ yvals;
         res = fillit[xvars,hues,yvals,style,nbl,shift]
];
   fiops = FilterOptions[Graphics,optsWithDefaults];
   Show[res, fiops,AspectRatio->(1/GoldenRatio)]
];

QListArgColorPlot[sillyArgument_,opts___] :=
    Message[LACP::listform]

SetAttributes[QSpinorPlot,HoldAll];

QSpinorPlot[{func1_,func2_},{x_Symbol,xmin_,xmax_},opts___] :=
	Module[{},
		Show[
			QArgColorPlot[Sqrt[Abs[func1]^2 + Abs[func2]^2], {x,xmin,xmax},
				QSaturation->0, QBrightness->1/4, opts],
			QArgColorPlot[func1, {x,xmin,xmax}, opts],
			QArgColorPlot[func2, {x,xmin,xmax},
				QPlotDown->True, opts, QSaturation->1/2],
			Evaluate[FilterOptions[Graphics,opts]]
  		]
  	]/;NumberQ[N[xmin]] && NumberQ[N[xmax]]

QListSpinorPlot[{list1_List,list2_List},opts___] :=
	Module[{},
		If[Length[list1]!=Length[list2],Message[LSP::unequal];Return[]];
		If[!VectorQ[list1],Message[LSP::listform];Return[]];
		If[!VectorQ[list2],Message[LSP::listform];Return[]];
		Show[
			QListArgColorPlot[Sqrt[Abs[list1]^2 + Abs[list2]^2],
				QSaturation->0, QBrightness->1/4, opts],
			QListArgColorPlot[list1, opts],
			QListArgColorPlot[list2,
				QPlotDown->True, opts, QSaturation->1/2],
			Evaluate[FilterOptions[Graphics,opts]]
		]
	]

QListSpinorPlot[list_, opts___?OptionQ] := 
	QListSpinorPlot[list//Transpose, opts
	]/;((Dimensions[list] == {Length[list],2}) && Length[list]>2 && VectorQ[Last/@list])

QListSpinorPlot[list_,opts___?OptionQ] := 
		Show[
			QListArgColorPlot[{First /@ list,
					Sqrt[Abs[First /@ Last /@ list]^2 +
					Abs[Last  /@ Last /@ list]^2]}//Transpose,
				QSaturation->0, QBrightness->1/4, opts],
			QListArgColorPlot[{First /@ list,First /@ Last /@ list}//Transpose, opts],
			QListArgColorPlot[{First /@ list,Last  /@ Last /@ list}//Transpose,
				QPlotDown->True, opts, QSaturation->1/2],
			Evaluate[FilterOptions[Graphics,opts]]
		]/;MatchQ[ N[list],{{_?NumberQ,{_?NumberQ,_?NumberQ}},__}]

LSP::unequal = "The two lists of complex numbers must be of equal length";
LACP::listform = "The first argument of QListArgColorPlot must be a one-dimensional list
of complex numbers zi, or a list of the form {{x1,z1},{x2,z2},...}";
LSP::listform = "One of the arguments is not a one-dimensional list of complex numbers";

Options[QCombinedPlot] = Join[Options[QArgColorPlot], {QCurveStyle->Automatic,
                                                       PlotRange -> Automatic}];
Options[QListCombinedPlot] = Join[Options[QListArgColorPlot], {QCurveStyle->Automatic}];

QCombinedPlot[{func1_,func2_},
             {x_Symbol,xmin_,xmax_},
             opts___Rule] :=
	Module[{ps = QCurveStyle/.{opts}/.Options[QCombinedPlot], auxopts},
		auxopts = JoinOptions[PlotStyle->ps, opts];
		Show[
			QArgColorPlot[func1, {x,xmin,xmax},
				Evaluate[FilterOptions[QArgColorPlot,opts]] ],
			Plot[func2, {x,xmin,xmax},
				Evaluate[FilterOptions[Plot,auxopts]] ],
            PlotRange -> Automatic,
			Evaluate[FilterOptions[QCombinedPlot,opts]]]
     ] /; NumberQ[N[xmin]] && NumberQ[N[xmax]]


QListCombinedPlot[{list_List,func_},
	{x_Symbol,xmin_,xmax_}, opts___Rule] :=
Module[{ran,ps,auxopts},
		{ran,ps} = {QHorizontalRange,QCurveStyle}/.{opts}/.Options[QListCombinedPlot];
		If[MatchQ[N[ran],{_?NumberQ,_?NumberQ}],
			ran = {ran,{xmin,xmax}}];
		If[Not[MatchQ[N[ran],{{_?NumberQ,_?NumberQ},{_?NumberQ,_?NumberQ}}]],
			ran = {{xmin,xmax},{xmin,xmax}} ];
		auxopts = JoinOptions[PlotStyle->ps, opts];
		Show[
          QListArgColorPlot[list,
            	Evaluate[FilterOptions[QListArgColorPlot,
            		JoinOptions[QHorizontalRange->ran,opts]]] ],          
            Plot[func, {x,xmin,xmax},
				Evaluate[FilterOptions[Plot,auxopts]] ],
                PlotRange -> Automatic,
		    Evaluate[FilterOptions[Graphics,opts]]
		]
   ]/;NumberQ[N[xmin]] && NumberQ[N[xmax]]

QSpinorCombinedPlot[{func1_,func2_},
             {x_Symbol,xmin_,xmax_},
             opts___Rule] :=
	Module[{ps = QCurveStyle/.{opts}/.Options[QCombinedPlot], auxopts},
		auxopts = JoinOptions[PlotStyle->ps, opts];
		Show[
			QSpinorPlot[func1, {x,xmin,xmax},
				Evaluate[FilterOptions[QArgColorPlot,opts]] ],
			Plot[func2, {x,xmin,xmax},
				Evaluate[FilterOptions[Plot,auxopts]] ],
			Evaluate[FilterOptions[Graphics,opts]]]
     ] /; NumberQ[N[xmin]] && NumberQ[N[xmax]]

QListSpinorCombinedPlot[{list_List,func_},
	{x_Symbol,xmin_,xmax_}, opts___Rule] :=
Module[{ran,ps,auxopts},
		{ran,ps} = {QHorizontalRange,QCurveStyle}/.{opts}/.Options[QListCombinedPlot];
		If[MatchQ[N[ran],{_?NumberQ,_?NumberQ}],
			ran = {ran,{xmin,xmax}}];
		If[Not[MatchQ[N[ran],{{_?NumberQ,_?NumberQ},{_?NumberQ,_?NumberQ}}]],
			ran = {{xmin,xmax},{xmin,xmax}} ];
		auxopts = JoinOptions[PlotStyle->ps, opts];
		Show[
            QListSpinorPlot[list,
            	Evaluate[FilterOptions[QListArgColorPlot,
            		JoinOptions[QHorizontalRange->ran,opts]]] ],
            Plot[func, {x,xmin,xmax},
				Evaluate[FilterOptions[Plot,auxopts]] ],
		    Evaluate[FilterOptions[Graphics,opts]]
		]
   ]/;NumberQ[N[xmin]] && NumberQ[N[xmax]]


(* auxiliary functions *)

TrimList[list_,{xmin_,xmax_}] :=
	Module[{x1 = First[list][[1]], xn = Last[list][[1]], auxlist=list},
		auxlist = Select[list, (xmin <= First[#] <= xmax)& ];
		If[xmin < x1, auxlist = Join[{{xmin,0},{x1,0}},auxlist] ];
		If[xmax > xn, auxlist = Join[auxlist,{{xn,0},{xmax,0}}] ];
		auxlist ]/;(Dimensions[list] == {Length[list],2})

QNiceTicks[xmin_,xmax_,dx_,n_:8] :=
   convertticks[LinearScale[xmin,xmax,n],xmin,dx] /;
   NumberQ[N[xmin]] && NumberQ[N[xmax]] && NumberQ[N[dx]]

convertticks[list_,xmin_,dx_]:=
   Module[{res=list,a,myind},
          myind[x_] = ind[x,xmin,dx];
          Do[
             a = res[[i]];
             If[Head[a]===List,
                res[[i]]=
                If[Head[a[[2]]]===List,
                   res[[i]]=Join[{myind[a[[1]]],a[[1]]},Take[a,1-Length[a]]],
                   res[[i]]=ReplacePart[a,myind[a[[1]]],1]],
                res[[i]]={myind[a],a}
             ],
          {i,Length[res]}];
          res
   ]

ind[x_,xmin_,dx_] := N[((x-xmin)/dx)+1];

(* This plots the filled curve: *)

(* the original fillit from 5.2 *)
fillit[xvars_,hues_List,values_,style_,bl_,sh_] :=
   Module[{nullv,xpts,valpts,lines,fills,shvar},
      nullv = Table[bl,{Length[xvars]}];
      shvar = Table[sh,{Length[xvars]}];
      If[bl==0.,nullv=shvar];
      xpts = {xvars,nullv}//Transpose;
      valpts = {xvars,values+shvar}//Transpose;
      If[style === Automatic,
          lines = Line[valpts], (*else*)
          lines = Flatten[{style, Line[valpts]}]
      ];
      fills =
         {Drop[hues,-1],
             Map[Polygon,
               { Drop[xpts,-1], Drop[valpts,-1], Drop[valpts,1], Drop[xpts, 1]
                    }//Transpose
             ]
         }//Transpose;
      Graphics[ {fills,lines} ]
   ]

fillit[xvars_,colorfun_Function,values_,style_,bl_,sh_, opts___?OptionQ] :=
   Module[{epilogs,nullv,xpts,valpts,lines,fills,shvar,cfun},
      nullv = Table[bl,{Length[xvars]}];
      shvar = Table[sh,{Length[xvars]}];
      If[bl==0.,nullv=shvar];
      valpts = {xvars,values+shvar}//Transpose;      

      If[style === Automatic,
          lines = Line[valpts], (*else*)
          lines = Flatten[{style, Line[valpts]}]
      ];
       (* since there is a problem with multiple Epilog primitives in 6.0.1 *)
      If[style === None, lines = Sequence[]];

     (* since there might other Epilog settings from QListArgColorPlot, we merge them here *)
     epilogs = DeleteCases[{lines, Epilog /. {opts}}, Epilog];
     ListLinePlot[valpts, 
                          Filling -> (bl+sh), ColorFunctionScaling -> False,
                          Epilog ->  epilogs,   (*FillingStyle->Automatic,*) (*RM: should work with this option ... *)
                          ColorFunction ->colorfun,
                          FilterOptions[ListLinePlot, opts]
                           
     ]
	]


(* ------ dealing with options ------ *)

(* joining two lists of rules with list1 having precedence
over list2: *)

JoinOptions[list1_,list2___]:=
    Module[{namelist=First /@ Flatten[{list1}]},
        Sequence @@
        Sort[Join[Flatten[{list1}],
            Select[Flatten[{list2}],
                !MemberQ[namelist,First[#]]&]]
        ]
    ]


(* ------ copied from 5.2 Graphics`Graphics ------- *)
LinearScale::usage =
"LinearScale[xmin, xmax] gives a list of \"nice\" values between xmin and xmax
suitable for use as tick mark positions.  LinearScale[xmin, xmax, n] attempts
to find n such values.";

LinearScale[min_, max_, n_Integer:8] :=
    Module[{spacing, t, nmin=N[min], nmax=N[max]},
        (spacing = TickSpacing[nmax-nmin, n, {1, 2, 2.5, 5, 10}] ;
        t = Range[Ceiling[nmin/spacing - 0.05] spacing, max, spacing] ;
        Map[{#, If[Round[#]==#, Round[#], #]}&, t])
    /; nmin < nmax
    ];

TickSpacing[dx_, n_, prefs_List, method_:GreaterEqual] :=
    Module[ { dist=N[dx/n], scale, prefsdelta, min, pos } ,
        scale = 10.^Floor[Log[10., dist]] ;
        dist /= scale ;
        If[dist < 1, dist *= 10 ; scale /= 10] ;
        If[dist >= 10, dist /= 10 ; scale *= 10] ;
        scale * Switch[method,
            GreaterEqual,
            (* "nice" tick spacing is greater than or equal to
                requested tick spacing *)
            First[Select[prefs, (dist <= #)&]],
            LessEqual,
            (* "nice" tick spacing is less than or equal to
                                requested tick spacing *)
            First[Select[Reverse[prefs], (dist >= #)&]],
            Nearest,
            (* "nice" tick spacing is the taken from the
                element of "prefs" nearest to "dist" *)
            prefsdelta = Map[Abs[#-dist]&, prefs];
            min = Min[prefsdelta];
            pos = Position[prefsdelta, min][[1, 1]];
            prefs[[pos]]
        ]
    ];


(*-----------------------------------*)
End[];      (* end `Private` Context *)
(*-----------------------------------*)

Protect@@VQM`ArgColorPlot`Private`Symbols;

(*-----------------------------------*)
EndPackage[]; (* end package Context *)
(*-----------------------------------*)

If[VQMmsgon, On[General::"spell"]];
If[VQMmsgon1, On[General::"spell1"]];



(* :Examples:

(* Turn off error messages *)
msgonindet = Head[General::"indet"]=!=$Off;
Off[Arg::"indet" ];
		
    QArgColorPlot[Exp[-I 6 x - x^2/2], {x,-4,4}];

    QArgColorPlot[Exp[-I 6 x - x^2/2], {x,-4,4},
           PlotStyle -> {Thickness[0.025],GrayLevel[0.5]},
           QSaturation -> .5, QBrightness->.7,
           PlotRange -> {-.2,1.2}, Frame -> True,
           Axes->{True,False}]
    
    mytab = Table[Sin[Pi x] Exp[2 I Pi x], {x,-2,2,.1}];
    QListArgColorPlot[mytab, Axes->True, QHorizontalRange->{-2,2}];
    QListArgColorPlot[mytab, Axes->True, QHorizontalRange->{{-2,2},{-1,0}}];
    QListArgColorPlot[mytab, Axes->True, QHorizontalRange->{{-2,2},{-4,4}}];
    
    QCombinedPlot[{Sin[Pi x],Sin[Pi x]}, {x,-4,3}];

    mytab = Table[Sin[Pi x], {x,-4,3,.1}];
    QListCombinedPlot[{mytab,Sin[Pi x]}, {x,-4,3}];
	QListCombinedPlot[{mytab,Sin[Pi x]}, {x,0,3}];
    QListCombinedPlot[{mytab,Sin[Pi x]}, {x,0,3}, QHorizontalRange->{-4,3}];
    QListCombinedPlot[{mytab,Sin[Pi x]}, {x,0,3}, QHorizontalRange->{{-4,3},{1,2}}];
    
    Show[
	QArgColorPlot[Sin[x+I], {x,-3,3},
		QBottomLine->1, Epilog->Line[{{-3,1},{3,1}}],
		DisplayFunction->Identity],
	QArgColorPlot[Cos[x+I] - Exp[I Arg[Cos[x+I]]], {x,-3,3},
		DisplayFunction->Identity],
	PlotRange->{-.2,1.7}, Frame->True,
	DisplayFunction->$DisplayFunction]
	
	spinorfunction[x_] := {Cos[3 x] Exp[-(x-1/4)^2/3+I  x ],
      	Sin[4x]  Exp[-(x+1/4)^2/2+2 I x ]};
	QSpinorPlot[Evaluate[spinorfunction[x]], {x,-3.,3.},
		PlotRange->All, Frame->True, Axes->{True,False}]
	
	spinorlist = Table[{x,spinorfunction[x]},{x,-3,3,.1}];
	QListSpinorCombinedPlot[{spinorlist,Sin[4 x]},
		{x,-2,3}, PlotRange->All,
		QCurveStyle->{RGBColor[0,.7,.9],Thickness[0.03]},
		PlotStyle->{RGBColor[.2,.2,.5],Thickness[0.02]}]
	
	If[msgonindet, On[Arg::"indet"]]; (* turn on if it was on before *)
*)

