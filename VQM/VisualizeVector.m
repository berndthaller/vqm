(* ::Package:: *)

(* :Title: VisualizeVector *)

(* :Name: VQM`VisualizeVector` *)

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
VQM`VisualizeVector` is a package for 'Visual Quantum Mechanics'.
It helps to visualize a vector through various types of arrows in three dimensions
and defines auxiliary graphics elements.
*)

(* :Date:	2007-07-20 *)

(* :Package Version: 2.0 *)

(* :Mathematica Version: 6.0.1 *)

(* :History:
    1.0 for Visual Quantum Mechanics, Book Two, 1st ed.
    1.1 fixed bug concerning QArrowShape
    2.0 adaption to Mathematica 6
*)

VQMmsgon  = Head[General::"spell"]  =!= $Off;
VQMmsgon1 = Head[General::"spell1"] =!= $Off;

Off[General::spell1, General::spell];
    

(*-----------------------------------*)
BeginPackage[ "VQM`VisualizeVector`" ];
(*-----------------------------------*)

VQM`VisualizeVector`Private`Symbols=Hold[
	QVectorToArrow, QVisualizeVector,
	QArrowHead, QArrowShaft, QArrowShape, QArrowScale, QNeedleStyle,
	QMinLength, QLinePointSize, QHeadColor, QShaftColor,
	QDrawUnitSphere, QDrawAxes, QCoordinateCube, QCoordinateCircles, QCoordinateCirclesColor,
	polyDisk, polyCone, polyCylinder, arrowCone, arrowWithShaft, doubleHead,
	unitSphere, xaxisLine, yaxisLine, zaxisLine, axesLines,
	coordinateCube, coordinateCircles, 	
	graphicElements, QRoundSphere];

Unprotect@@VQM`VisualizeVector`Private`Symbols;
ClearAll@@VQM`VisualizeVector`Private`Symbols;
    
QVectorToArrow::usage = "QVectorToArrow[pt1,pt2] gives a collection of lines representing
	a three-dimensional arrow from point pt1 to point pt2. pt1 is optional, default is {0,0,0}.
	If the vector is shorter than QMinLength, then it is represented by a point.
	Use with Graphics3D. Example: Show[Graphics3D[QVectorToArrow[{0,0,0},{1,1,1}]]].
	The following options control the appearance of the arrow:
	QArrowHead, QArrowShaft, QArrowShape, QArrowScale, QNeedleStyle,
	QHeadColor, QShaftColor, QLinePointSize, QMinLength.
	Package: VQM`VisualizeVector`.";
	
QArrowHead::usage = "QArrowHead is an option for QVectorToArrow. QArrowHead->True (default) causes
	the vector to be drawn with an arrowhead. QArrowHead->Automatic draws the vector as a line segment
	and with a point of size 2*QLinePointSize instead of the arrowhead. Otherwise the vector is just
	represented by a line segment.
	Package: VQM`VisualizeVector`.";
	
QArrowShaft::usage = "QArrowShaft is an option for QVectorToArrow. QArrowShaft->True causes
	the vector to be drawn with a shaft. With QArrowShaft->False (default) the shaft is drawn
	as a line with thickness QLinePointSize.
	Package: VQM`VisualizeVector`.";

QArrowShape::usage = "QArrowShape->{n,hfac,rfac,sfac} is an option for QVectorToArrow that controls the
	appearance of the arrow's head. Default value of QArrowShape is {6,1/4,1/(2 GoldenRatio),1/2},
	or {6,1/2,1/4,1/2} if QNeedleStyle->True.
	The head of the arrow is a cone with half height h = hfac.length(vector) and radius r=rfac*h
	drawn using n polygons. If hfac (rfac) is negative, then h=-hfac (r=-rfac). If the option
	QArrowShaft->True is set, then the shaft is a cylinder (approx. by n polygons) with radius r*sfac.
	If QNeedleStyle->True, then the shaft is a cone with apex at the origin.
	Package: VQM`VisualizeVector`.";

QArrowScale::usage = "QArrowScale->k is an option for QVectorToArrow. It scales the length of the
	arrow representing the vector by the factor k.
	Package: VQM`VisualizeVector`.";

QNeedleStyle::usage = "QNeedleStyle is an option for QVectorToArrow. QNeedleStyle->True causes
	the vector to be drawn as a needle. The needle is a double-cone whose style is controlled
	by the options QArrowShape, QShaftColor, QHeadColor.
	Package: VQM`VisualizeVector`.";

QMinLength::usage = "QMinLength is an option for QVectorToArrow. QMinLength->0.1 (default) specifies
	that vectors of length less than 0.1 are to be drawn as points of size QLinePointSize.
	Package: VQM`VisualizeVector`.";

QLinePointSize::usage = "QLinePointSize->0.001 is an option for QVectorToArrow. It sets the thickness of
	the lines and points representing the arrow.
	Package: VQM`VisualizeVector`.";

QHeadColor::usage = "QHeadColor->GrayLevel[0] is an option for QVectorToArrow. It sets the color for
	the head of the arrow.
	Package: VQM`VisualizeVector`.";

QShaftColor::usage = "QShaftColor->GrayLevel[0] is an option for QVectorToArrow. It sets the color
	for the shaft of the arrow, in case that VectorShaft->True.
	Package: VQM`VisualizeVector`.";

QVisualizeVector::usage = "QVisualizeVector[3vector] converts a vector into an arrow graphics
	and displays it together with other graphics elements whose appearance is controlled by the options
	QDrawUnitSphere, QDrawAxes, QCoordinateCube, QCoordinateCircles.
	Package: VQM`VisualizeVector`.";

QDrawUnitSphere::usage = "QDrawUnitSphere->n is an option for QVisualizeVector.
	It draws the outline of the unit sphere by plotting n (default: 20) circles.
	Package: VQM`VisualizeVector`.";

QDrawAxes::usage = "QDrawAxes->True, QDrawAxes->{True, False, True}. Option for QVisualizeVector etc.
	Adds gray coordinate axes to the visualization.
	Package: VQM`VisualizeVector`.";

QCoordinateCube::usage = "QCoordinateCube->False. Option for QVisualizeVector etc. QCoordinateCube->True adds a green
	cube indicating the components of the vector.
	Package: VQM`VisualizeVector`.";

QCoordinateCircles::usage = "QCoordinateCircles->True. Option for QVisualizeVector etc. Adds blue
	coordinate lines indicating the polar coordinates of the vector. QCoordinateCirclesColor->color draws
	the coordinate lines in the given color.
	Package: VQM`VisualizeVector`.";

QCoordinateCirclesColor::usage = "QCoordinateCircles->RGBColor[0,0,1]. Option for QVisualizeVector etc.
	Draws the coordinate circles in the given color (provided the option QCoordinateCircles is
	set to True.
	Package: VQM`VisualizeVector`.";
	
polyDisk::usage = "polyDisk[r] is a regular octagon with radius r in the xy-plane.
	polyDisk[r,n] is a regular polygon with n sides.
	Package: VQM`VisualizeVector`.";

polyCone::usage = "polyCone[r, h] is a set of 8 triangles approximating the shape of a half-cone
	with apex at the origin and pointing in the z-direction. polyCone[r,h,n] approximates
	the half-cone by n triangles.
	Package: VQM`VisualizeVector`.";
	
polyCylinder::usage = "polyCylinder[r1, r2, h, n(optional)] represents a cylindrical shape,
	symmetric around the z-axis. r1 is the radius at the bottom, r2 the radius at the top.
	h is the height. The cylinder is approximated by n polygons (default: 8).
	Package: VQM`VisualizeVector`.";

arrowCone::usage = "arrowCone[{pt1,pt2}, c, opts] represents a cone with color c, translated
	and rotated, so that it has apex at pt2 and points in the direction from pt1 to pt2. Its
	appearance is controlled by the option QArrowShape.
	Package: VQM`VisualizeVector`.";

arrowWithShaft::usage = "arrowWithShaft[{pt1,pt2}, colorhead, colorshaft, opts] represents
	an arrow from pt1 to pt2 (a collection of polygons). Its appearance is controlled by
	the option QArrowShape.
	Package: VQM`VisualizeVector`.";

doubleHead::usage = "doubleHead[{pt1,pt2}, color1, color2, opts] is actually a double-cone
	connecting the points pt1 with pt2. The part closer to pt1 has color1. The shape is
	controlled by the option QArrowShape.
	Package: VQM`VisualizeVector`.";
	
unitSphere::usage = "unitSphere[theta] represents a sphere with radius 1 and center at the origin
	as a collection of circles parallel to the xy-plane. The polar angles of the circles
	are multiples of theta.
	Package: VQM`VisualizeVector`.";

QRoundSphere::usage = "QRoundSphere[r,n,m,opts] is a sphere with radus r, represented by
	as a wireframe. QRoundSphere[r,n,m] is similar to Sphere[r,n,m] in the
	standard package Graphics`Shapes`, but QRoundSphere uses circles (generated by ParametricPlot3D)
	instead of polygons. The options are passed to ParametricPlot3D (useful for PlotPoints).
	Package: VQM`VisualizeVector`.";

xaxisLine::usage = "A gray line segment from (-1,0,0) to (1,0,0).
	Package: VQM`VisualizeVector`.";

yaxisLine::usage = "A gray line segment from (0,-1,0) to (0,1,0).
	Package: VQM`VisualizeVector`.";

zaxisLine::usage = "A gray line segment from (0,0,-1) to (0,0,1).
	Package: VQM`VisualizeVector`.";

axesLines::usage = "{xaxisLine,yaxisLine,zaxisLine}.
	Package: VQM`VisualizeVector`.";

coordinateCube::usage = "coordinateCube[pt] represents a rectangular shape with edges parallel to the
	coordinate axes. It has one point at the origin and one at the point pt.
	Package: VQM`VisualizeVector`.";
	
coordinateCircles::usage = "coordinateCircles[pt] are two circular arcs representing the polar
	angle and the azimuthal angle of a point pt.
	Package: VQM`VisualizeVector`.";

graphicElements::usage = "graphicElements[pt, opts] defines the graphic elements for the visualization
	of the position pt according to the values of the options QDrawUnitSphere, QDrawAxes,
	QCoordinateCube, QCoordinateCircles, QCoordinateCirclesColor.
	Package: VQM`VisualizeVector`.";


(*-----------------------------------*)
Begin["`Private`"];
(*-----------------------------------*)

(* RM: Version 6 changes *)
(* according to Compatibility/Tutorials/Geometry/Rotations *) 
Rotate3D[xyz_List,phi_,theta_,psi_]:=
(RotationMatrix[Pi-psi,{0,0,1}].RotationMatrix[theta,{1,0,0}].RotationMatrix[Pi-phi,{0,0,1}]).xyz;

(* A replacement for Sphere in Graphics`Shapes`: *)

QRoundSphere[r_, n_?IntegerQ, m_?IntegerQ, opts___] :=
	Module[{sp1,sp2, t, a}, 
		sp1 =
		ParametricPlot3D[
			Evaluate[Table[
			Abs[r]*{Sin[a]*Cos[t], Sin[a]*Sin[t], Cos[a]},
			{a, Pi/m, (m-1)*(Pi/m), Pi/m}]],
		{t, 0, 2*Pi}, opts][[1]];
		sp2 =
		ParametricPlot3D[
			Evaluate[Table[
			r*{Sin[t]*Cos[a], Sin[t]*Sin[a], Cos[t]},
			{a, 0, 2*Pi-2*Pi/n, 2*Pi/n}]],
		{t,Pi/m,(m-1)*Pi/m}, opts][[1]];
		{sp1,sp2}
	]



Options[QVectorToArrow] = {QArrowHead->True, QArrowShaft->False, QLinePointSize->0.001, QMinLength->0.1,
QArrowShape->Automatic, QHeadColor->(*Black*)White, QArrowScale->1,
QShaftColor->Black, QNeedleStyle->False };

(* some graphics elements, mainly for internal use by this package *)

(* a half-cone with apex at the origin, pointing in the z-direction *)

polyCone[r_, h_, no_:8] :=
  Table[Polygon[{{0,0,0}, {r Cos[2 k Pi/ no], r Sin[2 k Pi/no],-h},
      {r Cos[2 (k+1) Pi/ no], r Sin[2 (k+1) Pi/no],-h}}],{k,0,no-1}];

(* a polygonal cylinder *)

polyCylinder[r1_,r2_,h_,no_:8] :=
	Table[Polygon[{{r1 Cos[2 k Pi/ no], r1 Sin[2 k Pi/no], 0},
	{r2 Cos[2 k Pi/ no], r2 Sin[2 k Pi/no],h},
	{r2 Cos[2 (k+1) Pi/ no], r2 Sin[2 (k+1) Pi/no],h},
	{r1 Cos[2 (k+1) Pi/ no], r1 Sin[2 (k+1) Pi/no],0}}], {k,0,no-1}];

(* a polygonal disk, serving as the top and bottom surface of the cylinder or cone *)

polyDisk[r_,no_:8] :=
	Polygon[Table[{r Cos[2 k Pi/ no], r Sin[2 k Pi/no], 0},{k,0,no-1}]]

(* A cone, translated and rotated, so that it can serve as the head of an arrow.
	The cone has apex at pt2 and points in the direction from pt1 to pt2 *)

arrowCone[{pt1_?VectorQ,pt2_?VectorQ}, color_, opts___?OptionQ] :=
	Module[{ len = Sqrt[(pt2-pt1).(pt2-pt1)], coords=CoordinatesFromCartesian[pt2-pt1, Spherical],
			no, lenfac, headlength, headwidth, radfac, test},
			test = QArrowShape/.{opts}/.Options[QVectorToArrow];
			Which[
				test===Automatic || test===True,
					{no,lenfac,radfac,sfac}={6,1/4,1/(2 GoldenRatio),1/2},
				IntegerQ[test],
					no = test;{lenfac,radfac}={1/4,1/(2 GoldenRatio)},
				Length[test]==1,
					{no}=test;{lenfac,radfac}={1/4,1/(2 GoldenRatio)},
				Length[test]==2,
					{no,lenfac}=test;radfac=1/(2 GoldenRatio),
				Length[test]==3,
					{no,lenfac,radfac}=test,
				Length[test]==4,
					{no,lenfac,radfac,sfac}=test,
				TRUE,
					{no,lenfac,radfac,sfac}={6,1/4,1/(2 GoldenRatio),1/2}
			];
			If[lenfac < 0, headlength = -lenfac, headlength = len*lenfac, headlength = 0];
			If[radfac < 0, headwidth = -radfac, headwidth = headlength*radfac, headwidth = 0];
			{ color, Map[ (Rotate3D[#,0,coords[[2]],-coords[[3]]+Pi/2] + pt2)&,
					polyCone[headwidth, headlength, no],{3}]}
	];

(* an arrow from pt1 to pt2, with a cylinder forming the shaft, and a cone forming the head *)

arrowWithShaft[{pt1_?VectorQ,pt2_?VectorQ}, colorhead_, colorshaft_, opts___?OptionQ] :=
	Module[{ len = Sqrt[(pt2-pt1).(pt2-pt1)], coords=CoordinatesFromCartesian[pt2-pt1, Spherical],
			no, lenfac, headlength, headwidth, radfac, sfac, test, shaftwidth, shaftlength},
			test = QArrowShape/.{opts}/.Options[QVectorToArrow];
			light = Lighting ->( Lighting /. {opts}/.Options[QVectorToArrow] );
			Which[
				test===Automatic || test===True,
					{no,lenfac,radfac,sfac}={6,1/4,1/(2 GoldenRatio),1/2},
				IntegerQ[test],
					no = test;{lenfac,radfac,sfac}={1/4,1/(2 GoldenRatio),1/2},
				Length[test]==1,
					{no}=test;{lenfac,radfac,sfac}={1/4,1/(2 GoldenRatio),1/2},
				Length[test]==2,
					{no,lenfac}=test;{radfac,sfac}={1/(2 GoldenRatio),1/2},
				Length[test]==3,
					{no,lenfac,radfac}=test; sfac = 1/2,
				Length[test]==4,
					{no,lenfac,radfac,sfac}=test,
				TRUE,
					{no,lenfac,radfac,sfac}={6,1/4,1/(2 GoldenRatio),1/2}
			];
			If[lenfac < 0, headlength = -lenfac, headlength = len*lenfac, headlength = 0];
			If[radfac < 0, headwidth = -radfac, headwidth = headlength*radfac, headwidth = 0];
			If[sfac < 0, shaftwidth = -sfac, shaftwidth = sfac headwidth, shaftwidth = 0];
			If[shaftwidth > headwidth, shaftwidth = headwidth];
			shaftlength = len-headlength;
			{{colorhead, Map[ (Rotate3D[#,0,coords[[2]],-coords[[3]]+Pi/2] + pt2)&,
					polyCone[headwidth, headlength, no],{3}]},
			{colorshaft, Map[ (Rotate3D[#,0,coords[[2]],-coords[[3]]+Pi/2]+ pt1)&,
					polyCylinder[shaftwidth, shaftwidth, shaftlength, no],{3}]},
			{colorshaft, Map[ (Rotate3D[#,0,coords[[2]],-coords[[3]]+Pi/2]+ pt1)&,
					polyDisk[shaftwidth, no],{2}]},
			{colorhead, Map[ (Rotate3D[#,0,coords[[2]],-coords[[3]]+Pi/2] + pt1 + shaftlength (pt2-pt1)/len)&,
					polyDisk[headwidth, no],{2}]}}
	];

doubleHead[{pt1_?VectorQ,pt2_?VectorQ}, colorhead_, colorshaft_, opts___?OptionQ] :=
	Module[{ len = Sqrt[(pt2-pt1).(pt2-pt1)], coords=CoordinatesFromCartesian[pt2-pt1, Spherical],
			no, lenfac, headlength, headwidth, radfac, sfac, test, shaftwidth},
			test = QArrowShape/.{opts}/.Options[QVectorToArrow];
			
			Which[
				test===Automatic || test===True,
					{no,lenfac,radfac,sfac}={6,1/2,1/4,1/2},
				IntegerQ[test],
					no = test;{lenfac,radfac,sfac}={1/2,1/4,1/2},
				Length[test]==1,
					{no}=test;{lenfac,radfac,sfac}={1/2,1/4,1/2},
				Length[test]==2,
					{no,lenfac}=test;{radfac,sfac}={1/4,1/2},
				Length[test]==3,
					{no,lenfac,radfac}=test; sfac = 1/2,
				Length[test]==4,
					{no,lenfac,radfac,sfac}=test,
				TRUE,
					{no,lenfac,radfac,sfac}={6,1/2,1/4,1/2}
			];
			If[lenfac < 0, headlength = -lenfac, headlength = len*lenfac, headlength = 0];
			If[radfac < 0, headwidth = -radfac, headwidth = headlength*radfac, headwidth = 0];
			If[sfac < 0, shaftwidth = -sfac, shaftwidth = sfac headwidth, shaftwidth = 0];
			If[shaftwidth > headwidth, shaftwidth = headwidth];
			{{colorhead, Map[ (Rotate3D[#,0,coords[[2]],-coords[[3]]+Pi/2] + pt2)&,
					polyCone[headwidth, headlength, no],{3}]},
			{colorshaft, Map[ (Rotate3D[#,0,coords[[2]],-coords[[3]]+Pi/2] + pt1)&,
					polyCone[headwidth, -(len-headlength), no],{3}]}}
	];

QVectorToArrow[pt1_?VectorQ, pt2_?VectorQ, opts___?OptionQ] :=
	Module[{glow,lighting,  len1 = Sqrt[(pt2-pt1).(pt2-pt1)],
			len, arrowh, arrshft, thickn, minlen, colorhead, colorshaft,
			scale, asneedle, pt3},
		{arrowh, arrshft, thickn, minlen, colorhead, colorshaft, scale, asneedle} =
			{(*Glow, Lighting,*) QArrowHead, QArrowShaft, QLinePointSize, QMinLength, QHeadColor,
			QShaftColor, QArrowScale, QNeedleStyle}/.{opts}/.Options[QVectorToArrow];
		len = scale*len1;
		pt3 = pt1+scale*(pt2-pt1);

		Which[
			len<=minlen,Return[{(*GrayLevel[0.5],*)PointSize[thickn],Point[pt1]}],
			asneedle===True,
				Return[{ doubleHead[{pt1,pt3},  colorhead, colorshaft, opts]}],
			arrowh===True && arrshft===False,
				Return[{{Thickness[thickn],Line[{pt1,pt3}]},
					arrowCone[{pt1,pt3}, colorhead,  opts]}],
			arrowh===True && arrshft===True,
				Return[{arrowWithShaft[{pt1,pt3}, colorhead, colorshaft,  opts]}],
			arrowh===Automatic,
				Return[{{ colorhead,Thickness[thickn],Line[{pt1,pt3}]},
					{ colorhead,PointSize[1.75*thickn],Point[pt3]}}],
			True, {  colorhead,Thickness[thickn],Line[{pt1,pt3}]}
		]
	];

QVectorToArrow[pt_?VectorQ, opts___Rule] := QVectorToArrow[{0,0,0},pt, opts]/;Length[pt]==3;


(* Visualize Vector. Show vector in unit sphere with additional visual clues *)

myoptions = {QDrawUnitSphere->18, QDrawAxes->True, QCoordinateCube->False,
QCoordinateCircles->True, QCoordinateCirclesColor->Blue};

Options[QVisualizeVector] = Join[myoptions,Options[QVectorToArrow]];

QVisualizeVector[vec_?VectorQ, opts___Rule] :=
	Module[{myvec = Re[Chop[vec]]},  (* make sure that this is a real vector *)
		Show[Graphics3D[{ QVectorToArrow[myvec, opts], graphicElements[myvec,opts] }],
		FilterOptions[Graphics3D,opts], Boxed->False, Axes->False]];

(* Graphicelements for the visualization of vectors in the unit sphere: *)

(*RM: since ParemetricPlot3D adds Hue's to the Line primitives, get rid of them like this: *)
unitSphere[step_]:= Cases[
	ParametricPlot3D[
		Evaluate[
			Table[{Cos[t]Sin[the],Sin[t]Sin[the],Cos[the]},
				{the, step, Pi-step, step}]],
		{t,0,2 Pi},
		PlotPoints->200
	][[1]], _Line,-1];

xaxisLine := {GrayLevel[0.4],Thickness[0.005],
          Line[{{-1.0,0,0}, {1.0,0,0}}]};

yaxisLine := {GrayLevel[0.4],Thickness[0.005],
          Line[{{0,-1.0,0}, {0,1.0,0}}]};

zaxisLine := {GrayLevel[0.4],Thickness[0.005],
          Line[{{0,0,-1.0}, {0,0,1.0}}]};

axesLines := {xaxisLine,yaxisLine,zaxisLine};

coordinateCube[pt_] := {RGBColor[0,.6,.4],
	Line[{ {0,0,0}, {pt[[1]],0,0}, {pt[[1]],pt[[2]],0}, pt,
		{pt[[1]],0,pt[[3]]}, {0,0,pt[[3]]} }],
	Line[{ {0,0,0}, {0,pt[[2]],0}, {pt[[1]],pt[[2]],0}, pt,
		{0,pt[[2]],pt[[3]]}, {0,0,pt[[3]]} }],
	Line[{ {pt[[1]],0,0}, {pt[[1]],0,pt[[3]]} }],
	Line[{ {0,pt[[2]],0}, {0,pt[[2]],pt[[3]]} }]};

coordinateCircles[pt_,color_] :=
	Module[{coords=Simplify[CoordinatesFromCartesian[pt, Spherical]],the,phi,t},
		rad = coords[[1]];
		the = coords[[2]];
		phi = coords[[3]];
		circleattheta = Cases[
			ParametricPlot3D[rad*{Cos[t] Sin[the], Sin[t] Sin[the], Cos[the]},
				{t,0,2 Pi}, PlotPoints->100][[1]], _Line,-1];
		circleatphi = Cases[
			ParametricPlot3D[rad*{Cos[phi] Sin[t], Sin[phi] Sin[t], Cos[t]},
				{t,0,Pi}, PlotPoints->50][[1]], _Line,-1];
	{PointSize[.02],Thickness[0.003], color, Point[pt], circleattheta, circleatphi}
	];
graphicElements//Clear;
graphicElements[vec_, opts___Rule] := 
Module[{bl, ax, cc, co, ccol, graphiclist},
		{bl, ax, cc, co, ccol} = {QDrawUnitSphere, QDrawAxes,
			QCoordinateCube,QCoordinateCircles,QCoordinateCirclesColor}/.{opts}/.Options[QVisualizeVector];
		graphiclist = {};
    	If[bl===True, graphiclist =
    		Append[graphiclist, {GrayLevel[.5], unitSphere[Pi/18] }]
          ];
    	If[bl>0 && IntegerQ[bl], graphiclist =
    		Append[graphiclist, {GrayLevel[.5], unitSphere[Pi/bl]}]];
    	If[ax==True, graphiclist =
    		Append[graphiclist, axesLines]];
    	If[VectorQ[ax] && ax[[1]]==True, graphiclist =
    		Append[graphiclist, xaxisLine]];
    	If[VectorQ[ax] && ax[[2]]==True, graphiclist =
    		Append[graphiclist, yaxisLine]];
    	If[VectorQ[ax] && ax[[3]]==True, graphiclist =
    		Append[graphiclist, zaxisLine]];
    	If[cc==True, graphiclist =
    		Append[graphiclist, coordinateCube[vec]]];
    	If[co==True, graphiclist =
    		Append[graphiclist, coordinateCircles[vec,ccol]]];
    	Return[graphiclist]
    ];

(*-----------------------------------*)
End[];      (* end `Private` Context *)
(*-----------------------------------*)

Protect@@VQM`VisualizeVector`Private`Symbols;

(*-----------------------------------*)
EndPackage[]; (* end package Context *)
(*-----------------------------------*)

If[VQMmsgon, On[General::"spell"]];
If[VQMmsgon1, On[General::"spell1"]];

(* Examples:
	Show[Graphics3D[
		QVectorToArrow[{1,-.5,.7}, QArrowShaft->True,
			QHeadColor->RGBColor[1,0,0],
			QShaftColor->RGBColor[0,1,0]]
			], Lighting->None]

	Show[Graphics3D[
	QVectorToArrow[{1, .2, 0.7},
		QNeedleStyle -> True, QArrowShape->{10,1/5,1/6,1/4},
    	QHeadColor -> RGBColor[1, 0, 0],
    	QShaftColor -> RGBColor[0, 1, 1]]], 
	Lighting -> None, PlotRange -> All]
	  
	Graphics3D[polyCylinder[1,.5,1,13]]
	
	QVisualizeVector[{.6,-.2,.5},
		QCoordinateCube->True,
  		QArrowShape->{8,1/3,1/3,1/7}, QArrowShaft->True,
		QDrawUnitSphere->12]

        SphericalPlot3D[1, {\[Theta], 0.3, Pi - 0.3}, {\[Phi], 0, 2 Pi}, 
        PlotStyle -> {Opacity[0.5], Glow[White]}, Axes -> False, Mesh -> 7]
	Graphics3D[QRoundSphere[1,7,10]]
*)
