(* ::Package:: *)

If[$Notebooks,
	CurrentValue[$FrontEndSession, ShowAutoSpellCheck] = False
];

Off[CompiledFunction::cfsa];
Off[CompiledFunction::cfse];
Off[CompiledFunction::cfn];
Off[CompiledFunction::cfex];
Off[DeclarePackage::aldec];

(* RM: added $VQMDirectory *)
(* RM: added dPrint *)


(* VQM.m Master package *)



(* First, define the context of this package directory. *)
BeginPackage["VQM`"];

$VQMDirectory::usage="$VQMDirectory is the installation directory of VQM.";

(*
$VQMDirectory = If[TrueQ[$Notebooks]
 , NotebookDirectory @ EvaluationNotebook[]
 , DirectoryName[$InputFileName]
 ];
 *)
$VQMDirectory = DirectoryName[$InputFileName]

	dPrint ~ SetAttributes ~ HoldAll;

	dPrint[a__] := 
	If[$Notebooks =!= True
	  ,	
	  Print[ Style[DateString[{"Hour",":", "Minute",":", "SecondExact"}],"Text"], ": ", a ]
	  ,
	  NotebookWrite[
	    MessagesNotebook[], 
	      Cell[BoxData @ ToBoxes @ Row[{ Style[DateString[{"Hour",":", "Minute",":", "SecondExact"}],"Text"], ": ", a}]
	      ,"Print"]]
	];

(* dPrint["$VQMDirectory = ", $VQMDirectory];*)

EndPackage[];

(* Now declare all symbols that occur in any of the subpackages so the appropriate package is loaded automatically when the symbol is used *)
DeclarePackage["VQM`ArgColorPlot`", {"QArgColorPlot", "QListArgColorPlot", "QCombinedPlot",
"QListCombinedPlot", "QSpinorPlot", "QListSpinorPlot", "QSpinorCombinedPlot",
"QListSpinorCombinedPlot", "QNiceTicks", "QSaturation", "QBrightness", "QBottomLine",
"QShiftPlot", "QHorizontalRange", "QPlotDown", "QSquared", "QCurveStyle"}];

Get["VQM`ArgColorPlot`"]

DeclarePackage["VQM`ColorMaps`", {"QComplexToColor", "QRGBValues", "QComplexToRGBValues",
"QComplexToRGBColor", "QComplexToColorMap", "QValueRange", "QLightnessRange", "QSphereRadius",
"$QComplexToColorMap", "QSaturationFromLightness", "QBrightnessFromLightness",
"QHueFromArgument", "QLightnessFromModulus", "QProcessColorMapOptions", "QVectorToColor",
"QSpinorToColor"}];

DeclarePackage["VQM`ComplexPlot`", {"QComplexPlot3D", "QListComplexPlot3D",
"QComplexSurfacePlot", "QListComplexSurfacePlot", "QComplexDensityPlot", "QSpinorDensityPlot",
"QListComplexDensityPlot", "QComplexContourPlot", "QListComplexContourPlot", "QColorArrayPlot",
"QColorDensityGraphics", "QValueChecking", "QScaledValues", "QHighlights"}];

DeclarePackage["VQM`Coulomb`", {"QPrincipalQuantumNumber", "QRadialQuantumNumber",
"QCoulombEnergy", "QCoulombTimePeriod", "QCoulombFunction", "QEffectiveCoulombPotential",
"QCoulombHamiltonian", "QRadialCoulombFunction", "QRadialPositionAmplitude",
"QCoulombSpaceDimension", "QCoulombCoupling", "$QCoulombSpaceDimension", "$QCoulombCoupling"}];

DeclarePackage["VQM`Dirac1D`", {"QRelativisticEnergy", "QRelativisticVelocity", "QSignedEnergy",
"QSignedMomentum", "aplus", "aminus", "QDiracMatrix1D", "QProjectPositiveEnergy",
"QProjectNegativeEnergy", "QMomentumSpacePropagator", "QBaseSpinorRight", "QBaseSpinorLeft",
"QBaseSpinorPos", "QBaseSpinorNeg", "QDiracPlaneWaveRight", "QDiracPlaneWaveLeft",
"QDiracPlaneWavePos", "QDiracPlaneWaveNeg", "QFWMatrix", "QInverseFWMatrix", "$c", "$m",
"QDiracEquation1D", "QDiracOperator1D", "QTridiagonalBlockSolve", "QDiracTimeStep1D",
"QInitializeDirac1D", "QMomentumSpaceSpinorPos", "QMomentumSpaceSpinorNeg",
"QMomentumSpaceSpinorRight", "QMomentumSpaceSpinorLeft", "QMomentumSpaceSpinorPosRight",
"QMomentumSpaceSpinorPosLeft", "QMomentumSpaceSpinorNegRight", "QMomentumSpaceSpinorNegLeft",
"QMomentumSpaceSpinor", "QPositionSpaceSpinor", "QPositionSpaceSpinorPos",
"QPositionSpaceSpinorNeg", "QPositionSpaceSpinorRight", "QPositionSpaceSpinorLeft",
"QPositionSpaceSpinorPosRight", "QPositionSpaceSpinorPosLeft",
"QPositionSpaceSpinorNegRight", "QPositionSpaceSpinorNegLeft", "QPositionSpaceGrid",
"QPositionSpaceInterval", "QPositionSpaceStep", "QComputeMeanPosition", "QInnerProductList",
"QSpinorFT", "QInverseSpinorFT"}];

DeclarePackage["VQM`FastFourier`", {"QFourierList", "QInverseFourierList",
"QFourierListArgColorPlot", "QInverseFourierListArgColorPlot", "QStepSize", "QFourierRange",
"QGrid", "QLeftBorder", "QRightBorder", "QSpaceStep", "QSpaceInterval", "QIndexPosition",
"QFourierGrid", "QFourierLeftBorder", "QFourierRightBorder", "QFourierStep", "QFourierInterval"}];

DeclarePackage["VQM`Free`", {"QFreeHamiltonian1D", "QFreeGaussian", "QFreeGaussian1D",
"QFreeGaussian2D", "QFreeGaussian3D", "QGaussian1D", "QGaussian2D", "QGaussian3D",
"QFreeFourierGaussian1D", "QFourierGaussian1D", "QEnergyGaussian1D", "QFreeFallGaussian1D",
"QFreeMass", "QFreeSpaceDimension", "$QFreeMass", "$QFreeSpaceDimension"}];

DeclarePackage["VQM`Oscillator`", {"QOscillatorHamiltonian", "QOscillatorEnergy",
"QOscillatorFunction", "QOscillatorFunctionT", "QOscillatorGaussian", "QOscillatorBarDiagram",
"QOscillatorFrequency", "QOscillatorMass", "$QOscillatorFrequency", "$QOscillatorMass"}];

DeclarePackage["VQM`QGraphics2D`", {"QPrepareOptions", "QExtractPart", "QGetAndDensityPlot",
"QGetAndComplexDensityPlot", "QGetSpinorAndDensityPlot", "QMakeTable", "QZeroTable", "QParameters",
"QGetSpinorAndDensityPlotTwo", "QGetAndSpinorToColorPlot", "QGetAndSpinorToColorPlotTwo"}];

DeclarePackage["VQM`QuantumKernel`", {"QuantumKernel", "QuantumLink", "QNewFunction",
"QDisposeFunction", "QGetArray", "QGetFunctionInfo", "QGetColorArray", "QGetGrayArray",
"QGetRedBlueArray", "QGetBlackWhiteArray", 
(* RM2018*)"QGetAbsList", 
"QGetAbsArray", 
(*RM2018*) "QGetList", "QInfo", "QSchroedinger1D",
"QSchroedinger2D", "QSchroedinger3D", "QPauli2D", "QPauli3D", "QDirac2D", "QDirac3D",
"QDisposeOperator", "QGetOperatorInfo", "QTimeEvolution", "QBeginMovie", "QEndMovie",
"QFunctionObject", "QGetWindowInfo", "QHideWindow", "QOperatorObject", "QShowWindow"}];

DeclarePackage["VQM`Rectangular`", {"QPlaneWaveToRight", "QPlaneWaveToLeft",
"QReflectionCoefficientJump", "QTransmissionCoefficientJump", "QReflectionCoefficientWell",
"QTransmissionCoefficientWell", "QTransitionMatrix", "QSolutionWellToRight", "QPsiEvenWell",
"QPsiOddWell", "QDetEvenWell", "QDetOddWell", "QCriticalRadiusEvenWell", "QCriticalRadiusOddWell",
"QCriticalDepthEvenWell", "QCriticalDepthOddWell"}];

DeclarePackage["VQM`Spinors`", {"QNorm", "QNormalize", "$QSpinBasis", "QxBasis", "QyBasis",
"QzBasis", "QSetSpinBasis", "QSpinBasis", "QUseBasis", "QSpinorToComponents", "QComponentsToSpinor",
"QProjectUp", "QProjectDown", "QProbabilityUp", "QProbabilityDown", "QSpinorBarDiagram",
"QPauliSigma1", "QPauliSigma2", "QPauliSigma3", "QPauliSigmaV", "QIdentity2",
"QVectorToHermitianMatrix", "QVectorToDensityMatrix", "QSpinHamiltonian", "QSpinorUp",
"QSpinorDown", "QConjSpinorUp", "QConjSpinorDown", "QSpinorHarmonicUp", "QSpinorHarmonicDown",
"QSpinorToVector", "QVectorToSpinor", "QVectorLength", "QExtractPhase",
"QHermitianMatrixToVector", "QDensityMatrixToVector", "QSpinorToArrow", "QRotationSO3",
"QRotationSU2", "QVisualizeSpinor", "QVisualizeDensityMatrix"}];

DeclarePackage["VQM`VisualizeVector`", {"QVectorToArrow", "QArrowHead", "QArrowShaft", "QArrowShape",
"QArrowScale", "QNeedleStyle", "QMinLength", "QLinePointSize", "QHeadColor", "QShaftColor",
"QVisualizeVector", "QDrawUnitSphere", "QDrawAxes", "QCoordinateCube", "QCoordinateCircles",
"QCoordinateCirclesColor", "polyDisk", "polyCone", "polyCylinder", "arrowCone", "arrowWithShaft",
"doubleHead", "unitSphere", "QRoundSphere", "xaxisLine", "yaxisLine", "zaxisLine", "axesLines",
"coordinateCube", "coordinateCircles", "graphicElements"}];

DeclarePackage["VQM`VisualizeFunction1D`", {"QShowComplexPoint", "QShowComplexPointPolar", "QPlotRe", "QPlotIm",
"QPlotReIm", "QPlotAbs", "QPlotArg", "QPlotAbsArg", "QComplexFunctionGraph", "QBoxSize", "QProjectionAt"}];

(*added 20180621 *)
DeclarePackage["VQM`VisualizeFourierSum1D`", { "QAmplitudeFactor", "QBorders", "QEpilog", "QFourierSumPlot", 
 "QFourierSumWaveNumberPlot", "QFourierSumWaveStackPlot", "QHalfPeriod", "QKPlotLabel", 
 "QKRegion", "QLineSize", "QSeparationLines", "QSPlotLabel", "QStackWavesGraph", 
 "QVisualizeSuperposition", "QWaveNumberPlot", "QWaveStackPlot", 
 "QWaveStackWaveNumberPlot", "QXPlotLabel", "QXRegion"}];
