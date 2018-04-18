(*
Evaluate this notebook if you want to use the old symbol names (without "Q")
from the packages ArgColorPlot.m, ComplexPlot.m, as defined in
Visual Quantum Mechanics ("Book One").

For Advanced Visual Quantum Mechanics ("Book Two"), the context of all
packages is VQM`.

Visual Quantum Mechanics (Springer-Verlag 2000) came with three packages
Graphics`ArgColorPlot`
Graphics`ComplexPlot`

Most of the symbols defined in these packages are now replaced by symbols
defined in the following contexts
VQM`ArgColorPlot`
VQM`ComplexPlot`
VQM`ColorMaps`
The new symbols now have a "Q" in front of their name, in order to prevent
possible conflicts that may arise in the further development of Mathematica
or other third-party packages.

As a rule: Use the packages from Book One for the notebooks from Book One,
and use the packages from Book Two with the notebooks from Book Two.

If you absolutely want to use old symbols together with the new packages,
you can get rid of the "Q" in front of those names, that have been defined
for Book One:

<<VQM`; (* load new package declarations *)
<<VQM`Compatibility`; (* define names without Q *)

Note: Book One also had the package
QuantumMechanics`QuantumKernel`
whose symbols already followed the new naming convention.
This package is now replaced by
VQM`QuantumKernel`

*)

ComplexPlot3D = QComplexPlot3D;
ListComplexPlot3D = QListComplexPlot3D;
ListComplexSurfacePlot = QListComplexSurfacePlot;
ComplexDensityPlot = QComplexDensityPlot;
ListComplexDensityPlot = QListComplexDensityPlot;
ComplexContourPlot = QComplexContourPlot;
ListComplexContourPlot = QListComplexContourPlot;
ColorArrayPlot = QColorArrayPlot;
ColorDensityGraphics = QColorDensityGraphics;
ComplexToColor = QComplexToColor;
ComplexToColorMap = QComplexToColorMap;
ValueRange = QValueRange;
LightnessRange = QLightnessRange;
SphereRadius = QSphereRadius;
Highlights = QHighlights;
ValueChecking = QValueChecking;
ScaledValues = QScaledValues;
$ComplexToColorMap = $QComplexToColorMap;
ArgColorPlot = QArgColorPlot;
ListArgColorPlot = QListArgColorPlot;
CombinedPlot = QCombinedPlot;
ListCombinedPlot = QListCombinedPlot;
NiceTicks = QNiceTicks;
Saturation = QSaturation;
Brightness = QBrightness;
