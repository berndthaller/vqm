(* :Title:	QuantumKernel *)

(* :Name:	VQM`QuantumKernel` *)

(* :Copyright: ©2002 Bernd Thaller *)

(* :Author:	Bernd Thaller,
			Institute of Mathematics,
			University of Graz,
			A-8010 Graz
			bernd.thaller@uni-graz.at
	based on "QuantumKernel.m" by Manfred Liebmann  *)

(* :Source:
	M. Suzuki, J. Math. Phys. 32 (1991), p. 400,
	Manfred Liebmann: Diploma Thesis, University of Graz, Austria, 1998,
	Bernd Thaller: Visual Quantum Mechanics, Springer-Verlag 2000, 2003
*)

(* :Summary:
This package starts the auxiliary C++ program "QuantumKernel".
QuantumKernel implements an algorithm by M. Suzuki
(J. Math. Phys. 32 (1991), p. 400) for computing numerical solutions of the Schroedinger
and Dirac equation. This algorithm uses an operator splitting method based on
a variant of the Trotter product formula. The communication between QuantumKernel
and Mathematica is done via the MathLink protocol.
QuantumKernel receives numerical data describing
the initial function, the scalar and magnetic vector potentials, etc.,
from Mathematica and returns the data describing the solution at a
certain time (after a prescribed number of time steps).
*)

(* :Date:    Jul 20, 2007 *)

(* :Package Version: 2.2 *)

(* :Mathematica Version: 11.3.0 *)

(* :Keywords:
    Schroedinger equation, Dirac equation, numerical solution, Wavefunction
*)


(*-----------------------------------*)
BeginPackage["VQM`QuantumKernel`"];
(*-----------------------------------*)

(* QuantumKernel *)

QuantumKernel::usage = "QuantumKernel.m provides a MathLink to an external C++ program.
It serves to compute numerical solutions of the Schroedinger, Pauli, and Dirac equations."

QuantumLink::usage = "QuantumLink";

(* TFunction Objects *)

QNewFunction::usage = "QNewFunction[A,..] generates a function-object (of type
'TFunction') for QuantumKernel from a list of arrays of real numbers.
The function object contains the numerical data
representing the numerical discretization of vector-valued functions.
The idea is that these data are changed by some numerical computations
performed by QuantumKernel. You can read the changed values of arrays by
QGetArray[f]. Here f is the expression returned by
QNewFunction[arrays]. This expression is something like 'QFunctionObject[number]'
and it refers to the corresponding data structure in the program 'QuantumKernel'
(an object of type 'TFunction').
The dimensions of the arrays in the argument of QNewFunction
depend on the dimension of the numerical domain.
The number of arrays depends on the dimension of the data.
Usually, real- and imaginary parts of each component are expected.
Hence a complex-valued function
is represented by QNewFunction[Re[complexarray],Im[complexarray]] and a real-valued
function is obtained by QNewFunction[realarray, nullarray] (nullarray has the same
dimensions as realarray and has all elements 0.). All elements of the arrays should
be real numbers (integer values are converted to floats).
Combining arrays with different dimensions or with non-numerical elements
produces the error message 'QuantumKernel::err: out of sequence.'.
Package: VQM`QuantumKernel`.";

QDisposeFunction::usage = "QDisposeFunction[f] deletes the numerical
data representing the function f from QuantumKernel. Here 'function'
refers to a function-object of type TFunction (a data structure of
QuantumKernel). An expression
suitable as argument for QDisposeFunction is returned by QNewFunction.
Package: VQM`QuantumKernel`.";

QGetArray::usage = "QGetArray[f] returns the numerical values contained
in the function object f. ('function object' is obtained when calling
QNewFunction). Package: VQM`QuantumKernel`.";

QGetFunctionInfo::usage = "QGetFunctionInfo[f] gives some information
about the function object refered to by f. This reference is obtained
from 'QNewFunction'. Package: VQM`QuantumKernel`.";
	
QGetColorArray::usage = "QGetColorArray[f]. This returns an array of
RGBcolor-values. These RGB values are obtained from the function-object 
f. Here 'f' refers to the data structure holding an array of complex values.
These values are transformed to colors via a colormap like the one defined in VQM`ColorMaps`. 
An expression suitable as an argument for QGetColorArray is returned by QNewFunction.
Package: VQM`QuantumKernel`.";

QGetGrayArray::usage = "QGetGrayArray[f]. Package: VQM`QuantumKernel`.";
QGetRedBlueArray::usage = "QGetRedBlueArray[f]. Package: VQM`QuantumKernel`.";
QGetBlackWhiteArray::usage = "QGetBlackWhiteArray[f]. Package: VQM`QuantumKernel`.";
QGetAbsArray::usage = "QGetAbsArray[f]. Package: VQM`QuantumKernel`.";

QInfo::usage = "QInfo[] returns informations about the state of QuantumKernel.
It lists informations about all TFunction and TOperator objects. Package: VQM`QuantumKernel`.";

(* TOperator Objects*)

QSchroedinger1D::usage = "QSchroedinger1D[V, m, dx] generates a
data structure for QuantumKernel (of type 'TOperator'). It needs a scalar
potential V (a complex function object). The reference 'V' is obtained
by executing the command QNewFunction[Re[list],Im[list]]. The real number 'm' defines
the mass of the particle used in the Schroedinger operator. 'dx' is the size
of the space grid (only uniform grids are supported). QSchroedinger1D returns
a reference to the 'operator-object'. This is needed to specify the time evolution,
see QTimeEvolution. Package: VQM`QuantumKernel`.";
	
QSchroedinger2D::usage = "QSchroedinger2D[V, A, Dom, mass, charge, dx]
generates a data structure for QuantumKernel (a 'TOperator'-object) that represents
a Schroedinger operator in two dimensions. 'V' refers to a complex function
object (complex scalar potential), 'A' is a vectorfield with two
components (vector potential), 'Dom' is a real scalar field
whose positive values describe
the domain of the simulation. 'dx' is the step-size of the spatial grid.
'mass' and 'charge' are real-valued parameters. Package: VQM`QuantumKernel`.";
	
QSchroedinger3D::usage = "QSchroedinger3D[V, A, Dom, mass, charge, dx]
generates a data structure for QuantumKernel (a 'TOperator'-object) that represents
a Schroedinger operator in three dimensions. 'V' refers to a complex function
object (complex scalar potential), 'A' is a vectorfield with three
components (vector potential), 'Dom' is a real scalar field
whose positive values describe
the domain of the simulation. 'dx' is the step-size of the spatial grid.
'mass' and 'charge' are real-valued parameters. The fields are generated with
QNewFunction. Package: VQM`QuantumKernel`.";
	
QPauli2D::usage = "QPauli2D[V, A, Dom, mass, charge, dx]
generates a data structure for QuantumKernel (a 'TOperator'-object) that represents
a Pauli operator in two dimensions. 'V' refers to a complex function
object (complex scalar potential), 'A' is a vectorfield with two
components (vector potential), 'Dom' is a real scalar field
whose positive values describe
the domain of the simulation. 'dx' is the step-size of the spatial grid.
'mass' and 'charge' are real-valued parameters. The fields are generated with
QNewFunction. Package: VQM`QuantumKernel`.";

QPauli3D::usage = "QPauli3D[V, A, Dom, mass, charge, dx]
generates a data structure for QuantumKernel (a 'TOperator'-object) that represents
a Pauli operator in three dimensions. 'V' refers to a complex function
object (complex scalar potential), 'A' is a vectorfield with three
components (vector potential), 'Dom' is a real scalar field
whose positive values describe
the domain of the simulation. 'dx' is the step-size of the spatial grid.
'mass' and 'charge' are real-valued parameters. The fields are generated with
QNewFunction. Package: VQM`QuantumKernel`.";

QDirac2D::usage = "QDirac2D[V, A, Dom, mass, charge, dx]
generates a data structure for QuantumKernel (a 'TOperator'-object) that represents
a Dirac operator in two dimensions. 'V' refers to a complex function
object (complex Lorentz-scalar potential), 'A' is a vectorfield with three
components (electromagnetic vector potential), 'Dom' is a real scalar field whose
positive values describe
the domain of the simulation. 'dx' is the step-size of the spatial grid.
'mass' and 'charge' are real-valued parameters. The fields are generated with
QNewFunction. Package: VQM`QuantumKernel`.";

QDirac3D::usage = "QDirac3D[V, A, Dom, mass, charge, dx]
generates a data structure for QuantumKernel (a 'TOperator'-object) that represents
a Dirac operator in three dimensions. 'V' refers to a complex function
object (complex Lorentz scalar potential), 'A' is a vectorfield with four
components (vector potential), 'Dom' is a real scalar field
whose positive values describe
the domain of the simulation. 'dx' is the step-size of the spatial grid.
'mass' and 'charge' are real-valued parameters. The fields are generated with
QNewFunction. Package: VQM`QuantumKernel`.";

QDisposeOperator::usage = "QDisposeOperator[operator] deletes the data describing the
TOperator object 'operator'. Package: VQM`QuantumKernel`.";
	
QGetOperatorInfo::usage = "QGetOperatorInfo[operator] gives information about the
TOperator object 'operator'. Package: VQM`QuantumKernel`.";

QTimeEvolution::usage = "QTimeEvolution[operator, function, timestep, order, steps]
computes the time evolution generated by 'operator' (which refers to a TOperator-object)
for the initial function 'function' (which refers to a TFunction-object). You
can use any of the operators defined by the QuantumKernel package with the
corresponding compatible wave function. (For example, 'function' must have
two complex components or four real components in case of QPauli3D or QDirac2D.)
The real number 'timestep' is the length of one time step. The integer 'order'
describes the order of the method (a Trotter-Suzuki decomposition of the
exponential operator, see M.Suzuki, J.Math.Phys.32 (1991), 410). The integer
'step' describes the number of time steps to be performed. Package: VQM`QuantumKernel`.";


(*-----------------------------------*)
Begin["`Private`"];
(*-----------------------------------*)

SetDirectory @ DirectoryName[$InputFileName];

	QuantumLink = 
		Install[
			Which[
(* RM: fresh compiled *)
				$SystemID === "Linux",
					"QuantumKernel32",
(* RM: fresh compiled *)
				$SystemID === "Linux-x86-64",
					"QuantumKernel64",
(* RM: changed again to QuantumKernelX*)
				$OperatingSystem == "MacOSX", 
					"QuantumKernelX",
(* RM: fresh compiled *)
 				$OperatingSystem === "Windows",
                                        "QuantumKernel.exe"
                             ]
		];

(*-----------------------------------*)
End[];      (* end `Private` Context *)
(*-----------------------------------*)

(*-----------------------------------*)
EndPackage[]; (* end package Context *)
(*-----------------------------------*)

