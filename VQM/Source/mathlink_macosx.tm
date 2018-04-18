// ===========================================================================
//	mathlink_macosx.tm
//	Project: Visual Quantum Mechanics.
//  Supervision: Bernd Thaller.
//  Programming: Manfred Liebmann 1998, Wolfgang Thaller 2000
//
//	Changes by Wolfgang Thaller on 15.08.98:
//		- Added QSchroedinger1D function
//	Changes by Rolf Mertig  on 28.6.07:
// ===========================================================================
//
//	Mathlink templates

:Evaluate: Print["QuantumKernel 1.2 MacOS X, Copyright \[Copyright]1996-98 Manfred Liebmann, \[Copyright]1998-2000 Wolfgang Thaller, Updates 2007 by GluonVision"];
:Evaluate: QuantumKernel::err = "`1`.";

//	TFunction

:Begin:
:Function: QNewFunction
:Pattern: QNewFunction[ arrays__ ]
:Arguments: { { arrays } }
:ArgumentTypes: { Manual }
:ReturnType: Manual
:End:

:Begin:
:Function: QDisposeFunction
:Pattern: QDisposeFunction[ function_ ]
:Arguments: { function }
:ArgumentTypes: { Manual }
:ReturnType: Manual
:End:

:Begin:
:Function: QGetArray
:Pattern: QGetArray[ function_ ]
:Arguments: { function }
:ArgumentTypes: { Manual }
:ReturnType: Manual
:End:

:Begin:
:Function: QGetFunctionInfo
:Pattern: QGetFunctionInfo[ function_ ]
:Arguments: { function }
:ArgumentTypes: { Manual }
:ReturnType: Manual
:End:

:Begin:
:Function: QGetColorArray
:Pattern: QGetColorArray[ function_ ]
:Arguments: { function }
:ArgumentTypes: { Manual }
:ReturnType: Manual
:End:

:Begin:
:Function: QGetGrayArray
:Pattern: QGetGrayArray[ function_ ]
:Arguments: { function }
:ArgumentTypes: { Manual }
:ReturnType: Manual
:End:

:Begin:
:Function: QGetRedBlueArray
:Pattern: QGetRedBlueArray[ function_ ]
:Arguments: { function }
:ArgumentTypes: { Manual }
:ReturnType: Manual
:End:

:Begin:
:Function: QGetBlackWhiteArray
:Pattern: QGetBlackWhiteArray[ function_ ]
:Arguments: { function }
:ArgumentTypes: { Manual }
:ReturnType: Manual
:End:

:Begin:
:Function: QGetAbsArray
:Pattern: QGetAbsArray[ function_ ]
:Arguments: { function }
:ArgumentTypes: { Manual }
:ReturnType: Manual
:End:

:Begin:
:Function: QInfo
:Pattern: QInfo[ ]
:Arguments: { }
:ArgumentTypes: { }
:ReturnType: Manual
:End:

//	TOperator

:Begin:
:Function: QSchroedinger1D
:Pattern: QSchroedinger1D[scalar_:None,mass_,dx_]
:Arguments: { scalar, mass, dx  }
:ArgumentTypes: { Manual }
:ReturnType: Manual
:End:

:Begin:
:Function: QSchroedinger2D
:Pattern: QSchroedinger2D[scalar_:None,vector_:None,domain_:None,mass_:1.,charge_:1.,units_:1.]
:Arguments: { scalar, vector, domain, mass, charge, units }
:ArgumentTypes: { Manual }
:ReturnType: Manual
:End:

:Begin:
:Function: QSchroedinger3D
:Pattern: QSchroedinger3D[scalar_:None,vector_:None,domain_:None,mass_:1.,charge_:1.,units_:1.]
:Arguments: { scalar, vector, domain, mass, charge, units }
:ArgumentTypes: { Manual }
:ReturnType: Manual
:End:

:Begin:
:Function: QPauli2D
:Pattern: QPauli2D[scalar_:None,vector_:None,domain_:None,mass_:1.,charge_:1.,units_:1.]
:Arguments: { scalar, vector, domain, mass, charge, units }
:ArgumentTypes: { Manual }
:ReturnType: Manual
:End:

:Begin:
:Function: QPauli3D
:Pattern: QPauli3D[scalar_:None,vector_:None,domain_:None,mass_:1.,charge_:1.,units_:1.]
:Arguments: { scalar, vector, domain, mass, charge, units }
:ArgumentTypes: { Manual }
:ReturnType: Manual
:End:

:Begin:
:Function: QDirac2D
:Pattern: QDirac2D[scalar_:None,vector_:None,domain_:None,mass_:1.,charge_:1.,units_:1.]
:Arguments: { scalar, vector, domain, mass, charge, units }
:ArgumentTypes: { Manual }
:ReturnType: Manual
:End:

:Begin:
:Function: QDirac3D
:Pattern: QDirac3D[scalar_:None,vector_:None,domain_:None,mass_:1.,charge_:1.,units_:1.]
:Arguments: { scalar, vector, domain, mass, charge, units }
:ArgumentTypes: { Manual }
:ReturnType: Manual
:End:

:Begin:
:Function: QDisposeOperator
:Pattern: QDisposeOperator[ operator_ ]
:Arguments: { operator }
:ArgumentTypes: { Manual }
:ReturnType: Manual
:End:

:Begin:
:Function: QGetOperatorInfo
:Pattern: QGetOperatorInfo[ operator_ ]
:Arguments: { operator }
:ArgumentTypes: { Manual }
:ReturnType: Manual
:End:

:Begin:
:Function: QTimeEvolution
:Pattern: QTimeEvolution[ operator_, function_, timestep_, fractal_:4, steps_:1 ]
:Arguments: { operator, function, timestep, fractal, steps }
:ArgumentTypes: { Manual }
:ReturnType: Manual
:End:

