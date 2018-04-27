(* Wolfram Language Test file *)
Needs["VQM`QuantumKernel`"]

Print["$Version = ", $Version];

(* QGetArray Test ********************** *)
(* QGetAbsArray Test ********************** *)
Test[
	Block[{startTime, numleft, numright, dx, p0, x0, x1, psi0, psi, psitemp, abspsi, testRes}, 
      Print["QGetArray test started"]; startTime = AbsoluteTime[];
      numleft = -20;      (* left border of the domain *)
	  numright = 20;      (* right border of the domain *)
	  dx = 0.02;          (* space step *)
	  p0 = 8.;            (* initial momentum *)
	  x0 = -7.;           (* initial position *)
	  psi0 = Table[ Module[ {x=x1-x0}, (1./Pi)^(1/4) Exp[I*p0*x-x*x/2] ], {x1,numleft,numright,dx}]; 
	  psi = QNewFunction[Re[psi0], Im[psi0]];
	  psitemp = QGetArray[psi];
	  abspsi = QGetAbsList[psi];
Print["abspsi = ",abspsi];
	  testRes = Union[ Flatten @ Chop[#, 10^-7]& @ ( {Re[psi0], Im[psi0]} - psitemp ) ] === {0};
      Print["QGetArray test finished ", testRes, " time used = ", AbsoluteTime[] - startTime];
      testRes
	]
	,
	True
	,
	TestID->"QGetArray-20180426-N8H1X8"
]