(* Wolfram Language Test file *)
Needs["VQM`QuantumKernel`"]

Print["$Version = ", $Version];

(* QGetArray Test ********************** *)
(* QGetAbsArray Test ********************** *)
Test[
	Block[{startTime, numleft, numright, dx, p0, x0, x1, psi0, psi, psitemp, abspsi, testRes}, 
		startTime = AbsoluteTime[];

pos = {-3., -3.}; mom = {4., 4.}; a = 1.; 
Gauss[x0_, y0_, kx_, ky_, a_] := Compile @@ 
          {{x, y}, Simplify[(a/Pi)^(1/2)*Exp[I*(kx*x + ky*y)]*
                   Exp[-((a*((x - x0)^2 + (y - y0)^2))/2)]]}; 
f = Gauss[pos[[1]], pos[[2]], mom[[1]], mom[[2]], 1]; 
numleft = {-7., -7.}; numright = {7., 7.}; 
dx = 0.14; 
psi0 = Table[
   f[x, y], {y, numleft[[2]] + dx/2, numright[[2]] - dx/2, dx}, 
          {x, numleft[[1]] + dx/2, numright[[1]] - dx/2, dx}]; 
psi = QNewFunction[Re[psi0], Im[psi0]];
Emean = mom.mom/2;
potfnc = Compile[{x, y}, Emean/Cosh[2 (x + y)]]; 
V0 = Table[potfnc[x, y],
        {y, numleft[[2]] + dx/2, numright[[2]] - dx/2, dx},
        {x, numleft[[1]] + dx/2, numright[[1]] - dx/2, dx}];
V = QNewFunction[Re[V0], Im[V0]];
H = QSchroedinger2D[V, None, None, 1., 1., dx];
Print["H = ", H];
dt = 0.01; ordr = 4; reps = 20;
QTimeEvolution[H, psi, dt, ordr, 80];

psi1 = QGetArray[psi];
Print["psi1[[{1,2},44,44]] = " ,psi1[[{1,2},44,44]] ];
Print["Norm @ psi1[[{1,2},44,44]] = " , Norm[psi1[[{1,2},44,44]] ]];

(*Print["********** psitemp = ", MatrixQ[psitemp, MachineNumberQ], " *******************"];*)
	  abspsiar = QGetAbsArray[psi];
Print["Dimensions[abspsiar] = ", abspsiar//Dimensions];
Print["abspsiar[[44,44]] = ", abspsiar[[44,44]]];
Print["abspsiar[[50,50]] = ", abspsiar[[50,50]]];
      Print["QGetArray test finished ", testRes, " time used = ", AbsoluteTime[] - startTime];
      testRes
	]
	,
	True
	,
	TestID->"QGetArray-20180426-N8H1X8"
]