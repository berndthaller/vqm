/*
 * This file automatically produced by mprep.exe from:
 *	mathlink_windows.tm
 * mprep Revision 18 Copyright (c) Wolfram Research, Inc. 1990-2013
 */

#define MPREP_REVISION 18


#include "mathlink.h"

int MLAbort = 0;
int MLDone  = 0;
long MLSpecialCharacter = '\0';
HANDLE MLInstance = (HANDLE)0;
HWND MLIconWindow = (HWND)0;

MLINK stdlink = 0;
MLEnvironment stdenv = 0;
#if MLINTERFACE >= 3
MLYieldFunctionObject stdyielder = (MLYieldFunctionObject)0;
MLMessageHandlerObject stdhandler = (MLMessageHandlerObject)0;
#else
MLYieldFunctionObject stdyielder = 0;
MLMessageHandlerObject stdhandler = 0;
#endif /* MLINTERFACE >= 3 */

#include <windows.h>

#if defined(__GNUC__)

#	ifdef TCHAR
#		undef TCHAR
#	endif
#	define TCHAR char

#	ifdef PTCHAR
#		undef PTCHAR
#	endif
#	define PTCHAR char *

#	ifdef __TEXT
#		undef __TEXT
#	endif
#	define __TEXT(arg) arg

#	ifdef _tcsrchr
#		undef _tcsrchr
#	endif
#	define _tcsrchr strrchr

#	ifdef _tcscat
#		undef _tcscat
#	endif
#	define _tcscat strcat

#	ifdef _tcsncpy
#		undef _tcsncpy
#	endif
#	define _tcsncpy _fstrncpy
#else
#	include <tchar.h>
#endif

#include <stdlib.h>
#include <string.h>
#if (WIN32_MATHLINK || WIN64_MATHLINK || __GNUC__) && !defined(_fstrncpy)
#       define _fstrncpy strncpy
#endif

#ifndef CALLBACK
#define CALLBACK FAR PASCAL
typedef LONG LRESULT;
typedef unsigned int UINT;
typedef WORD WPARAM;
typedef DWORD LPARAM;
#endif


LRESULT CALLBACK MLEXPORT
IconProcedure( HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam);

LRESULT CALLBACK MLEXPORT
IconProcedure( HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
	switch( msg){
	case WM_CLOSE:
		MLDone = 1;
		MLAbort = 1;
		break;
	case WM_QUERYOPEN:
		return 0;
	}
	return DefWindowProc( hWnd, msg, wParam, lParam);
}


#ifdef _UNICODE
#define _APISTR(i) L ## #i
#else
#define _APISTR(i) #i
#endif

#define APISTR(i) _APISTR(i)

HWND MLInitializeIcon( HINSTANCE hInstance, int nCmdShow)
{
	TCHAR path_name[260];
	PTCHAR icon_name;

	WNDCLASS  wc;
	HMODULE hdll;

	MLInstance = hInstance;
	if( ! nCmdShow) return (HWND)0;

	hdll = GetModuleHandle( __TEXT("ml32i" APISTR(MLINTERFACE)));

	(void)GetModuleFileName( hInstance, path_name, 260);

	icon_name = _tcsrchr( path_name, '\\') + 1;
	*_tcsrchr( icon_name, '.') = '\0';


	wc.style = 0;
	wc.lpfnWndProc = IconProcedure;
	wc.cbClsExtra = 0;
	wc.cbWndExtra = 0;
	wc.hInstance = hInstance;

	if( hdll)
		wc.hIcon = LoadIcon( hdll, __TEXT("MLIcon"));
	else
		wc.hIcon = LoadIcon( NULL, IDI_APPLICATION);

	wc.hCursor = LoadCursor( NULL, IDC_ARROW);
	wc.hbrBackground = (HBRUSH)( COLOR_WINDOW + 1);
	wc.lpszMenuName =  NULL;
	wc.lpszClassName = __TEXT("mprepIcon");
	(void)RegisterClass( &wc);

	MLIconWindow = CreateWindow( __TEXT("mprepIcon"), icon_name,
			WS_OVERLAPPEDWINDOW | WS_MINIMIZE, CW_USEDEFAULT,
			CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT,
			(HWND)0, (HMENU)0, hInstance, (void FAR*)0);

	if( MLIconWindow){
		ShowWindow( MLIconWindow, SW_MINIMIZE);
		UpdateWindow( MLIconWindow);
	}
	return MLIconWindow;
}


#if __BORLANDC__
#pragma argsused
#endif

#if MLINTERFACE >= 3
MLYDEFN( int, MLDefaultYielder, ( MLINK mlp, MLYieldParameters yp))
#else
MLYDEFN( devyield_result, MLDefaultYielder, ( MLINK mlp, MLYieldParameters yp))
#endif /* MLINTERFACE >= 3 */
{
	MSG msg;

#if !__BORLANDC__
	mlp = mlp; /* suppress unused warning */
	yp = yp; /* suppress unused warning */
#endif

	if( PeekMessage( &msg, (HWND)0, 0, 0, PM_REMOVE)){
		TranslateMessage( &msg);
		DispatchMessage( &msg);
	}
	return MLDone;
}


/********************************* end header *********************************/


# line 1 "mathlink_windows.tm"
// ===========================================================================
//	mathlink_windows.tm for 32bit Windows
//	Project: Visual Quantum Mechanics.
//  Supervision: Bernd Thaller.
//  Programming: Manfred Liebmann 1998, Wolfgang Thaller 2000
//
//	Changes by Wolfgang Thaller on 15.08.98:
//		- Added QSchroedinger1D function
//	Changes by Rolf Mertig, GluonVision.com  on 20.7.07:
// ===========================================================================
//
//	Mathlink templates

# line 198 "MathStubs_windows.cpp"


# line 17 "mathlink_windows.tm"
//	TFunction

# line 204 "MathStubs_windows.cpp"


# line 99 "mathlink_windows.tm"
//	TOperator

# line 210 "MathStubs_windows.cpp"


void QNewFunction P(( void));

#if MLPROTOTYPES
static int _tr0( MLINK mlp)
#else
static int _tr0(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	QNewFunction();

	res = 1;

	return res;
} /* _tr0 */


void QDisposeFunction P(( void));

#if MLPROTOTYPES
static int _tr1( MLINK mlp)
#else
static int _tr1(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	QDisposeFunction();

	res = 1;

	return res;
} /* _tr1 */


void QGetArray P(( void));

#if MLPROTOTYPES
static int _tr2( MLINK mlp)
#else
static int _tr2(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	QGetArray();

	res = 1;

	return res;
} /* _tr2 */


void QGetFunctionInfo P(( void));

#if MLPROTOTYPES
static int _tr3( MLINK mlp)
#else
static int _tr3(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	QGetFunctionInfo();

	res = 1;

	return res;
} /* _tr3 */


void QGetColorArray P(( void));

#if MLPROTOTYPES
static int _tr4( MLINK mlp)
#else
static int _tr4(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	QGetColorArray();

	res = 1;

	return res;
} /* _tr4 */


void QGetGrayArray P(( void));

#if MLPROTOTYPES
static int _tr5( MLINK mlp)
#else
static int _tr5(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	QGetGrayArray();

	res = 1;

	return res;
} /* _tr5 */


void QGetRedBlueArray P(( void));

#if MLPROTOTYPES
static int _tr6( MLINK mlp)
#else
static int _tr6(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	QGetRedBlueArray();

	res = 1;

	return res;
} /* _tr6 */


void QGetBlackWhiteArray P(( void));

#if MLPROTOTYPES
static int _tr7( MLINK mlp)
#else
static int _tr7(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	QGetBlackWhiteArray();

	res = 1;

	return res;
} /* _tr7 */


void QGetAbsArray P(( void));

#if MLPROTOTYPES
static int _tr8( MLINK mlp)
#else
static int _tr8(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	QGetAbsArray();

	res = 1;

	return res;
} /* _tr8 */


void QInfo P(( void));

#if MLPROTOTYPES
static int _tr9( MLINK mlp)
#else
static int _tr9(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if ( ! MLNewPacket(mlp) ) goto L0;
	if( !mlp) return res; /* avoid unused parameter warning */

	QInfo();

	res = 1;

L0:	return res;
} /* _tr9 */


void QSchroedinger1D P(( void));

#if MLPROTOTYPES
static int _tr10( MLINK mlp)
#else
static int _tr10(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	QSchroedinger1D();

	res = 1;

	return res;
} /* _tr10 */


void QSchroedinger2D P(( void));

#if MLPROTOTYPES
static int _tr11( MLINK mlp)
#else
static int _tr11(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	QSchroedinger2D();

	res = 1;

	return res;
} /* _tr11 */


void QSchroedinger3D P(( void));

#if MLPROTOTYPES
static int _tr12( MLINK mlp)
#else
static int _tr12(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	QSchroedinger3D();

	res = 1;

	return res;
} /* _tr12 */


void QPauli2D P(( void));

#if MLPROTOTYPES
static int _tr13( MLINK mlp)
#else
static int _tr13(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	QPauli2D();

	res = 1;

	return res;
} /* _tr13 */


void QPauli3D P(( void));

#if MLPROTOTYPES
static int _tr14( MLINK mlp)
#else
static int _tr14(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	QPauli3D();

	res = 1;

	return res;
} /* _tr14 */


void QDirac2D P(( void));

#if MLPROTOTYPES
static int _tr15( MLINK mlp)
#else
static int _tr15(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	QDirac2D();

	res = 1;

	return res;
} /* _tr15 */


void QDirac3D P(( void));

#if MLPROTOTYPES
static int _tr16( MLINK mlp)
#else
static int _tr16(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	QDirac3D();

	res = 1;

	return res;
} /* _tr16 */


void QDisposeOperator P(( void));

#if MLPROTOTYPES
static int _tr17( MLINK mlp)
#else
static int _tr17(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	QDisposeOperator();

	res = 1;

	return res;
} /* _tr17 */


void QGetOperatorInfo P(( void));

#if MLPROTOTYPES
static int _tr18( MLINK mlp)
#else
static int _tr18(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	QGetOperatorInfo();

	res = 1;

	return res;
} /* _tr18 */


void QTimeEvolution P(( void));

#if MLPROTOTYPES
static int _tr19( MLINK mlp)
#else
static int _tr19(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	QTimeEvolution();

	res = 1;

	return res;
} /* _tr19 */


static struct func {
	int   f_nargs;
	int   manual;
	int   (*f_func)P((MLINK));
	const char  *f_name;
	} _tramps[20] = {
		{ 0, 2, _tr0, "QNewFunction" },
		{ 0, 2, _tr1, "QDisposeFunction" },
		{ 0, 2, _tr2, "QGetArray" },
		{ 0, 2, _tr3, "QGetFunctionInfo" },
		{ 0, 2, _tr4, "QGetColorArray" },
		{ 0, 2, _tr5, "QGetGrayArray" },
		{ 0, 2, _tr6, "QGetRedBlueArray" },
		{ 0, 2, _tr7, "QGetBlackWhiteArray" },
		{ 0, 2, _tr8, "QGetAbsArray" },
		{ 0, 0, _tr9, "QInfo" },
		{ 0, 2, _tr10, "QSchroedinger1D" },
		{ 0, 2, _tr11, "QSchroedinger2D" },
		{ 0, 2, _tr12, "QSchroedinger3D" },
		{ 0, 2, _tr13, "QPauli2D" },
		{ 0, 2, _tr14, "QPauli3D" },
		{ 0, 2, _tr15, "QDirac2D" },
		{ 0, 2, _tr16, "QDirac3D" },
		{ 0, 2, _tr17, "QDisposeOperator" },
		{ 0, 2, _tr18, "QGetOperatorInfo" },
		{ 0, 2, _tr19, "QTimeEvolution" }
		};

static const char* evalstrs[] = {
	"Print[\"QuantumKernel 1.3 Windows, Copyright 1996-98 Manfred Lieb",
	"mann, Copyright 1998-2018 Wolfgang Thaller, Updates 2007 und 201",
	"8 by Rolf Mertig, GluonVision.com\"];",
	(const char*)0,
	"QuantumKernel::err = \"`1`.\";",
	(const char*)0,
	(const char*)0
};
#define CARDOF_EVALSTRS 2

static int _definepattern P(( MLINK, char*, char*, int));

static int _doevalstr P(( MLINK, int));

int  _MLDoCallPacket P(( MLINK, struct func[], int));


#if MLPROTOTYPES
int MLInstall( MLINK mlp)
#else
int MLInstall(mlp) MLINK mlp;
#endif
{
	int _res;
	_res = MLConnect(mlp);
	if (_res) _res = _doevalstr( mlp, 0);
	if (_res) _res = _doevalstr( mlp, 1);
	if (_res) _res = _definepattern(mlp, (char *)"QNewFunction[ arrays__ ]", (char *)"{ { arrays } }", 0);
	if (_res) _res = _definepattern(mlp, (char *)"QDisposeFunction[ function_ ]", (char *)"{ function }", 1);
	if (_res) _res = _definepattern(mlp, (char *)"QGetArray[ function_ ]", (char *)"{ function }", 2);
	if (_res) _res = _definepattern(mlp, (char *)"QGetFunctionInfo[ function_ ]", (char *)"{ function }", 3);
	if (_res) _res = _definepattern(mlp, (char *)"QGetColorArray[ function_ ]", (char *)"{ function }", 4);
	if (_res) _res = _definepattern(mlp, (char *)"QGetGrayArray[ function_ ]", (char *)"{ function }", 5);
	if (_res) _res = _definepattern(mlp, (char *)"QGetRedBlueArray[ function_ ]", (char *)"{ function }", 6);
	if (_res) _res = _definepattern(mlp, (char *)"QGetBlackWhiteArray[ function_ ]", (char *)"{ function }", 7);
	if (_res) _res = _definepattern(mlp, (char *)"QGetAbsArray[ function_ ]", (char *)"{ function }", 8);
	if (_res) _res = _definepattern(mlp, (char *)"QInfo[ ]", (char *)"{ }", 9);
	if (_res) _res = _definepattern(mlp, (char *)"QSchroedinger1D[scalar_:None,mass_,dx_]", (char *)"{ scalar, mass, dx  }", 10);
	if (_res) _res = _definepattern(mlp, (char *)"QSchroedinger2D[scalar_:None,vector_:None,domain_:None,mass_:1.,charge_:1.,units_:1.]", (char *)"{ scalar, vector, domain, mass, charge, units }", 11);
	if (_res) _res = _definepattern(mlp, (char *)"QSchroedinger3D[scalar_:None,vector_:None,domain_:None,mass_:1.,charge_:1.,units_:1.]", (char *)"{ scalar, vector, domain, mass, charge, units }", 12);
	if (_res) _res = _definepattern(mlp, (char *)"QPauli2D[scalar_:None,vector_:None,domain_:None,mass_:1.,charge_:1.,units_:1.]", (char *)"{ scalar, vector, domain, mass, charge, units }", 13);
	if (_res) _res = _definepattern(mlp, (char *)"QPauli3D[scalar_:None,vector_:None,domain_:None,mass_:1.,charge_:1.,units_:1.]", (char *)"{ scalar, vector, domain, mass, charge, units }", 14);
	if (_res) _res = _definepattern(mlp, (char *)"QDirac2D[scalar_:None,vector_:None,domain_:None,mass_:1.,charge_:1.,units_:1.]", (char *)"{ scalar, vector, domain, mass, charge, units }", 15);
	if (_res) _res = _definepattern(mlp, (char *)"QDirac3D[scalar_:None,vector_:None,domain_:None,mass_:1.,charge_:1.,units_:1.]", (char *)"{ scalar, vector, domain, mass, charge, units }", 16);
	if (_res) _res = _definepattern(mlp, (char *)"QDisposeOperator[ operator_ ]", (char *)"{ operator }", 17);
	if (_res) _res = _definepattern(mlp, (char *)"QGetOperatorInfo[ operator_ ]", (char *)"{ operator }", 18);
	if (_res) _res = _definepattern(mlp, (char *)"QTimeEvolution[ operator_, function_, timestep_, fractal_:4, steps_:1 ]", (char *)"{ operator, function, timestep, fractal, steps }", 19);
	if (_res) _res = MLPutSymbol( mlp, "End");
	if (_res) _res = MLFlush( mlp);
	return _res;
} /* MLInstall */


#if MLPROTOTYPES
int MLDoCallPacket( MLINK mlp)
#else
int MLDoCallPacket( mlp) MLINK mlp;
#endif
{
	return _MLDoCallPacket( mlp, _tramps, 20);
} /* MLDoCallPacket */

/******************************* begin trailer ********************************/

#ifndef EVALSTRS_AS_BYTESTRINGS
#	define EVALSTRS_AS_BYTESTRINGS 1
#endif

#if CARDOF_EVALSTRS
static int  _doevalstr( MLINK mlp, int n)
{
	long bytesleft, charsleft, bytesnow;
#if !EVALSTRS_AS_BYTESTRINGS
	long charsnow;
#endif
	char **s, **p;
	char *t;

	s = (char **)evalstrs;
	while( n-- > 0){
		if( *s == 0) break;
		while( *s++ != 0){}
	}
	if( *s == 0) return 0;
	bytesleft = 0;
	charsleft = 0;
	p = s;
	while( *p){
		t = *p; while( *t) ++t;
		bytesnow = (long)(t - *p);
		bytesleft += bytesnow;
		charsleft += bytesnow;
#if !EVALSTRS_AS_BYTESTRINGS
		t = *p;
		charsleft -= MLCharacterOffset( &t, t + bytesnow, bytesnow);
		/* assert( t == *p + bytesnow); */
#endif
		++p;
	}


	MLPutNext( mlp, MLTKSTR);
#if EVALSTRS_AS_BYTESTRINGS
	p = s;
	while( *p){
		t = *p; while( *t) ++t;
		bytesnow = (long)(t - *p);
		bytesleft -= bytesnow;
		MLPut8BitCharacters( mlp, bytesleft, (unsigned char*)*p, bytesnow);
		++p;
	}
#else
	MLPut7BitCount( mlp, (long_st)charsleft, (long_st)bytesleft);

	p = s;
	while( *p){
		t = *p; while( *t) ++t;
		bytesnow = t - *p;
		bytesleft -= bytesnow;
		t = *p;
		charsnow = bytesnow - MLCharacterOffset( &t, t + bytesnow, bytesnow);
		/* assert( t == *p + bytesnow); */
		charsleft -= charsnow;
		MLPut7BitCharacters(  mlp, charsleft, *p, bytesnow, charsnow);
		++p;
	}
#endif
	return MLError( mlp) == MLEOK;
}
#endif /* CARDOF_EVALSTRS */


static int  _definepattern( MLINK mlp, char *patt, char *args, int func_n)
{
	MLPutFunction( mlp, "DefineExternal", (long)3);
	  MLPutString( mlp, patt);
	  MLPutString( mlp, args);
	  MLPutInteger( mlp, func_n);
	return !MLError(mlp);
} /* _definepattern */


int _MLDoCallPacket( MLINK mlp, struct func functable[], int nfuncs)
{
#if MLINTERFACE >= 4
	int len;
#else
	long len;
#endif
	int n, res = 0;
	struct func* funcp;

	if( ! MLGetInteger( mlp, &n) ||  n < 0 ||  n >= nfuncs) goto L0;
	funcp = &functable[n];

	if( funcp->f_nargs >= 0
#if MLINTERFACE >= 4
	&& ( ! MLTestHead(mlp, "List", &len)
#else
	&& ( ! MLCheckFunction(mlp, "List", &len)
#endif
	     || ( !funcp->manual && (len != funcp->f_nargs))
	     || (  funcp->manual && (len <  funcp->f_nargs))
	   )
	) goto L0;

	stdlink = mlp;
	res = (*funcp->f_func)( mlp);

L0:	if( res == 0)
		res = MLClearError( mlp) && MLPutSymbol( mlp, "$Failed");
	return res && MLEndPacket( mlp) && MLNewPacket( mlp);
} /* _MLDoCallPacket */


mlapi_packet MLAnswer( MLINK mlp)
{
	mlapi_packet pkt = 0;
#if MLINTERFACE >= 4
	int waitResult;

	while( ! MLDone && ! MLError(mlp)
		&& (waitResult = MLWaitForLinkActivity(mlp),waitResult) &&
		waitResult == MLWAITSUCCESS && (pkt = MLNextPacket(mlp), pkt) &&
		pkt == CALLPKT)
	{
		MLAbort = 0;
		if(! MLDoCallPacket(mlp))
			pkt = 0;
	}
#else
	while( !MLDone && !MLError(mlp) && (pkt = MLNextPacket(mlp), pkt) && pkt == CALLPKT){
		MLAbort = 0;
		if( !MLDoCallPacket(mlp)) pkt = 0;
	}
#endif
	MLAbort = 0;
	return pkt;
}



/*
	Module[ { me = $ParentLink},
		$ParentLink = contents of RESUMEPKT;
		Message[ MessageName[$ParentLink, "notfe"], me];
		me]
*/

static int refuse_to_be_a_frontend( MLINK mlp)
{
	int pkt;

	MLPutFunction( mlp, "EvaluatePacket", 1);
	  MLPutFunction( mlp, "Module", 2);
	    MLPutFunction( mlp, "List", 1);
		  MLPutFunction( mlp, "Set", 2);
		    MLPutSymbol( mlp, "me");
	        MLPutSymbol( mlp, "$ParentLink");
	  MLPutFunction( mlp, "CompoundExpression", 3);
	    MLPutFunction( mlp, "Set", 2);
	      MLPutSymbol( mlp, "$ParentLink");
	      MLTransferExpression( mlp, mlp);
	    MLPutFunction( mlp, "Message", 2);
	      MLPutFunction( mlp, "MessageName", 2);
	        MLPutSymbol( mlp, "$ParentLink");
	        MLPutString( mlp, "notfe");
	      MLPutSymbol( mlp, "me");
	    MLPutSymbol( mlp, "me");
	MLEndPacket( mlp);

	while( (pkt = MLNextPacket( mlp), pkt) && pkt != SUSPENDPKT)
		MLNewPacket( mlp);
	MLNewPacket( mlp);
	return MLError( mlp) == MLEOK;
}


#if MLINTERFACE >= 3
int MLEvaluate( MLINK mlp, char *s)
#else
int MLEvaluate( MLINK mlp, charp_ct s)
#endif /* MLINTERFACE >= 3 */
{
	if( MLAbort) return 0;
	return MLPutFunction( mlp, "EvaluatePacket", 1L)
		&& MLPutFunction( mlp, "ToExpression", 1L)
		&& MLPutString( mlp, s)
		&& MLEndPacket( mlp);
}


#if MLINTERFACE >= 3
int MLEvaluateString( MLINK mlp, char *s)
#else
int MLEvaluateString( MLINK mlp, charp_ct s)
#endif /* MLINTERFACE >= 3 */
{
	int pkt;
	if( MLAbort) return 0;
	if( MLEvaluate( mlp, s)){
		while( (pkt = MLAnswer( mlp), pkt) && pkt != RETURNPKT)
			MLNewPacket( mlp);
		MLNewPacket( mlp);
	}
	return MLError( mlp) == MLEOK;
} /* MLEvaluateString */


#if __BORLANDC__
#pragma argsused
#endif

#if MLINTERFACE >= 3
MLMDEFN( void, MLDefaultHandler, ( MLINK mlp, int message, int n))
#else
MLMDEFN( void, MLDefaultHandler, ( MLINK mlp, unsigned long message, unsigned long n))
#endif /* MLINTERFACE >= 3 */
{
#if !__BORLANDC__
	mlp = (MLINK)0; /* suppress unused warning */
	n = 0;          /* suppress unused warning */
#endif

	switch (message){
	case MLTerminateMessage:
		MLDone = 1;
	case MLInterruptMessage:
	case MLAbortMessage:
		MLAbort = 1;
	default:
		return;
	}
}



#if MLINTERFACE >= 3
static int _MLMain( char **argv, char **argv_end, char *commandline)
#else
static int _MLMain( charpp_ct argv, charpp_ct argv_end, charp_ct commandline)
#endif /* MLINTERFACE >= 3 */
{
	MLINK mlp;
#if MLINTERFACE >= 3
	int err;
#else
	long err;
#endif /* MLINTERFACE >= 3 */

#if MLINTERFACE >= 4
	if( !stdenv)
		stdenv = MLInitialize( (MLEnvironmentParameter)0);
#else
	if( !stdenv)
		stdenv = MLInitialize( (MLParametersPointer)0);
#endif

	if( stdenv == (MLEnvironment)0) goto R0;

	if( !stdyielder)
#if MLINTERFACE >= 3
		stdyielder = (MLYieldFunctionObject)MLDefaultYielder;
#else
		stdyielder = MLCreateYieldFunction( stdenv,
			NewMLYielderProc( MLDefaultYielder), 0);
#endif /* MLINTERFACE >= 3 */


#if MLINTERFACE >= 3
	if( !stdhandler)
		stdhandler = (MLMessageHandlerObject)MLDefaultHandler;
#else
	if( !stdhandler)
		stdhandler = MLCreateMessageHandler( stdenv,
			NewMLHandlerProc( MLDefaultHandler), 0);
#endif /* MLINTERFACE >= 3 */


	mlp = commandline
		? MLOpenString( stdenv, commandline, &err)
#if MLINTERFACE >= 3
		: MLOpenArgcArgv( stdenv, (int)(argv_end - argv), argv, &err);
#else
		: MLOpenArgv( stdenv, argv, argv_end, &err);
#endif
	if( mlp == (MLINK)0){
		MLAlert( stdenv, MLErrorString( stdenv, err));
		goto R1;
	}

	if( MLIconWindow){
#define TEXTBUFLEN 64
		TCHAR textbuf[TEXTBUFLEN];
		PTCHAR tmlname;
		const char *mlname;
		size_t namelen, i;
		int len;
		len = GetWindowText(MLIconWindow, textbuf, 62 );
		mlname = MLName(mlp);
		namelen = strlen(mlname);
		tmlname = (PTCHAR)malloc((namelen + 1)*sizeof(TCHAR));
		if(tmlname == NULL) goto R2;

		for(i = 0; i < namelen; i++){
			tmlname[i] = mlname[i];
		}
		tmlname[namelen] = '\0';
		
#if defined(_MSC_VER) && (_MSC_VER >= 1400)
		_tcscat_s( textbuf + len, TEXTBUFLEN - len, __TEXT("("));
		_tcsncpy_s(textbuf + len + 1, TEXTBUFLEN - len - 1, tmlname, TEXTBUFLEN - len - 3);
		textbuf[TEXTBUFLEN - 2] = '\0';
		_tcscat_s(textbuf, TEXTBUFLEN, __TEXT(")"));
#else
		_tcscat( textbuf + len, __TEXT("("));
		_tcsncpy( textbuf + len + 1, tmlname, TEXTBUFLEN - len - 3);
		textbuf[TEXTBUFLEN - 2] = '\0';
		_tcscat( textbuf, __TEXT(")"));
#endif
		textbuf[len + namelen + 2] = '\0';
		free(tmlname);
		SetWindowText( MLIconWindow, textbuf);
	}

	if( MLInstance){
		if( stdyielder) MLSetYieldFunction( mlp, stdyielder);
		if( stdhandler) MLSetMessageHandler( mlp, stdhandler);
	}

	if( MLInstall( mlp))
		while( MLAnswer( mlp) == RESUMEPKT){
			if( ! refuse_to_be_a_frontend( mlp)) break;
		}

R2:	MLClose( mlp);
R1:	MLDeinitialize( stdenv);
	stdenv = (MLEnvironment)0;
R0:	return !MLDone;
} /* _MLMain */


#if MLINTERFACE >= 3
int MLMainString( char *commandline)
#else
int MLMainString( charp_ct commandline)
#endif /* MLINTERFACE >= 3 */
{
#if MLINTERFACE >= 3
	return _MLMain( (char **)0, (char **)0, commandline);
#else
	return _MLMain( (charpp_ct)0, (charpp_ct)0, commandline);
#endif /* MLINTERFACE >= 3 */
}

int MLMainArgv( char** argv, char** argv_end) /* note not FAR pointers */
{   
	static char FAR * far_argv[128];
	int count = 0;
	
	while(argv < argv_end)
		far_argv[count++] = *argv++;
		 
#if MLINTERFACE >= 3
	return _MLMain( far_argv, far_argv + count, (char *)0);
#else
	return _MLMain( far_argv, far_argv + count, (charp_ct)0);
#endif /* MLINTERFACE >= 3 */

}

#if MLINTERFACE >= 3
int MLMain( int argc, char **argv)
#else
int MLMain( int argc, charpp_ct argv)
#endif /* MLINTERFACE >= 3 */
{
#if MLINTERFACE >= 3
 	return _MLMain( argv, argv + argc, (char *)0);
#else
 	return _MLMain( argv, argv + argc, (charp_ct)0);
#endif /* MLINTERFACE >= 3 */
}
 
