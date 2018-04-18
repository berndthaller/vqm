/*
 * This file automatically produced by /usr/local/Wolfram/Mathematica/6.0/SystemFiles/Links/MathLink/DeveloperKit/Linux/CompilerAdditions/mprep from:
 *	mathlink_linux32.tm
 * mprep Revision 12 Copyright (c) Wolfram Research, Inc. 1990-2006
 */

#define MPREP_REVISION 12

#include "mathlink.h"

int MLAbort = 0;
int MLDone  = 0;
long MLSpecialCharacter = '\0';

MLINK stdlink = 0;
MLEnvironment stdenv = 0;
#if MLINTERFACE >= 3
MLYieldFunctionObject stdyielder = (MLYieldFunctionObject)0;
MLMessageHandlerObject stdhandler = (MLMessageHandlerObject)0;
#else
MLYieldFunctionObject stdyielder = 0;
MLMessageHandlerObject stdhandler = 0;
#endif /* MLINTERFACE >= 3 */

/********************************* end header *********************************/


# line 1 "mathlink_linux32.tm"
// ===========================================================================
//	mathlink_linux32.tm for 32bit Linux
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

# line 43 "MathStubs_linux32.cpp"


# line 17 "mathlink_linux32.tm"
//	TFunction

# line 49 "MathStubs_linux32.cpp"


# line 99 "mathlink_linux32.tm"
//	TOperator

# line 55 "MathStubs_linux32.cpp"


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
	char  *f_name;
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

static char* evalstrs[] = {
	"Print[\"QuantumKernel 1.2 Linux, Copyright 1996-98 Manfred Liebma",
	"nn, Copyright 1998-2000 Wolfgang Thaller, Updates 2007 by GluonV",
	"ision.com\"];",
	(char*)0,
	"QuantumKernel::err = \"`1`.\";",
	(char*)0,
	(char*)0
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
	if (_res) _res = _definepattern(mlp, "QNewFunction[ arrays__ ]", "{ { arrays } }", 0);
	if (_res) _res = _definepattern(mlp, "QDisposeFunction[ function_ ]", "{ function }", 1);
	if (_res) _res = _definepattern(mlp, "QGetArray[ function_ ]", "{ function }", 2);
	if (_res) _res = _definepattern(mlp, "QGetFunctionInfo[ function_ ]", "{ function }", 3);
	if (_res) _res = _definepattern(mlp, "QGetColorArray[ function_ ]", "{ function }", 4);
	if (_res) _res = _definepattern(mlp, "QGetGrayArray[ function_ ]", "{ function }", 5);
	if (_res) _res = _definepattern(mlp, "QGetRedBlueArray[ function_ ]", "{ function }", 6);
	if (_res) _res = _definepattern(mlp, "QGetBlackWhiteArray[ function_ ]", "{ function }", 7);
	if (_res) _res = _definepattern(mlp, "QGetAbsArray[ function_ ]", "{ function }", 8);
	if (_res) _res = _definepattern(mlp, "QInfo[ ]", "{ }", 9);
	if (_res) _res = _definepattern(mlp, "QSchroedinger1D[scalar_:None,mass_,dx_]", "{ scalar, mass, dx  }", 10);
	if (_res) _res = _definepattern(mlp, "QSchroedinger2D[scalar_:None,vector_:None,domain_:None,mass_:1.,charge_:1.,units_:1.]", "{ scalar, vector, domain, mass, charge, units }", 11);
	if (_res) _res = _definepattern(mlp, "QSchroedinger3D[scalar_:None,vector_:None,domain_:None,mass_:1.,charge_:1.,units_:1.]", "{ scalar, vector, domain, mass, charge, units }", 12);
	if (_res) _res = _definepattern(mlp, "QPauli2D[scalar_:None,vector_:None,domain_:None,mass_:1.,charge_:1.,units_:1.]", "{ scalar, vector, domain, mass, charge, units }", 13);
	if (_res) _res = _definepattern(mlp, "QPauli3D[scalar_:None,vector_:None,domain_:None,mass_:1.,charge_:1.,units_:1.]", "{ scalar, vector, domain, mass, charge, units }", 14);
	if (_res) _res = _definepattern(mlp, "QDirac2D[scalar_:None,vector_:None,domain_:None,mass_:1.,charge_:1.,units_:1.]", "{ scalar, vector, domain, mass, charge, units }", 15);
	if (_res) _res = _definepattern(mlp, "QDirac3D[scalar_:None,vector_:None,domain_:None,mass_:1.,charge_:1.,units_:1.]", "{ scalar, vector, domain, mass, charge, units }", 16);
	if (_res) _res = _definepattern(mlp, "QDisposeOperator[ operator_ ]", "{ operator }", 17);
	if (_res) _res = _definepattern(mlp, "QGetOperatorInfo[ operator_ ]", "{ operator }", 18);
	if (_res) _res = _definepattern(mlp, "QTimeEvolution[ operator_, function_, timestep_, fractal_:4, steps_:1 ]", "{ operator, function, timestep, fractal, steps }", 19);
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
#if MLPROTOTYPES
static int  _doevalstr( MLINK mlp, int n)
#else
static int  _doevalstr( mlp, n)
	 MLINK mlp; int n;
#endif
{
	long bytesleft, charsleft, bytesnow;
#if !EVALSTRS_AS_BYTESTRINGS
	long charsnow;
#endif
	char **s, **p;
	char *t;

	s = evalstrs;
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
		bytesnow = t - *p;
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
		bytesnow = t - *p;
		bytesleft -= bytesnow;
		MLPut8BitCharacters( mlp, bytesleft, (unsigned char*)*p, bytesnow);
		++p;
	}
#else
	MLPut7BitCount( mlp, charsleft, bytesleft);
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


#if MLPROTOTYPES
static int  _definepattern( MLINK mlp, char *patt, char *args, int func_n)
#else
static int  _definepattern( mlp, patt, args, func_n)
	MLINK  mlp;
	char  *patt, *args;
	int    func_n;
#endif
{
	MLPutFunction( mlp, "DefineExternal", (long)3);
	  MLPutString( mlp, patt);
	  MLPutString( mlp, args);
	  MLPutInteger( mlp, func_n);
	return !MLError(mlp);
} /* _definepattern */


#if MLPROTOTYPES
int _MLDoCallPacket( MLINK mlp, struct func functable[], int nfuncs)
#else
int _MLDoCallPacket( mlp, functable, nfuncs)
	MLINK mlp;
	struct func functable[];
	int nfuncs;
#endif
{
	long len;
	int n, res = 0;
	struct func* funcp;

	if( ! MLGetInteger( mlp, &n) ||  n < 0 ||  n >= nfuncs) goto L0;
	funcp = &functable[n];

	if( funcp->f_nargs >= 0
	&& ( ! MLCheckFunction(mlp, "List", &len)
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


#if MLPROTOTYPES
mlapi_packet MLAnswer( MLINK mlp)
#else
mlapi_packet MLAnswer( mlp)
	MLINK mlp;
#endif
{
	mlapi_packet pkt = 0;

	while( !MLDone && !MLError(mlp) && (pkt = MLNextPacket(mlp), pkt) && pkt == CALLPKT){
		MLAbort = 0;
		if( !MLDoCallPacket(mlp)) pkt = 0;
	}
	MLAbort = 0;
	return pkt;
} /* MLAnswer */



/*
	Module[ { me = $ParentLink},
		$ParentLink = contents of RESUMEPKT;
		Message[ MessageName[$ParentLink, "notfe"], me];
		me]
*/

#if MLPROTOTYPES
static int refuse_to_be_a_frontend( MLINK mlp)
#else
static int refuse_to_be_a_frontend( mlp)
	MLINK mlp;
#endif
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


#if MLPROTOTYPES
#if MLINTERFACE >= 3
int MLEvaluate( MLINK mlp, char *s)
#else
int MLEvaluate( MLINK mlp, charp_ct s)
#endif /* MLINTERFACE >= 3 */
#else
int MLEvaluate( mlp, s)
	MLINK mlp;
#if MLINTERFACE >= 3
	char *s;
#else
	charp_ct s;
#endif /* MLINTERFACE >= 3 */
#endif
{
	if( MLAbort) return 0;
	return MLPutFunction( mlp, "EvaluatePacket", 1L)
		&& MLPutFunction( mlp, "ToExpression", 1L)
		&& MLPutString( mlp, s)
		&& MLEndPacket( mlp);
} /* MLEvaluate */


#if MLPROTOTYPES
#if MLINTERFACE >= 3
int MLEvaluateString( MLINK mlp, char *s)
#else
int MLEvaluateString( MLINK mlp, charp_ct s)
#endif /* MLINTERFACE >= 3 */
#else
int MLEvaluateString( mlp, s)
	MLINK mlp;
#if MLINTERFACE >= 3
	char *s;
#else
	charp_ct s;
#endif /* MLINTERFACE >= 3 */
#endif
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


#if MLINTERFACE >= 3
#if MLPROTOTYPES
void MLDefaultHandler( MLINK mlp, int message, int n)
#else
void MLDefaultHandler( mlp, message, n)
	MLINK mlp;
	int message, n;
#endif
#else
#if MLPROTOTYPES
void MLDefaultHandler( MLINK mlp, unsigned long message, unsigned long n)
#else
void MLDefaultHandler( mlp, message, n)
	MLINK mlp;
	unsigned long message, n;
#endif
#endif /* MLINTERFACE >= 3 */
{
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

#if MLPROTOTYPES
#if MLINTERFACE >= 3
static int _MLMain( char **argv, char **argv_end, char *commandline)
#else
static int _MLMain( charpp_ct argv, charpp_ct argv_end, charp_ct commandline)
#endif /* MLINTERFACE >= 3 */
#else
static int _MLMain( argv, argv_end, commandline)
#if MLINTERFACE >= 3
  char **argv, argv_end;
  char *commandline;
#else
  charpp_ct argv, argv_end;
  charp_ct commandline;
#endif /* MLINTERFACE >= 3 */
#endif
{
	MLINK mlp;
#if MLINTERFACE >= 3
	int err;
#else
	long err;
#endif /* MLINTERFACE >= 3 */

	if( !stdenv)
		stdenv = MLInitialize( (MLParametersPointer)0);
	if( stdenv == (MLEnvironment)0) goto R0;

#if MLINTERFACE >= 3
	if( !stdhandler)
		stdhandler = (MLMessageHandlerObject)MLDefaultHandler;
#else
	if( !stdhandler)
		stdhandler = MLCreateMessageHandler( stdenv, MLDefaultHandler, 0);
#endif /* MLINTERFACE >= 3 */


	mlp = commandline
		? MLOpenString( stdenv, commandline, &err)
		: MLOpenArgcArgv( stdenv, (int)(argv_end - argv), argv, &err);
	if( mlp == (MLINK)0){
		MLAlert( stdenv, MLErrorString( stdenv, err));
		goto R1;
	}

	if( stdyielder) MLSetYieldFunction( mlp, stdyielder);
	if( stdhandler) MLSetMessageHandler( mlp, stdhandler);

	if( MLInstall( mlp))
		while( MLAnswer( mlp) == RESUMEPKT){
			if( ! refuse_to_be_a_frontend( mlp)) break;
		}

	MLClose( mlp);
R1:	MLDeinitialize( stdenv);
	stdenv = (MLEnvironment)0;
R0:	return !MLDone;
} /* _MLMain */


#if MLPROTOTYPES
#if MLINTERFACE >= 3
int MLMainString( char *commandline)
#else
int MLMainString( charp_ct commandline)
#endif /* MLINTERFACE >= 3 */
#else
#if MLINTERFACE >= 3
int MLMainString( commandline)  char *commandline;
#else
int MLMainString( commandline)  charp_ct commandline;
#endif /* MLINTERFACE >= 3 */
#endif
{
	return _MLMain( (charpp_ct)0, (charpp_ct)0, commandline);
}

#if MLPROTOTYPES
int MLMainArgv( char** argv, char** argv_end) /* note not FAR pointers */
#else
int MLMainArgv( argv, argv_end)  char **argv, **argv_end;
#endif
{   
	static char FAR * far_argv[128];
	int count = 0;
	
	while(argv < argv_end)
		far_argv[count++] = *argv++;
		 
	return _MLMain( far_argv, far_argv + count, (charp_ct)0);

}

#if MLPROTOTYPES
#if MLINTERFACE >= 3
int MLMain( int argc, char **argv)
#else
int MLMain( int argc, charpp_ct argv)
#endif /* MLINTERFACE >= 3 */
#else
#if MLINTERFACE >= 3
int MLMain( argc, argv) int argc; char **argv;
#else
int MLMain( argc, argv) int argc; charpp_ct argv;
#endif /* MLINTERFACE >= 3 */
#endif
{
#if MLINTERFACE >= 3
 	return _MLMain( argv, argv + argc, (char *)0);
#else
 	return _MLMain( argv, argv + argc, (charp_ct)0);
#endif /* MLINTERFACE >= 3 */
}
 
