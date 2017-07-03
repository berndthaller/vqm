// ===========================================================================
//	TypeDefinition.h
//	Project: Visual Quantum Mechanics.
//  Supervision: Bernd Thaller.
//  Programming: Manfred Liebmann 1998
// ===========================================================================

#ifndef _H_TypeDefinition
#define _H_TypeDefinition
#include <math.h>

//	Float type

	// Uncomment the following line to use single precision instead of double precision
#define __USE_FLOAT__

#ifdef __USE_FLOAT__
	typedef	float		Float;
#else
	#define MLPutFloatArray MLPutDoubleArray 
	#define MLGetFloat MLGetDouble 
	typedef	double		Float;
#endif

//	Integer types
typedef long		Int32;
typedef short		Int16;
typedef char		Int8;

//	Constants
/*#ifdef WINDOWS_MATHLINK
const double pi = 3.1415926535897932384626433832795028841971693993751;
#else
extern const double	pi;
#endif
const double pi = M_PI;*/
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795028841971693993751
#endif

#endif
