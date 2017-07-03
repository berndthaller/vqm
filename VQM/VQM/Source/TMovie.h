// ===========================================================================
//	TMovie.h
//	Project: Visual Quantum Mechanics.
//  Supervision: Bernd Thaller.
//  Programming: Manfred Liebmann 1998
// ===========================================================================

#ifndef _H_TMovie
#define _H_TMovie

#ifdef __APPLE_CC__
#include <Carbon/Carbon.h>
#include <QuickTime/QuickTime.h>
#else
#include <ImageCompression.h>
#include <Movies.h>
#endif
#include "TypeDefinition.h"


class	TMovie
{
	public:
				TMovie( GWorldPtr inOffGWorldP );
				~TMovie( void );
		Int32	GetFrame( void ) { return mFrame; };
		Int32	Begin( void );
		Int32	AddFrame( void );
		Int32	End( void );
		
		Int32	mID;

	private:
		void	NotifyUser( void );
		
		GWorldPtr				mOffGWorldP;
		SFReply 				mSFReply;
		FSSpec					mSpec;
		Int16					mResRefNum;
		Movie					mMovie;
		Track					mTrack;
		Media					mMedia;
		ImageSequence			mSequenceID;
		ImageDescriptionHandle	mImageDesc;
		Handle					mCompressedData;
		Int32					mFrame;
};

extern bool	gSwitchedIn;

#endif
