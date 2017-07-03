// ===========================================================================
//	TList.h
//	Project: Visual Quantum Mechanics.
//  Supervision: Bernd Thaller.
//  Programming: Manfred Liebmann 1998
// ===========================================================================
//
//	A template class for a list of items

#ifndef _H_TList
#define _H_TList
#pragma once


template <class T>
class	TList
{
	public:
				TList( void );
				~TList( void );
		bool	New( Int32 inSize );
		Int32	Insert( T* inItem );
		T*		Remove( Int32 inID );
		T*		Fetch( Int32 inID );
		Int32	GetSize( void ) { return mSize; };

	private:
		T**		mItemH;
		Int32	mSize;
		Int32	mID;

};



// ---------------------------------------------------------------------------
//		$ TList
// ---------------------------------------------------------------------------
//	Constructor

template <class T>
TList<T>::TList( void )
{
	mItemH = NULL;
	mSize = 0;
	mID = 0;
}



// ---------------------------------------------------------------------------
//		$ ~TList
// ---------------------------------------------------------------------------
//	Destructor

template <class T>
TList<T>::~TList( void )
{
	delete [] mItemH;
}



// ---------------------------------------------------------------------------
//		$ New
// ---------------------------------------------------------------------------
//	Create a new list

template <class T>
bool	TList<T>::New( Int32 inSize )
{
	Int32	i;


	mItemH = new T*[inSize];
	if( mItemH ) {
		for( i = 0; i < inSize; i++ )
			mItemH[i] = NULL;
		mSize = inSize;
		return true;
	}
	return false;
}



// ---------------------------------------------------------------------------
//		$ Insert
// ---------------------------------------------------------------------------
//	Insert the object pointer into the list

template <class T>
Int32	TList<T>::Insert( T* inItem )
{
	mID++;
	if( (mID >= 0) && (mID < mSize) ) {
		mItemH[mID] = inItem;
		return mID;
	}
	mID = mSize;
	return 0;
}



// ---------------------------------------------------------------------------
//		$ Remove
// ---------------------------------------------------------------------------
//	Remove the object pointer from the list

template <class T>
T*	TList<T>::Remove( Int32 inID )
{
	T*	item;

	if( (inID >= 0) && (inID < mSize) ) {
		item = mItemH[inID];
		mItemH[inID] = NULL;
		return item;
	}
	return NULL;
}



// ---------------------------------------------------------------------------
//		$ Fetch
// ---------------------------------------------------------------------------
//	Fetch the object pointer from the list

template <class T>
T*	TList<T>::Fetch( Int32 inID )
{
	if( (mID >= 0) && (mID < mSize) ) {
		return mItemH[inID];
	}
	return NULL;
}


#endif
