#ifndef FILE_STACK
#define FILE_STACK

/*****************************************************************************/
/*  File: stack.hh                                                           */
/*  Author: Wolfram Muehlhuber                                               */
/*  Date: September 98                                                       */
/*****************************************************************************/

/*
  
  Stack class, based on a resizable array

 */


// #include "array.hpp"

namespace netgen
{

///
template <class T> class STACK
{
public:
  ///
  inline STACK (INDEX asize = 0, INDEX ainc = 0);
  ///
  inline ~STACK ();

  ///
  inline void Push (const T & el);
  ///
  inline T & Pop ();
  ///
  const inline T & Top () const;
  ///
  inline int IsEmpty () const;
  ///
  inline void MakeEmpty ();

private:
  ///
  NgArray<T> elems;
  ///
  INDEX size;
};




/*
  
  Stack class, based on a resizable array

 */

template <class T>
inline STACK<T> :: STACK (INDEX asize, INDEX ainc)
  : elems(asize, ainc)
{
  size = 0;
}


template <class T>
inline STACK<T> :: ~STACK ()
{
  ;
}


template <class T> 
inline void STACK<T> :: Push (const T & el)
{
  if (size < elems.Size())
    elems.Elem(++size) = el;
  else
    {
      elems.Append(el);
      size++;
    }
}


template <class T> 
inline T & STACK<T> :: Pop ()
{
  return elems.Elem(size--);
}


template <class T>
const inline T & STACK<T> :: Top () const
{
  return elems.Get(size);
}

template <class T>
inline int STACK<T> :: IsEmpty () const
{
  return (size == 0);
}


template <class T>
inline void STACK<T> :: MakeEmpty ()
{
  size = 0;
}

}

#endif
