/**************************************************************************/
/* File:   bitarray.cc                                                    */
/* Autho: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/* 
   data type NgBitArray
*/

#include <mystdlib.h>
#include <myadt.hpp>


namespace netgen
{
  //using namespace netgen;

  NgBitArray :: NgBitArray ()
  {
    size = 0;
    data = NULL;
  }

  NgBitArray :: NgBitArray (int asize)
  {
    size = 0;
    data = NULL;
    SetSize (asize);
  }

  NgBitArray :: ~NgBitArray ()
  {
    delete [] data;
  }

  void NgBitArray :: SetSize (int asize)
  {
    if (size == asize) return;
    delete [] data;

    size = asize;
    data = new unsigned char [Addr (size)+1];
  }

  void NgBitArray :: Set ()
  {
    if (!size) return;
    for (int i = 0; i <= Addr (size); i++)
      data[i] = UCHAR_MAX;
  }

  void NgBitArray :: Clear ()
  {
    if (!size) return;
    for (int i = 0; i <= Addr (size); i++)
      data[i] = 0;
  }



  void NgBitArray :: Invert ()
  {
    if (!size) return;
    for (int i = 0; i <= Addr (size); i++)
      data[i] ^= 255;
  }

  void NgBitArray :: And (const NgBitArray & ba2)
  {
    if (!size) return;
    for (int i = 0; i <= Addr (size); i++)
      data[i] &= ba2.data[i];
  }


  void NgBitArray :: Or (const NgBitArray & ba2)
  {
    if (!size) return;
    for (int i = 0; i <= Addr (size); i++)
      data[i] |= ba2.data[i];
  }


}
