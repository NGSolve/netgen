#ifndef FILE_BitArray
#define FILE_BitArray

/**************************************************************************/
/* File:   bitarray.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

#include <limits.h>

namespace netgen
{


/**
   data type BitArray
   
   BitArray is a compressed array of Boolean information. By Set and Clear
   the whole array or one bit can be set or reset, respectively. 
   Test returns the state of the occurring bit.
   No range checking is done.

   index ranges from 0 to size-1
*/
class BitArray
{
  INDEX size;
  unsigned char * data;

public:
  BitArray ();
  ///
  BitArray (INDEX asize);
  ///
  ~BitArray ();

  /// 
  void SetSize (INDEX asize);
  ///
  INDEX Size () const
  {
    return size;
  }

  ///
  void Set ();
  ///
  void Set (INDEX i)
  {
    data[Addr(i)] |= Mask(i);
  }
  
  void Clear ();


  void Clear (INDEX i)
  {
    data[Addr(i)] &= ~Mask(i);
  }

  bool Test (INDEX i) const
  {
    return (data[i / CHAR_BIT] & (char(1) << (i % CHAR_BIT) ) ) ? true : false;
  }

  ///
  void Invert ();
  ///
  void And (const BitArray & ba2);
  ///
  void Or (const BitArray & ba2);
private:
  ///
  inline unsigned char Mask (INDEX i) const
  {
    return char(1) << (i % CHAR_BIT);
  }
  ///
  inline INDEX Addr (INDEX i) const
  {
  return (i / CHAR_BIT);
  }

  ///
  BitArray & operator= (BitArray &);
  ///
  BitArray (const BitArray &);
};



// print bitarray
inline ostream & operator<< (ostream & s, const BitArray & a)
{
  for (int i = 1; i <= a.Size(); i++)
    {
      s << int (a.Test(i));
      if (i % 40 == 0) s << "\n";
    }
  if (a.Size() % 40 != 0) s << "\n";
  return s;
}


/*
inline
INDEX BitArray :: Size () const
  {
  return size;
  }

inline
unsigned char BitArray :: Mask (INDEX i) const
  {
  return char(1) << (i % CHAR_BIT);
  }

inline
INDEX BitArray :: Addr (INDEX i) const
  {
  return (i / CHAR_BIT);
  }
inline
void BitArray :: Set (INDEX i)
  {
  data[Addr(i)] |= Mask(i);
  }

inline
void BitArray :: Clear (INDEX i)
  {
  data[Addr(i)] &= ~Mask(i);
  }


inline
int BitArray :: Test (INDEX i) const
  {
  return (data[i / CHAR_BIT] & (char(1) << (i % CHAR_BIT) ) ) ? 1 : 0;
  }

*/

}

#endif
