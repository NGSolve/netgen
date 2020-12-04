/**************************************************************************/
/* File:   bitarray.cpp                                                   */
/* Autho: Joachim Schoeberl                                               */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/*
   data type BitArray
*/

#include "bitarray.hpp"

namespace ngcore
{
  BitArray :: BitArray (size_t asize)
  {
    size = 0;
    data = NULL;
    SetSize (asize);
  }

  BitArray :: BitArray (size_t asize, LocalHeap & lh)
  {
    size = asize;
    data = new (lh) unsigned char [Addr (size)+1];
    owns_data = false;
  }

  BitArray :: BitArray (const BitArray & ba2)
  {
    size = 0;
    data = NULL;
    (*this) = ba2;
  }

  void BitArray :: SetSize (size_t asize)
  {
    if (size == asize) return;
    if (owns_data)
      {
        delete [] data;
        mt.Free(Addr(size)+1);
      }

    size = asize;
    data = new unsigned char [Addr (size)+1];
    mt.Alloc(Addr(size)+1);
  }

  BitArray & BitArray :: Set () throw()
  {
    if (!size) return *this;
    for (size_t i = 0; i <= Addr (size); i++)
      data[i] = UCHAR_MAX;
    return *this;
  }

  BitArray & BitArray :: Clear () throw()
  {
    if (!size) return *this;
    for (size_t i = 0; i <= Addr (size); i++)
      data[i] = 0;
    return *this;
  }

  BitArray & BitArray :: Invert ()
  {
    if (!size) return *this;
    for (size_t i = 0; i <= Addr (size); i++)
      data[i] ^= 255;
    return *this;
  }

  BitArray & BitArray :: And (const BitArray & ba2)
  {
    if (!size) return *this;
    for (size_t i = 0; i <= Addr (size); i++)
      data[i] &= ba2.data[i];
        return *this;
  }


  BitArray & BitArray :: Or (const BitArray & ba2)
  {
    if (!size) return *this;
    for (size_t i = 0; i <= Addr (size); i++)
      data[i] |= ba2.data[i];
    return *this;
  }

  bool BitArray :: operator==(const BitArray& other) const
  {
    if(size != other.Size())
      return false;
    for(auto i : Range(size/CHAR_BIT))
      if(data[i] != other.data[i])
        return false;
    for(auto i : Range(size%CHAR_BIT))
      if(Test(i + CHAR_BIT * (size/CHAR_BIT)) != other.Test(i + CHAR_BIT * (size/CHAR_BIT)))
        return false;
    return true;
  }

  BitArray & BitArray :: operator= (const BitArray & ba2)
  {
    SetSize (ba2.Size());
    if (!size)
      return *this;
    for (size_t i = 0; i <= Addr (size); i++)
      data[i] = ba2.data[i];
    return *this;
  }

  std::ostream & operator<<(std::ostream & s, const BitArray & ba)
  {
    size_t n = ba.Size();
    for (size_t i = 0; i < n; i++)
      {
	if (i % 50 == 0) s << i << ": ";
	s << int(ba[i]);
	if (i % 50 == 49) s << "\n";
      }
    s << std::flush;
    return s;
  }

  size_t BitArray :: NumSet () const
  {
    size_t cnt = 0;
    for (size_t i = 0; i < Size(); i++)
      if (Test(i)) cnt++;
    return cnt;
  }

  void BitArray :: DoArchive(Archive& archive)
  {
    if(archive.GetVersion("netgen") >= "v6.2.2007-62")
      {
        archive.NeedsVersion("netgen", "v6.2.2007-62");
        auto size = Size();
        archive & size;
        if(archive.Input())
          SetSize(size);
        if(archive.GetVersion("netgen") < "v6.2.2009-20")
          archive.Do(data, size/CHAR_BIT+1);
        else
          {
            archive.NeedsVersion("netgen", "v6.2.2009-20");
            archive.Do(data, size/CHAR_BIT);
            for(size_t i = 0; i < size%CHAR_BIT; i++)
              {
                size_t index =  CHAR_BIT * (size/CHAR_BIT) + i;
                bool b = Test(index);
                archive & b;
                b ? SetBit(index) : Clear(index);
              }
          }
      }
    else
      {
        if (archive.Output())
          {
            throw Exception("should not get here");
            archive << Size();
            for (size_t i = 0; i < Size(); i++)
              archive << (*this)[i];
          }
        else
          {
            size_t size;
            archive & size;
            SetSize (size);
            Clear();
            for (size_t i = 0; i < size; i++)
              {
                bool b;
                archive & b;
                if (b) SetBit(i);
              }
          }
      }
  }
} // namespace ngcore
