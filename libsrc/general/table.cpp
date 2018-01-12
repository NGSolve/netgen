/**************************************************************************/
/* File:   table.cpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/* 
   Abstract data type TABLE
*/

#include <mystdlib.h>
#include <myadt.hpp>

namespace netgen
{
  //using namespace netgen;

  BASE_TABLE :: BASE_TABLE (int size)
    : data(size)
  {
    for (int i = 0; i < size; i++)
      {
	data[i].maxsize = 0;
	data[i].size = 0;
	data[i].col = NULL;
      }
    oneblock = NULL;
  }

  BASE_TABLE :: BASE_TABLE (const FlatArray<int> & entrysizes, int elemsize)
    : data(entrysizes.Size())
  {
    size_t cnt = 0;
    size_t n = entrysizes.Size();

    for (size_t i = 0; i < n; i++)
      cnt += entrysizes[i];
    oneblock = new char[elemsize * cnt];
    // mem_total_alloc_table += elemsize * cnt;

    cnt = 0;
    for (size_t i = 0; i < n; i++)
      {
	data[i].maxsize = entrysizes[i];
	data[i].size = 0;

	data[i].col = &oneblock[elemsize * cnt];
	cnt += entrysizes[i];
      }
  }

  BASE_TABLE :: ~BASE_TABLE ()
  {
    if (oneblock)
      delete [] oneblock;
    else
      {
	for (int i = 0; i < data.Size(); i++)
	  delete [] (char*)data[i].col;
      }
  }
  
  void BASE_TABLE :: SetSize (int size)
  {
    for (int i = 0; i < data.Size(); i++)
      delete [] (char*)data[i].col;
    
    data.SetSize(size);
    for (int i = 0; i < size; i++)
      {
	data[i].maxsize = 0;
	data[i].size = 0;
	data[i].col = NULL;
      }    
  }
  
  void BASE_TABLE :: ChangeSize (int size)
  {
    int oldsize = data.Size();
    if (size == oldsize) 
      return;

    if (size < oldsize)
      for (int i = size; i < oldsize; i++)
	delete [] (char*)data[i].col;
    
    data.SetSize(size);

    for (int i = oldsize; i < size; i++)
      {
	data[i].maxsize = 0;
	data[i].size = 0;
	data[i].col = NULL;
      }    
  }

  void BASE_TABLE :: IncSize2 (int i, int elsize)
  {
#ifdef DEBUG
    if (i < 0 || i >= data.Size())
      {
	MyError ("BASE_TABLE::Inc: Out of range");
	return;
      }
#endif

    linestruct & line = data[i];
    if (line.size == line.maxsize)
      {
	void * p = new char [(line.maxsize+5) * elsize];
      
	memcpy (p, line.col, line.maxsize * elsize);
	delete [] (char*)line.col;

	line.col = p;
	line.maxsize += 5;
      }
  
    line.size++;
  }




  void BASE_TABLE :: SetEntrySize2 (int i, int newsize, int elsize)
  {
    linestruct & line = data[i];
    if (newsize > line.maxsize)
      {
	void * p = new char [newsize * elsize];
      
	memcpy (p, line.col, min2 (newsize, line.size) * elsize);
	delete [] (char*)line.col;

	line.col = p;
      }

    line.size = newsize;
  }





  /*
  void BASE_TABLE :: DecSize (int i)
  {
#ifdef DEBUG
    if (i < 0 || i >= data.Size())
      {
	MyError ("BASE_TABLE::Dec: Out of range");
	return;
      }
#endif

    linestruct & line = data[i];
  
#ifdef DEBUG
    if (line.size == 0)
      {
	MyError ("BASE_TABLE::Dec: EntrySize < 0");
	return;      
      }
#endif
  
    line.size--;
  }
  */



  void BASE_TABLE :: AllocateElementsOneBlock (int elemsize)
  {
    size_t cnt = 0;
    size_t n = data.Size();

    for (size_t i = 0; i < n; i++)
      cnt += data[i].maxsize;
    oneblock = new char[elemsize * cnt];

    cnt = 0;
    for (size_t i = 0; i < n; i++)
      {
	data[i].size = 0;
	data[i].col = &oneblock[elemsize * cnt];
	cnt += data[i].maxsize;
      }
  }



  size_t BASE_TABLE :: AllocatedElements () const
  {
    size_t els = 0;
    for (size_t i = 0; i < data.Size(); i++)
      els += data[i].maxsize;
    return els;
  }
  
  size_t BASE_TABLE :: UsedElements () const
  {
    size_t els = 0;
    for (size_t i = 0; i < data.Size(); i++)
      els += data[i].size;
    return els;
  }

  void BASE_TABLE :: SetElementSizesToMaxSizes ()
  {
    for (int i = 0; i < data.Size(); i++)
      data[i].size = data[i].maxsize;
  }

}
