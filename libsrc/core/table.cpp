/**************************************************************************/
/* File:   table.cpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   25. Mar. 2000                                                  */
/**************************************************************************/

/*
   Abstract data type Table
*/

#include "table.hpp"

namespace ngcore
{
  template <typename TI>
  size_t * TablePrefixSum2 (FlatArray<TI> entrysize)
  {
    size_t size  = entrysize.Size();
    size_t * index = new size_t[size+1];

    if (entrysize.Size() < 100)
      {
        size_t mysum = 0;
        for (size_t i = 0; i < entrysize.Size(); i++)
          {
            index[i] = mysum;
            mysum += entrysize[i];
          }
        index[entrysize.Size()] = mysum;
        return index;
      }

    
    Array<size_t> partial_sums(TaskManager::GetNumThreads()+1);
    partial_sums[0] = 0;
    ParallelJob
      ([&] (TaskInfo ti)
       {
         IntRange r = IntRange(size).Split(ti.task_nr, ti.ntasks);
         size_t mysum = 0;
         for (size_t i : r)
           mysum += entrysize[i];
         partial_sums[ti.task_nr+1] = mysum;
       });

    for (size_t i = 1; i < partial_sums.Size(); i++)
      partial_sums[i] += partial_sums[i-1];

    ParallelJob
      ([&] (TaskInfo ti)
       {
         IntRange r = IntRange(size).Split(ti.task_nr, ti.ntasks);
         size_t mysum = partial_sums[ti.task_nr];
         for (size_t i : r)
           {
             index[i] = mysum;
             mysum += entrysize[i];
           }
       });
    index[size] = partial_sums.Last();

    return index;
  }

  NGCORE_API size_t * TablePrefixSum32 (FlatArray<unsigned int> entrysize)
  { return TablePrefixSum2 (entrysize); }
  NGCORE_API size_t * TablePrefixSum64 (FlatArray<size_t> entrysize)
  { return TablePrefixSum2 (entrysize); }

  /*
  BaseDynamicTable :: BaseDynamicTable (int size)
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

  BaseDynamicTable :: BaseDynamicTable (const Array<int> & entrysizes, int elemsize)
    : data(entrysizes.Size())
  {
    int cnt = 0;
    int n = entrysizes.Size();

    for (int i = 0; i < n; i++)
      cnt += entrysizes[i];
    oneblock = new char[elemsize * cnt];

    cnt = 0;
    for (int i = 0; i < n; i++)
      {
	data[i].maxsize = entrysizes[i];
	data[i].size = 0;

	data[i].col = &oneblock[elemsize * cnt];
	cnt += entrysizes[i];
      }
  }

  BaseDynamicTable :: ~BaseDynamicTable ()
  {
    if (oneblock)
      delete [] oneblock;
    else
      for (int i = 0; i < data.Size(); i++)
	delete [] static_cast<char*> (data[i].col);
  }

  void BaseDynamicTable :: SetSize (int size)
  {
    for (int i = 0; i < data.Size(); i++)
      delete [] static_cast<char*> (data[i].col);

    data.SetSize(size);
    for (int i = 0; i < size; i++)
      {
	data[i].maxsize = 0;
	data[i].size = 0;
	data[i].col = NULL;
      }
  }

  void BaseDynamicTable :: IncSize (IndexType i, int elsize)
  {
    if (i < 0 || i >= data.Size())
      {
        std::cerr << "BaseDynamicTable::Inc: Out of range, i = " << i << ", size = " << data.Size() << std::endl;
	return;
      }

    linestruct & line = data[i];

    if (line.size == line.maxsize)
      {
	void * p = new char [(2*line.maxsize+5) * elsize];

	memcpy (p, line.col, line.maxsize * elsize);
	delete [] static_cast<char*> (line.col);
	line.col = p;
	line.maxsize = 2*line.maxsize+5;
      }

    line.size++;
  }

  void BaseDynamicTable :: DecSize (IndexType i)
  {
    if (i < 0 || i >= data.Size())
      {
        std::cerr << "BaseDynamicTable::Dec: Out of range" << std::endl;
	return;
      }

    linestruct & line = data[i];

    if (line.size == 0)
      {
        std::cerr << "BaseDynamicTable::Dec: EntrySize < 0" << std::endl;
	return;
      }

    line.size--;
  }
  */

  void FilteredTableCreator::Add (size_t blocknr, int data)
  {
    if (!takedofs||takedofs->Test(data))
      TableCreator<int>::Add(blocknr,data);
  }

  void FilteredTableCreator::Add (size_t blocknr, IntRange range)
  {
    for (size_t i=range.First(); i<range.Next();i++)
      if (!takedofs||takedofs->Test(i))
	TableCreator<int>::Add(blocknr,i);
  }  
  
  void FilteredTableCreator::Add (size_t blocknr, FlatArray<int> dofs)
  {
    for (size_t i = 0; i < dofs.Size(); i++)
      if (!takedofs||takedofs->Test(dofs[i]))
	TableCreator<int>::Add(blocknr,dofs[i]);
  }  

} // namespace ngcore
