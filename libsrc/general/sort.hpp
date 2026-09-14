#ifndef FILE_SORT
#define FILE_SORT

/**************************************************************************/
/* File:   sort.hh                                                        */
/* Author: Joachim Schoeberl                                              */
/* Date:   07. Jan. 00                                                    */
/**************************************************************************/

namespace netgen
{

// order(i) is sorted index of element i
extern void Sort (FlatArray<double> values, Array<int> & order);

extern void QuickSort (FlatArray<double> values, Array<int> & order);




template <class T>
inline void BubbleSort (int size, T * data)
{
  T hv;
  for (int i = 0; i < size; i++)
    for (int j = i+1; j < size; j++)
      if (data[i] > data[j])
	{
	  hv = data[i];
	  data[i] = data[j];
	  data[j] = hv;
	}
}

template <class T, class S>
inline void QuickSortPairRec (FlatArray<T> data, FlatArray<S> index,
			  int left, int right)
{
  int i = left;
  int j = right;
  T midval = data[(left+right)/2];

  do
    {
      while (data[i] < midval) i++;
      while (midval < data[j]) j--;

      if (i <= j)
	{
	  ngcore::Swap (data[i], data[j]);
	  ngcore::Swap (index[i], index[j]);
	  i++; j--;
	}
    }
  while (i <= j);
  if (left < j) QuickSortPairRec (data, index, left, j);
  if (i < right) QuickSortPairRec (data, index, i, right);
}

template <class T, class S>
inline void QuickSortPair (FlatArray<T> data, FlatArray<S> index)
{
  if (data.Size() > 1)
    QuickSortPairRec (data, index, 0, data.Size()-1);
}

  template <class T> 
  void Intersection (FlatArray<T> in1, FlatArray<T> in2,
		     Array<T> & out)
  {
    out.SetSize(0);
    for(int i=0; i<in1.Size(); i++)
      if(in2.Contains(in1[i]))
	out.Append(in1[i]);
  }
  template <class T> 
  void Intersection (FlatArray<T> in1, FlatArray<T> in2, FlatArray<T> in3,
		     Array<T> & out)
  {
    out.SetSize(0);
    for(int i=0; i<in1.Size(); i++)
      if(in2.Contains(in1[i]) && in3.Contains(in1[i]))
	out.Append(in1[i]);
  }

}

#endif
