/**************************************************************************/
/* File:   sort.cc                                                        */
/* Author: Joachim Schoeberl                                              */
/* Date:   07. Jan. 00                                                    */
/**************************************************************************/

/* 
   Sorting
*/


#include <algorithm>
#include <mystdlib.h>
#include <myadt.hpp>

namespace netgen
{

  void Sort (FlatArray<double> values, Array<int> & order)
  {
    int n = values.Size();
    int i, j;

    order.SetSize (n);

    for (i = 1; i <= n; i++)
      order[i-1] = i;
    for (i = 1; i <= n-1; i++)
      for (j = 1; j <= n-1; j++)
        if (values[order[j-1]-1] > values[order[j]-1])
          {
            Swap (order[j-1], order[j]);
          }
  }


  void QuickSortRec (FlatArray<double> values, FlatArray<int> order,
                     int left, int right)
  {
    int i, j;
    double midval;

    i = left;
    j = right;
    midval = values[order[(i+j)/2-1]-1];
  
    do
      {
        while (values[order[i-1]-1] < midval) i++;
        while (midval < values[order[j-1]-1]) j--;
      
        if (i <= j)
          {
            Swap (order[i-1], order[j-1]);
            i++; j--;
          }
      }
    while (i <= j);
    if (left < j) QuickSortRec (values, order, left, j);
    if (i < right) QuickSortRec (values, order, i, right);
  }

  void QuickSort (FlatArray<double> values, Array<int> & order)
  {
    int i, n = values.Size();
    order.SetSize (n);
    for (i = 1; i <= n; i++)
      order[i-1] = i;

    QuickSortRec (values, order, 1, order.Size());
  }
}
