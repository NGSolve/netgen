#ifndef FILE_VECTOR
#define FILE_VECTOR

/* *************************************************************************/
/* File:   vector.hpp                                                      */
/* Author: Joachim Schoeberl                                               */
/* Date:   01. Oct. 94                                                     */
/* *************************************************************************/

namespace netgen
{

template <typename T>
class TFlatVector
{
protected:
  int s;
  T * data;
public:
  TFlatVector () { ; }
  TFlatVector (int as, T * adata)
  { s = as; data = adata; }
    
  int Size () const
  { return s; }

  TFlatVector & operator= (const TFlatVector & v) 
  { memcpy (data, v.data, s*sizeof(T)); return *this; }

  TFlatVector & operator= (T scal) 
  {
    for (int i = 0; i < s; i++) data[i] = scal; 
    return *this;
  }

  T & operator[] (int i) { return data[i]; }
  const T & operator[] (int i) const { return data[i]; }
  T & operator() (int i) { return data[i]; }
  const T & operator() (int i) const { return data[i]; }

  // double & Elem (int i) { return data[i-1]; }
  // const double & Get (int i) const { return data[i-1]; }
  // void Set (int i, double val) { data[i-1] = val; }

  TFlatVector & operator*= (T scal)
  {
    for (int i = 0; i < s; i++) data[i] *= scal;
    return *this;
  }
};




class FlatVector
{
protected:
  int s;
  double *data;
public:
  FlatVector () { ; }
  FlatVector (int as, double * adata)
  { s = as; data = adata; }

  int Size () const
  { return s; }

  FlatVector & operator= (const FlatVector & v) 
  { memcpy (data, v.data, s*sizeof(double)); return *this; }

  FlatVector & operator= (double scal) 
  {
    for (int i = 0; i < s; i++) data[i] = scal; 
    return *this;
  }

  double & operator[] (int i) { return data[i]; }
  const double & operator[] (int i) const { return data[i]; }
  double & operator() (int i) { return data[i]; }
  const double & operator() (int i) const { return data[i]; }

  // double & Elem (int i) { return data[i-1]; }
  // const double & Get (int i) const { return data[i-1]; }
  // void Set (int i, double val) { data[i-1] = val; }

  FlatVector & operator*= (double scal)
  {
    for (int i = 0; i < s; i++) data[i] *= scal;
    return *this;
  }

  FlatVector & Add (double scal, const FlatVector & v2)
  {
    for (int i = 0; i < s; i++) 
      data[i] += scal * v2[i];
    return *this;
  }

  FlatVector & Set (double scal, const FlatVector & v2)
  {
    for (int i = 0; i < s; i++) 
      data[i] = scal * v2[i];
    return *this;
  }

  FlatVector & Set2 (double scal1, const FlatVector & v1,
		 double scal2, const FlatVector & v2)
  {
    for (int i = 0; i < s; i++) 
      data[i] = scal1 * v1[i] + scal2 * v2[i];
    return *this;
  }
  
  double L2Norm() const
  {
    double sum = 0;
    for (int i = 0; i < s; i++)
      sum += data[i] * data[i];
    return sqrt (sum);
  }

  operator TFlatVector<double> () const { return TFlatVector<double> (s, data); } 
  friend double operator* (const FlatVector & v1, const FlatVector & v2);
};




class Vector : public FlatVector
{
  bool ownmem;
public:
  Vector () 
  { s = 0; data = 0; ownmem = false; }
  Vector (int as)
  { s = as; data = new double[s]; ownmem = true; }
  Vector (int as, double * mem)
  { s = as; data = mem; ownmem = false; }
  ~Vector ()
  { if (ownmem) delete [] data; }

  template<typename ARCHIVE>
  void DoArchive(ARCHIVE& ar)
  {
    auto size = s;
    ar & ownmem & size;
    if(!ar.Output())
      SetSize(size);
    ar.Do(data, size);
  }
  Vector & operator= (const FlatVector & v) 
  { memcpy (data, &v(0), s*sizeof(double)); return *this; }

  Vector & operator= (double scal) 
  {
    for (int i = 0; i < s; i++) data[i] = scal; 
    return *this;
  }

  void SetSize (int as)
  {
    if (s != as)
      {
	s = as;
	if (ownmem) delete [] data;
	data = new double [s];
        ownmem = true;
      }
  }

  operator TFlatVector<double> () const { return TFlatVector<double> (s, data); } 
};

template <int S>
class VectorMem : public Vector
{
  double mem[S];
public:
  VectorMem () : Vector(S, &mem[0]) { ; }

  VectorMem & operator= (const FlatVector & v) 
  { memcpy (data, &v(0), S*sizeof(double)); return *this; }

  VectorMem & operator= (double scal) 
  {
    for (int i = 0; i < S; i++) data[i] = scal; 
    return *this;
  }
};





inline double operator* (const FlatVector & v1, const FlatVector & v2)
{
  double sum = 0;
  for (int i = 0; i < v1.s; i++)
    sum += v1.data[i] * v2.data[i];
  return sum;
}




inline ostream & operator<< (ostream & ost, const FlatVector & v)
{
  for (int i = 0; i < v.Size(); i++)
    ost << " " << setw(7) << v[i];
  return ost;
}

} //namespace netgen

#endif


