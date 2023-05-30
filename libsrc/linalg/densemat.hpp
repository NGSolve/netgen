#ifndef FILE_DENSEMAT
#define FILE_DENSEMAT

/**************************************************************************/
/* File:   densemat.hpp                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Oct. 94                                                    */
/**************************************************************************/

/** 
    Data type dense matrix
*/


class DenseMatrix
{
protected:
  int height;
  int width;
  double * data;

public:
  ///
  DLL_HEADER DenseMatrix ();
  ///
  DLL_HEADER DenseMatrix (int h, int w = 0);
  ///
  DLL_HEADER DenseMatrix (const DenseMatrix & m2);
  ///
  DLL_HEADER ~DenseMatrix ();

  ///
  DLL_HEADER void SetSize (int h, int w = 0);

  int Height() const { return height; }
  int Width() const {return width; }

  double & operator() (int i, int j) { return data[i*width+j]; }
  double operator() (int i, int j) const { return data[i*width+j]; }
  double & operator() (int i) { return data[i]; }
  double operator() (int i) const { return data[i]; }

  ///
  DLL_HEADER DenseMatrix & operator= (const DenseMatrix & m2);
  ///
  DLL_HEADER DenseMatrix & operator+= (const DenseMatrix & m2);
  ///
  DLL_HEADER DenseMatrix & operator-= (const DenseMatrix & m2);

  ///
  DLL_HEADER DenseMatrix & operator= (double v);
  ///
  DLL_HEADER DenseMatrix & operator*= (double v);

  ///
  DLL_HEADER void Mult (const FlatVector & v, FlatVector & prod) const
  {
    double sum;
    const double * mp, * sp;
    double * dp;
    
#ifdef DEBUG
    if (prod.Size() != height)
      {
	(*myerr) << "Mult: wrong vector size " << endl;
      }
    if (!height) 
      {
	cout << "DenseMatrix::Mult height = 0" << endl;
      }
    if (!width) 
      {
	cout << "DenseMatrix::Mult width = 0" << endl;
      }
    
    if (width != v.Size())
      {
	(*myerr) << "\nMatrix and Vector don't fit" << endl;
      }
    else if (Height() != prod.Size())
      {
	(*myerr) << "Base_Matrix::operator*(Vector): prod vector not ok" << endl;
      }
    else
#endif
      {      
	mp = data;
	dp = &prod(0);
        for (int i = 0; i < height; i++)
	  {
	    sum = 0;
	    sp = &v(0);
	    
	    for (int j = 0; j < width; j++)
	      {
		//        sum += Get(i,j) * v.Get(j);
		sum += *mp * *sp;
		mp++;
		sp++;
	      }
	    
	    *dp = sum;
	    dp++;
	  }
      }
  }

  ///
  DLL_HEADER void MultTrans (const Vector & v, Vector & prod) const;
  ///
  DLL_HEADER void Residuum (const Vector & x, const Vector & b, Vector & res) const;
  ///
  DLL_HEADER double Det () const;

  ///
  friend DenseMatrix operator* (const DenseMatrix & m1, const DenseMatrix & m2);
  ///
  friend DenseMatrix operator+ (const DenseMatrix & m1, const DenseMatrix & m2);

  /// 
  friend void Transpose (const DenseMatrix & m1, DenseMatrix & m2);
  ///
  friend void Mult (const DenseMatrix & m1, const DenseMatrix & m2, DenseMatrix & m3);
  ///
//  friend void CalcInverse (const DenseMatrix & m1, DenseMatrix & m2);
  ///
  friend void CalcAAt (const DenseMatrix & a, DenseMatrix & m2);
  ///
//  friend void CalcAtA (const DenseMatrix & a, DenseMatrix & m2);
  ///
  friend void CalcABt (const DenseMatrix & a, const DenseMatrix & b, DenseMatrix & m2);
  ///
  friend void CalcAtB (const DenseMatrix & a, const DenseMatrix & b, DenseMatrix & m2);
  ///
  DLL_HEADER void Solve (const Vector & b, Vector & x) const;
  ///
  void SolveDestroy (const Vector & b, Vector & x);
  ///
  const double & Get(int i, int j) const { return data[(i-1)*width+j-1]; }
  ///
  const double & Get(int i) const { return data[i-1]; }
  ///
  void Set(int i, int j, double v) { data[(i-1)*width+j-1] = v; }
  ///
  double & Elem(int i, int j) { return data[(i-1)*width+j-1]; }
  ///
  const double & ConstElem(int i, int j) const { return data[(i-1)*width+j-1]; }
};


extern ostream & operator<< (ostream & ost, const DenseMatrix & m);



template <int WIDTH, typename T = double>
class MatrixFixWidth
{
protected:
  int height;
  T * __restrict data;
  bool ownmem;
public:
  ///
  MatrixFixWidth () 
  { height = 0; data = 0; ownmem = false; }
  ///
  MatrixFixWidth (int h)
  { height = h; data = new T[WIDTH*height]; ownmem = true; }
  /// 
  MatrixFixWidth (int h, T * adata)
  { height = h; data = adata; ownmem = false; }
  ///
  MatrixFixWidth (const MatrixFixWidth & m2)
    : height(m2.height), ownmem(true)
  {
    data = new T[height*WIDTH];
    for (int i = 0; i < WIDTH*height; i++)
      data[i] = m2.data[i];
  }
  // : height(m2.height), data(m2.data), ownmem(false)
  //{ ; } 
  ///
  ~MatrixFixWidth ()
  { if (ownmem) delete [] data; }

  void SetSize (int h)
  {
    if (h != height)
      {
	if (ownmem) delete data;
	height = h;
	data = new T[WIDTH*height]; 
	ownmem = true;
      }
  }

  ///
  int Height() const { return height; }

  ///
  int Width() const { return WIDTH; }
  MatrixFixWidth & operator= (const MatrixFixWidth & m2)
  {
    for (int i = 0; i < height*WIDTH; i++)
      data[i] = m2.data[i];
  }
  ///
  MatrixFixWidth & operator= (T v)
  {
    for (int i = 0; i < height*WIDTH; i++)
      data[i] = v; 
    return *this;
  }

  /*  
  ///
  void Mult (const FlatVector & v, FlatVector & prod) const
  {
    T sum;
    const T * mp, * sp;
    T * dp;

    mp = data;
    dp = &prod[0];
    for (int i = 0; i < height; i++)
      {
	sum = 0;
	sp = &v[0];
	
	for (int j = 0; j < WIDTH; j++)
	  {
	    sum += *mp * *sp;
	    mp++;
	    sp++;
	  }
	    
	*dp = sum;
	dp++;
      }
  }
  */

  T & operator() (int i, int j)
  { return data[i*WIDTH+j]; }

  const T & operator() (int i, int j) const
  { return data[i*WIDTH+j]; }


  MatrixFixWidth & operator*= (T v)
  {
    if (data)
      for (int i = 0; i < height*WIDTH; i++)
        data[i] *= v;
    return *this;
  }



  const T & Get(int i, int j) const { return data[(i-1)*WIDTH+j-1]; }
  ///
  const T & Get(int i) const { return data[i-1]; }
  ///
  void Set(int i, int j, T v) { data[(i-1)*WIDTH+j-1] = v; }
  ///
  T & Elem(int i, int j) { return data[(i-1)*WIDTH+j-1]; }
  ///
  const T & ConstElem(int i, int j) const { return data[(i-1)*WIDTH+j-1]; }
};


template <int WIDTH>
class MatrixFixWidth<WIDTH,double>
{
protected:
  int height;
  double * data;
  bool ownmem;
public:
  ///
  MatrixFixWidth () 
  { height = 0; data = 0; ownmem = false; }
  ///
  MatrixFixWidth (int h)
  { height = h; data = new double[WIDTH*height]; ownmem = true; }
  
  MatrixFixWidth (const MatrixFixWidth & m2)
    : height(m2.height), ownmem(true)
  {
    data = new double[height*WIDTH];
    for (int i = 0; i < WIDTH*height; i++)
      data[i] = m2.data[i];
  }

  /// 
  MatrixFixWidth (int h, double * adata)
  { height = h; data = adata; ownmem = false; }
  ///
  ~MatrixFixWidth ()
  { if (ownmem) delete [] data; }

  void SetSize (int h)
  {
    if (h != height)
      {
	if (ownmem) delete data;
	height = h;
	data = new double[WIDTH*height]; 
	ownmem = true;
      }
  }

  ///
  int Height() const { return height; }

  ///
  int Width() const { return WIDTH; }

  ///
  MatrixFixWidth & operator= (double v)
  {
    for (int i = 0; i < height*WIDTH; i++)
      data[i] = v; 
    return *this;
  }

  ///
  void Mult (const FlatVector & v, FlatVector & prod) const
  {
    double sum;
    const double * mp, * sp;
    double * dp;

    /*    
    if (prod.Size() != height)
      {
	cerr << "MatrixFixWidth::Mult: wrong vector size " << endl;
	assert (1);
      }
    */    

    mp = data;
    dp = &prod[0];
    for (int i = 0; i < height; i++)
      {
	sum = 0;
	sp = &v[0];
	
	for (int j = 0; j < WIDTH; j++)
	  {
	    sum += *mp * *sp;
	    mp++;
	    sp++;
	  }
	    
	*dp = sum;
	dp++;
      }
  }

  double & operator() (int i, int j)
  { return data[i*WIDTH+j]; }

  const double & operator() (int i, int j) const
  { return data[i*WIDTH+j]; }


  MatrixFixWidth & operator*= (double v)
  {
    if (data)
      for (int i = 0; i < height*WIDTH; i++)
        data[i] *= v;
    return *this;
  }



  const double & Get(int i, int j) const { return data[(i-1)*WIDTH+j-1]; }
  ///
  const double & Get(int i) const { return data[i-1]; }
  ///
  void Set(int i, int j, double v) { data[(i-1)*WIDTH+j-1] = v; }
  ///
  double & Elem(int i, int j) { return data[(i-1)*WIDTH+j-1]; }
  ///
  const double & ConstElem(int i, int j) const { return data[(i-1)*WIDTH+j-1]; }
};





template <int WIDTH>
extern ostream & operator<< (ostream & ost, const MatrixFixWidth<WIDTH> & m)
{
  for (int i = 0; i < m.Height(); i++)
    {
      for (int j = 0; j < m.Width(); j++)
	ost << m.Get(i+1,j+1) << " ";
      ost << endl;
    }
  return ost;
};


extern DLL_HEADER void CalcAtA (const DenseMatrix & a, DenseMatrix & m2);
extern DLL_HEADER void CalcInverse (const DenseMatrix & m1, DenseMatrix & m2);


#endif
