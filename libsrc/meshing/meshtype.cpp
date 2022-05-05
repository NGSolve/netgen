#include <mystdlib.h>

#include "meshing.hpp"  

namespace netgen
{
  int MultiPointGeomInfo :: 
  AddPointGeomInfo (const PointGeomInfo & gi)
  {
    for (auto & pgi : mgi)
      if (pgi.trignum == gi.trignum)
	return 0;
  
    mgi.Append(gi);
    return 0;
  }
  


#ifdef PARALLEL
  MPI_Datatype MeshPoint :: MyGetMPIType ( )
  { 
    static MPI_Datatype type = MPI_DATATYPE_NULL;
    static MPI_Datatype htype = MPI_DATATYPE_NULL;
    if (type == MPI_DATATYPE_NULL)
      {
	MeshPoint hp;
	int blocklen[] = { 3, 1, 1 };
	MPI_Aint displ[] = { (char*)&hp.x[0] - (char*)&hp,
			     (char*)&hp.layer - (char*)&hp,
			     (char*)&hp.singular - (char*)&hp };
	MPI_Datatype types[] = { MPI_DOUBLE, MPI_INT, MPI_DOUBLE };
	// *testout << "displ = " << displ[0] << ", " << displ[1] << ", " << displ[2] << endl;
	// *testout << "sizeof = " << sizeof (MeshPoint) << endl;
	MPI_Type_create_struct (3, blocklen, displ, types, &htype);
	MPI_Type_commit ( &htype );
	MPI_Aint lb, ext;
	MPI_Type_get_extent (htype, &lb, &ext);
	// *testout << "lb = " << lb << endl;
	// *testout << "ext = " << ext << endl;
	ext = sizeof (MeshPoint);
	MPI_Type_create_resized (htype, lb, ext, &type);
	MPI_Type_commit ( &type );
	
      }
    return type;
  }


  MPI_Datatype Element2d :: MyGetMPIType ( )
  { 
    static MPI_Datatype type = MPI_DATATYPE_NULL;
    static MPI_Datatype htype = MPI_DATATYPE_NULL;
    if (type == MPI_DATATYPE_NULL)
      {
	Element2d hel;
	int blocklen[] = { ELEMENT2D_MAXPOINTS, 1, 1, 1 };
	MPI_Aint displ[] =
          { (char*)&hel.pnum[0] - (char*)&hel,
            (char*)&hel.index - (char*)&hel,
            (char*)&hel.typ - (char*)&hel,
            (char*)&hel.np - (char*)&hel
          };
	MPI_Datatype types[] = { GetMPIType<PointIndex>(), GetMPIType(hel.index),
                                 GetMPIType(hel.typ), GetMPIType(hel.np) };
	// *testout << "displ = " << displ[0] << ", " << displ[1] << ", " << displ[2] << endl;
	// *testout << "sizeof = " << sizeof (MeshPoint) << endl;
	MPI_Type_create_struct (4, blocklen, displ, types, &htype);
	MPI_Type_commit ( &htype );
	MPI_Aint lb, ext;
	MPI_Type_get_extent (htype, &lb, &ext);
	// *testout << "lb = " << lb << endl;
	// *testout << "ext = " << ext << endl;
	ext = sizeof (Element2d);
	MPI_Type_create_resized (htype, lb, ext, &type);
	MPI_Type_commit ( &type );
      }
    return type;
  }

  MPI_Datatype Element :: MyGetMPIType ( )
  {
    static MPI_Datatype type = MPI_DATATYPE_NULL;
    static MPI_Datatype htype = MPI_DATATYPE_NULL;
    if (type == MPI_DATATYPE_NULL)
      {
	Element hel;
	int blocklen[] = { ELEMENT_MAXPOINTS, 1, 1, 1 };
	MPI_Aint displ[] =
          { (char*)&hel.pnum[0] - (char*)&hel,
            (char*)&hel.index - (char*)&hel,
            (char*)&hel.typ - (char*)&hel,
            (char*)&hel.np - (char*)&hel
          };
	MPI_Datatype types[] = { GetMPIType<PointIndex>(), GetMPIType(hel.index),
                                 GetMPIType(hel.typ), GetMPIType(hel.np) };
	// *testout << "displ = " << displ[0] << ", " << displ[1] << ", " << displ[2] << endl;
	// *testout << "sizeof = " << sizeof (MeshPoint) << endl;
	MPI_Type_create_struct (4, blocklen, displ, types, &htype);
	MPI_Type_commit ( &htype );
	MPI_Aint lb, ext;
	MPI_Type_get_extent (htype, &lb, &ext);
	// *testout << "lb = " << lb << endl;
	// *testout << "ext = " << ext << endl;
	ext = sizeof (Element);
	MPI_Type_create_resized (htype, lb, ext, &type);
	MPI_Type_commit ( &type );
      }
    return type;
  }

  MPI_Datatype Segment :: MyGetMPIType ( )
  {
    static MPI_Datatype type = MPI_DATATYPE_NULL;
    static MPI_Datatype htype = MPI_DATATYPE_NULL;
    if (type == MPI_DATATYPE_NULL)
      {
	Segment hel;
	int blocklen[] = { 3, 1, 1, 1 };
	MPI_Aint displ[] =
          { (char*)&hel.pnums[0] - (char*)&hel,
            (char*)&hel.edgenr - (char*)&hel,
            (char*)&hel.cd2i - (char*)&hel,
            (char*)&hel.si - (char*)&hel
          };
	MPI_Datatype types[] = {
          GetMPIType<PointIndex>(), GetMPIType(hel.edgenr), GetMPIType(hel.cd2i), GetMPIType(hel.si)
        };
	// *testout << "displ = " << displ[0] << ", " << displ[1] << ", " << displ[2] << endl;
	// *testout << "sizeof = " << sizeof (MeshPoint) << endl;
	MPI_Type_create_struct (4, blocklen, displ, types, &htype);
	MPI_Type_commit ( &htype );
	MPI_Aint lb, ext;
	MPI_Type_get_extent (htype, &lb, &ext);
	// *testout << "lb = " << lb << endl;
	// *testout << "ext = " << ext << endl;
	ext = sizeof (Segment);
	MPI_Type_create_resized (htype, lb, ext, &type);
	MPI_Type_commit ( &type );
      }
    return type;
  }

#endif




  Segment :: Segment() 
    : is_curved(false)
  {
    pnums[0] = PointIndex::INVALID;
    pnums[1] = PointIndex::INVALID;
    edgenr = -1;

    singedge_left = 0.;
    singedge_right = 0.;
    seginfo = 0;

    si = -1;

    domin = -1;
    domout = -1;
    tlosurf = -1;

    surfnr1 = -1;
    surfnr2 = -1;
    pnums[2] = PointIndex::INVALID;
    meshdocval = 0;
    geominfo[0].trignum=-1;
    geominfo[1].trignum=-1;

    /*
      epgeominfo[0].edgenr = 1;
      epgeominfo[0].dist = 0;
      epgeominfo[1].edgenr = 1;
      epgeominfo[1].dist = 0;
    */
  }    

  Segment::Segment (const Segment & other)
    : 
    edgenr(other.edgenr),
    singedge_left(other.singedge_left),
    singedge_right(other.singedge_right),
    seginfo(other.seginfo),
    si(other.si),
    domin(other.domin),
    domout(other.domout),
    tlosurf(other.tlosurf),
    geominfo(),
    surfnr1(other.surfnr1),
    surfnr2(other.surfnr2),
    epgeominfo(),
    meshdocval(other.meshdocval),
    is_curved(other.is_curved),
    hp_elnr(other.hp_elnr)
  {
    for (int j = 0; j < 3; j++)
      pnums[j] = other.pnums[j];

    geominfo[0] = other.geominfo[0];
    geominfo[1] = other.geominfo[1];
    epgeominfo[0] = other.epgeominfo[0];
    epgeominfo[1] = other.epgeominfo[1];
  }

  Segment& Segment::operator=(const Segment & other)
  {
    if (&other != this)
      {
	pnums[0] = other[0];
	pnums[1] = other[1];
	edgenr = other.edgenr;
	singedge_left = other.singedge_left;
	singedge_right = other.singedge_right;
	seginfo = other.seginfo;
	si = other.si;
	domin = other.domin;
	domout = other.domout;
	tlosurf = other.tlosurf;
	geominfo[0] = other.geominfo[0];
	geominfo[1] = other.geominfo[1];
	surfnr1 = other.surfnr1;
	surfnr2 = other.surfnr2;
	epgeominfo[0] = other.epgeominfo[0];
	epgeominfo[1] = other.epgeominfo[1];
	pnums[2] = other.pnums[2];
	meshdocval = other.meshdocval;
	hp_elnr = other.hp_elnr;
        is_curved = other.is_curved;
      }
    
    return *this;
  }
  
  void Segment :: DoArchive (Archive & ar)
  {
    string * bcname_dummy = nullptr;
    ar & pnums[0] & pnums[1] & pnums[2]
      & edgenr & singedge_left & singedge_right
      & si & cd2i & domin & domout & tlosurf
      & surfnr1 & surfnr2
      & bcname_dummy // keep this for backward compatibility
      & epgeominfo[0].edgenr & epgeominfo[1].edgenr;
  }


  ostream & operator<<(ostream  & s, const Segment & seg)
  {
    s << seg[0] << "(gi=" << seg.geominfo[0].trignum << ") - "
      << seg[1] << "(gi=" << seg.geominfo[1].trignum << ")"
      << " domin = " << seg.domin << ", domout = " << seg.domout 
      << " si = " << seg.si << ", edgenr = " << seg.edgenr;
    return s;
  }

  // needed, e.g. for MPI communication
  Element2d :: Element2d ()
  {
    for (int i = 0; i < ELEMENT2D_MAXPOINTS; i++)
      {
	pnum[i] = 0;
	geominfo[i].trignum = 0;
      }
    np = 3;
    index = 0;
    badel = 0;
    deleted = 0;
    visible = 1;
    typ = TRIG;
    orderx = ordery = 1;
    refflag = 1;
    strongrefflag = false;
    is_curved = false;
  } 

  Element2d :: Element2d (int anp)
  { 
    for (int i = 0; i < ELEMENT2D_MAXPOINTS; i++)
      {
        pnum[i].Invalidate();
	geominfo[i].trignum = 0;
      }
    np = anp;
    index = 0;
    badel = 0;
    deleted = 0;
    visible = 1;
    switch (np)
      {
      case 3: typ = TRIG; break;
      case 4: typ = QUAD; break;
      case 6: typ = TRIG6; break;
      case 8: typ = QUAD8; break;
      }
    orderx = ordery = 1;
    refflag = 1;
    strongrefflag = false;
    is_curved = (np >= 4); // false;
  } 

  Element2d :: Element2d (ELEMENT_TYPE atyp)
  { 
    for (int i = 0; i < ELEMENT2D_MAXPOINTS; i++)
      {
        pnum[i].Invalidate();
	geominfo[i].trignum = 0;
      }

    SetType (atyp);

    index = 0;
    badel = 0;
    deleted = 0;
    visible = 1;
    orderx = ordery = 1;
    refflag = 1;
    strongrefflag = false;
    is_curved = (np >= 4); // false;
  } 



  Element2d :: Element2d (int pi1, int pi2, int pi3)
  {
    pnum[0] = pi1;
    pnum[1] = pi2;
    pnum[2] = pi3;
    np = 3;
    typ = TRIG;
    
    for (int i = 3; i < ELEMENT2D_MAXPOINTS; i++)
      pnum[i].Invalidate();
  
    for (int i = 0; i < ELEMENT2D_MAXPOINTS; i++)
      geominfo[i].trignum = 0;
    index = 0;
    badel = 0;
    refflag = 1;
    strongrefflag = false;
    deleted = 0;
    visible = 1;
    orderx = ordery = 1;
    is_curved = false;
  }

  Element2d :: Element2d (int pi1, int pi2, int pi3, int pi4)
  {
    pnum[0] = pi1;
    pnum[1] = pi2;
    pnum[2] = pi3;
    pnum[3] = pi4;
    np = 4;
    typ = QUAD;

    pnum[4] = 0;
    pnum[5] = 0;
  
    for (int i = 0; i < ELEMENT2D_MAXPOINTS; i++)
      geominfo[i].trignum = 0;
    index = 0;
    badel = 0;
    refflag = 1;
    strongrefflag = false;
    deleted = 0;
    visible = 1;
    orderx = ordery = 1;
    is_curved = true;
  }


  /*
    void Element2d :: SetType (ELEMENT_TYPE atyp)
    {
    typ = atyp;
    switch (typ)
    {
    case TRIG: np = 3; break;
    case QUAD: np = 4; break;
    case TRIG6: np = 6; break;
    case QUAD6: np = 6; break;
    default:
    PrintSysError ("Element2d::SetType, illegal type ", typ);
    }
    }
  */


  void Element2d :: GetBox (const T_POINTS & points, Box3d & box) const
  {
    box.SetPoint (points[pnum[0]]);
    for (unsigned i = 1; i < np; i++)
      box.AddPoint (points[pnum[i]]);
  }

  bool Element2d :: operator==(const Element2d & el2) const
  {
    bool retval = (el2.GetNP() == np);
    for(int i= 0; retval && i<np; i++)
      retval = (el2[i] == (*this)[i]);

    return retval;
  }


  void Element2d :: Invert2()
  {
    switch (typ)
      {
      case TRIG:
        {
          Swap (pnum[1], pnum[2]);
          break;
        }
      case TRIG6:
        {
          Swap (pnum[1], pnum[2]);
          Swap (pnum[4], pnum[5]);
          break;
        }
      case QUAD:
        {
          Swap (pnum[0], pnum[3]);
          Swap (pnum[1], pnum[2]);
          break;
        }
      default:
        {
          cerr << "Element2d::Invert2, illegal element type " << int(typ) << endl;
        }
      }
  }

  int Element2d::HasFace(const Element2d & el) const
  {
    //nur für tets!!! hannes
    for (int i = 1; i <= 3; i++)
      {
        if (PNumMod(i)   == el[0] && 
            PNumMod(i+1) == el[1] && 
            PNumMod(i+2) == el[2])
          {
            return 1;
          }
      }
    return 0;
  }

  void Element2d :: NormalizeNumbering2 ()
  {
    if (GetNP() == 3)
      {
        if (PNum(1) < PNum(2) && PNum(1) < PNum(3))
          return;
        else
          {
            if (PNum(2) < PNum(3))
              {
                PointIndex pi1 = PNum(2);
                PNum(2) = PNum(3);
                PNum(3) = PNum(1);
                PNum(1) = pi1;
              }
            else
              {
                PointIndex pi1 = PNum(3);
                PNum(3) = PNum(2);
                PNum(2) = PNum(1);
                PNum(1) = pi1;
              }
          }
      }
    else
      {
        int mini = 1;
        for (int i = 2; i <= GetNP(); i++)
          if (PNum(i) < PNum(mini)) mini = i;
      
        Element2d hel = (*this);
        for (int i = 1; i <= GetNP(); i++)
          PNum(i) = hel.PNumMod (i+mini-1);
      }
  }




  NgArray<IntegrationPointData*> ipdtrig;
  NgArray<IntegrationPointData*> ipdquad;


  int Element2d :: GetNIP () const
  {
    int nip;
    switch (np)
      {
      case 3: nip = 1; break;
      case 4: nip = 4; break;
      default: nip = 0; break;
      }
    return nip;
  }

  void Element2d :: 
  GetIntegrationPoint (int ip, Point<2> & p, double & weight) const
  {
    static double eltriqp[1][3] =
      {
        { 1.0/3.0, 1.0/3.0, 0.5 }
      };

    static double elquadqp[4][3] =
      { 
        { 0, 0, 0.25 },
        { 0, 1, 0.25 },
        { 1, 0, 0.25 },
        { 1, 1, 0.25 }
      };
  
    double * pp = 0;
    switch (typ)
      {
      case TRIG: pp = &eltriqp[0][0]; break;
      case QUAD: pp = &elquadqp[ip-1][0]; break;
      default:
        PrintSysError ("Element2d::GetIntegrationPoint, illegal type ", int(typ));
      }

    p[0] = pp[0];
    p[1] = pp[1];
    weight = pp[2];
  }

  void Element2d :: 
  GetTransformation (int ip, const NgArray<Point<2>> & points,
                     DenseMatrix & trans) const
  {
    int np = GetNP();
    DenseMatrix pmat(2, np), dshape(2, np);
    pmat.SetSize (2, np);
    dshape.SetSize (2, np);

    Point<2> p;
    double w;

    GetPointMatrix (points, pmat);
    GetIntegrationPoint (ip, p, w);
    GetDShape (p, dshape);
  
    CalcABt (pmat, dshape, trans);

    /*
      (*testout) << "p = " << p  << endl
      << "pmat = " << pmat << endl
      << "dshape = " << dshape << endl
      << "tans = " << trans << endl;
    */
  }

  void Element2d :: 
  GetTransformation (int ip, class DenseMatrix & pmat,
                     class DenseMatrix & trans) const
  {
    //  int np = GetNP();

#ifdef DEBUG
    if (pmat.Width() != np || pmat.Height() != 2)
      {
        (*testout) << "GetTransofrmation: pmat doesn't fit" << endl;
        return;
      }
#endif

    ComputeIntegrationPointData ();
    DenseMatrix * dshapep = NULL;
    switch (typ)
      {
      case TRIG: dshapep = &ipdtrig.Get(ip)->dshape; break;
      case QUAD: dshapep = &ipdquad.Get(ip)->dshape; break;
      default:
        PrintSysError ("Element2d::GetTransformation, illegal type ", int(typ));
      }
  
    CalcABt (pmat, *dshapep, trans);
  }


  void Element2d :: GetShape (const Point<2> & p, Vector & shape) const
  {
    if (shape.Size() != GetNP())
      {
        cerr << "Element::GetShape: Length not fitting" << endl;
        return;
      }

    switch (typ)
      {
      case TRIG:
        shape(0) = 1 - p[0] - p[1];
        shape(1) = p[0];
        shape(2) = p[1];
        break;
      case QUAD:
        shape(0) = (1-p[0]) * (1-p[1]);
        shape(1) = p[0] * (1-p[1]);
        shape(2) = p[0] * p[1];
        shape(3) = (1-p[0]) * p[1];
        break;
      default:
        PrintSysError ("Element2d::GetShape, illegal type ", int(typ));
      }
  }



  void Element2d :: GetShapeNew (const Point<2> & p, FlatVector & shape) const
  {
    switch (typ)
      {
      case TRIG:
        {
          shape(0) = p(0);
          shape(1) = p(1);
          shape(2) = 1-p(0)-p(1);
          break;
        }

      case QUAD:
        {
          shape(0) = (1-p(0))*(1-p(1));
          shape(1) =    p(0) *(1-p(1));
          shape(2) =    p(0) *   p(1) ;
          shape(3) = (1-p(0))*   p(1) ;
          break;
        }

      default:
        throw NgException ("illegal element type in GetShapeNew");
      }
  }

  template <typename T>
  void Element2d :: GetShapeNew (const Point<2,T> & p, TFlatVector<T> shape) const
  {
    switch (typ)
      {
      case TRIG:
        {
          shape(0) = p(0);
          shape(1) = p(1);
          shape(2) = 1-p(0)-p(1);
          break;
        }

      case QUAD:
        {
          shape(0) = (1-p(0))*(1-p(1));
          shape(1) =    p(0) *(1-p(1));
          shape(2) =    p(0) *   p(1) ;
          shape(3) = (1-p(0))*   p(1) ;
          break;
        }
      default:
        throw NgException ("illegal element type in GetShapeNew");
      }
  }









  void Element2d :: 
  GetDShape (const Point<2> & p, DenseMatrix & dshape) const
  {
#ifdef DEBUG
    if (dshape.Height() != 2 || dshape.Width() != np)
      {
        PrintSysError ("Element::DShape: Sizes don't fit");
        return;
      }
#endif

    switch (typ)
      {
      case TRIG:
        dshape.Elem(1, 1) = -1;
        dshape.Elem(1, 2) = 1;
        dshape.Elem(1, 3) = 0;
        dshape.Elem(2, 1) = -1;
        dshape.Elem(2, 2) = 0;
        dshape.Elem(2, 3) = 1;
        break;
      case QUAD:
        dshape.Elem(1, 1) = -(1-p[1]);
        dshape.Elem(1, 2) = (1-p[1]);
        dshape.Elem(1, 3) = p[1];
        dshape.Elem(1, 4) = -p[1];
        dshape.Elem(2, 1) = -(1-p[0]);
        dshape.Elem(2, 2) = -p[0];
        dshape.Elem(2, 3) = p[0];
        dshape.Elem(2, 4) = (1-p[0]);
        break;

      default:
        PrintSysError ("Element2d::GetDShape, illegal type ", int(typ));
      }
  }


  template <typename T>
  void Element2d :: 
  GetDShapeNew (const Point<2,T> & p, MatrixFixWidth<2,T> & dshape) const
  {
    switch (typ)
      {
      case TRIG:
        {
          dshape = T(0.0);
          dshape(0,0) = 1;
          dshape(1,1) = 1;
          dshape(2,0) = -1;
          dshape(2,1) = -1;
          break;
        }
      case QUAD:
        {
          dshape(0,0) = -(1-p(1));
          dshape(0,1) = -(1-p(0));

          dshape(1,0) =  (1-p(1));
          dshape(1,1) =  -p(0);

          dshape(2,0) = p(1);
          dshape(2,1) = p(0);

          dshape(3,0) = -p(1);
          dshape(3,1) = (1-p(0));
          break;
        }
      default:
        throw NgException ("illegal element type in GetDShapeNew");
      }
  }





  void Element2d :: 
  GetPointMatrix (const NgArray<Point<2>> & points,
                  DenseMatrix & pmat) const
  {
    int np = GetNP();

#ifdef DEBUG
    if (pmat.Width() != np || pmat.Height() != 2)
      {
        cerr << "Element::GetPointMatrix: sizes don't fit" << endl;
        return;
      }
#endif
  
    for (int i = 1; i <= np; i++)
      {
        const auto& p = points.Get(PNum(i));
        pmat.Elem(1, i) = p[0];
        pmat.Elem(2, i) = p[1];
      }
  }





  double Element2d :: CalcJacobianBadness (const NgArray<Point<2>> & points) const
  {
    int i, j;
    int nip = GetNIP();
    DenseMatrix trans(2,2);
    DenseMatrix pmat;
  
    pmat.SetSize (2, GetNP());
    GetPointMatrix (points, pmat);

    double err = 0;
    for (i = 1; i <= nip; i++)
      {
        GetTransformation (i, pmat, trans);

        // Frobenius norm
        double frob = 0;
        for (j = 1; j <= 4; j++)
          frob += sqr (trans.Get(j));
        frob = sqrt (frob);
        frob /= 2;

        double det = trans.Det();

        if (det <= 0)
          err += 1e12;
        else
          err += frob * frob / det;
      }

    err /= nip;
    return err;
  }



  static const int qip_table[4][4] =
    { { 0, 1, 0, 3 },
      { 0, 1, 1, 2 },
      { 3, 2, 0, 3 },
      { 3, 2, 1, 2 }
    };

  double Element2d :: 
  CalcJacobianBadnessDirDeriv (const NgArray<Point<2>> & points,
                               int pi, Vec<2> & dir, double & dd) const
  {
    if (typ == QUAD)
      {
        Mat<2,2> trans, dtrans;
        Mat<2,4> vmat, pmat;
      
        for (int j = 0; j < 4; j++)
          {
            const auto& p = points.Get( (*this)[j] );
            pmat(0, j) = p[0];
            pmat(1, j) = p[1];
          }

        vmat = 0.0;
        vmat(0, pi-1) = dir[0];
        vmat(1, pi-1) = dir[1];
      
        double err = 0;
        dd = 0;

        for (int i = 0; i < 4; i++)
          {
            int ix1 = qip_table[i][0];
            int ix2 = qip_table[i][1];
            int iy1 = qip_table[i][2];
            int iy2 = qip_table[i][3];
	      
            trans(0,0) = pmat(0, ix2) - pmat(0,ix1);
            trans(1,0) = pmat(1, ix2) - pmat(1,ix1);
            trans(0,1) = pmat(0, iy2) - pmat(0,iy1);
            trans(1,1) = pmat(1, iy2) - pmat(1,iy1);

            double det = trans(0,0)*trans(1,1)-trans(1,0)*trans(0,1);

            if (det <= 0)
              {
                dd = 0;
                return 1e12;
              }
	  
            dtrans(0,0) = vmat(0, ix2) - vmat(0,ix1);
            dtrans(1,0) = vmat(1, ix2) - vmat(1,ix1);
            dtrans(0,1) = vmat(0, iy2) - vmat(0,iy1);
            dtrans(1,1) = vmat(1, iy2) - vmat(1,iy1);


            // Frobenius norm
            double frob = 0;
            for (int j = 0; j < 4; j++) 
              frob += sqr (trans(j));
            frob = sqrt (frob);
	  
            double dfrob = 0;
            for (int j = 0; j < 4; j++)
              dfrob += trans(j) * dtrans(j);
            dfrob = dfrob / frob;
	  
            frob /= 2;      
            dfrob /= 2;
	  
	  
            // ddet = \sum_j det (m_j)   with m_j = trans, except col j = dtrans
            double ddet 
              = dtrans(0,0) * trans(1,1) - trans(0,1) * dtrans(1,0)
              + trans(0,0) * dtrans(1,1) - dtrans(0,1) * trans(1,0);
	  
            err += frob * frob / det;
            dd += (2 * frob * dfrob * det - frob * frob * ddet) / (det * det);
          }
      
        err /= 4;
        dd /= 4;
        return err;
      }

    int nip = GetNIP();
    DenseMatrix trans(2,2), dtrans(2,2);
    DenseMatrix pmat, vmat;
  
    pmat.SetSize (2, GetNP());
    vmat.SetSize (2, GetNP());

    GetPointMatrix (points, pmat);
  
    vmat = 0.0;
    vmat.Elem(1, pi) = dir[0];
    vmat.Elem(2, pi) = dir[1];


    double err = 0;
    dd = 0;

    for (int i = 1; i <= nip; i++)
      {
        GetTransformation (i, pmat, trans);
        GetTransformation (i, vmat, dtrans);

        // Frobenius norm
        double frob = 0;
        for (int j = 1; j <= 4; j++)
          frob += sqr (trans.Get(j));
        frob = sqrt (frob);
      
        double dfrob = 0;
        for (int j = 1; j <= 4; j++)
          dfrob += trans.Get(j) * dtrans.Get(j);
        dfrob = dfrob / frob;
      
        frob /= 2;      
        dfrob /= 2;
      
        double det = trans(0,0)*trans(1,1)-trans(1,0)*trans(0,1);

        // ddet = \sum_j det (m_j)   with m_j = trans, except col j = dtrans
        double ddet 
          = dtrans(0,0) * trans(1,1) - trans(0,1) * dtrans(1,0)
          + trans(0,0) * dtrans(1,1) - dtrans(0,1) * trans(1,0);

        if (det <= 0)
          err += 1e12;
        else
          {
            err += frob * frob / det;
            dd += (2 * frob * dfrob * det - frob * frob * ddet) / (det * det);
          }
      }

    err /= nip;
    dd /= nip;
    return err;
  }



  double Element2d :: 
  CalcJacobianBadness (const T_POINTS & points, const Vec<3> & n) const
  {
    int i, j;
    int nip = GetNIP();
    DenseMatrix trans(2,2);
    DenseMatrix pmat;
  
    pmat.SetSize (2, GetNP());

    Vec<3> t1, t2;
    t1 = n.GetNormal();
    t2 = Cross (n, t1);

    for (i = 1; i <= GetNP(); i++)
      {
        const auto& p = points[PNum(i)];
        pmat.Elem(1, i) = p[0] * t1(0) + p[1] * t1(1) + p[2] * t1(2);
        pmat.Elem(2, i) = p[0] * t2(0) + p[1] * t2(1) + p[2] * t2(2);
      }

    double err = 0;
    for (i = 1; i <= nip; i++)
      {
        GetTransformation (i, pmat, trans);

        // Frobenius norm
        double frob = 0;
        for (j = 1; j <= 4; j++)
          frob += sqr (trans.Get(j));
        frob = sqrt (frob);
        frob /= 2;

        double det = trans.Det();
        if (det <= 0)
          err += 1e12;
        else
          err += frob * frob / det;
      }

    err /= nip;
    return err;
  }




  void Element2d :: ComputeIntegrationPointData () const
  {
    switch (np)
      {
      case 3: if (ipdtrig.Size()) return; break;
      case 4: if (ipdquad.Size()) return; break;
      }

    for (int i = 1; i <= GetNIP(); i++)
      {
        IntegrationPointData * ipd = new IntegrationPointData;
        Point<2> hp;
        GetIntegrationPoint (i, hp, ipd->weight);
        ipd->p(0) = hp[0];
        ipd->p(1) = hp[1];
        ipd->p(2) = 0;

        ipd->shape.SetSize(GetNP());
        ipd->dshape.SetSize(2, GetNP());

        GetShape (hp, ipd->shape);
        GetDShape (hp, ipd->dshape);

        switch (np)
          {
          case 3: ipdtrig.Append (ipd); break;
          case 4: ipdquad.Append (ipd); break;
          }
      }
  }







  ostream & operator<<(ostream  & s, const Element0d & el)
  {
    s << el.pnum << ", index = " << el.index;
    return s;
  }


  ostream & operator<<(ostream  & s, const Element2d & el)
  {
    s << "np = " << el.GetNP();
    for (int j = 1; j <= el.GetNP(); j++)
      s << " " << el.PNum(j);
    return s;
  }


  ostream & operator<<(ostream  & s, const Element & el)
  {
    s << "np = " << el.GetNP();
    for (int j = 0; j < el.GetNP(); j++)
      s << " " << int(el[j]);
    return s;
  }

  /*
  Element :: Element ()
  {
    typ = TET;
    np = 4;
    for (int i = 0; i < ELEMENT_MAXPOINTS; i++)
      pnum[i] = 0;
    index = 0;
    flags.marked = 1;
    flags.badel = 0;
    flags.reverse = 0;
    flags.illegal = 0;
    flags.illegal_valid = 0;
    flags.badness_valid = 0;
    flags.refflag = 1;
    flags.strongrefflag = false;
    flags.deleted = 0;
    flags.fixed = 0;
    orderx = ordery = orderz = 1;
    is_curved = false;
#ifdef PARALLEL
    partitionNumber = -1;
#endif
  }
  */

  Element :: Element (int anp)
  {
    np = anp;
    for (int i = 0; i < ELEMENT_MAXPOINTS; i++)
        pnum[i].Invalidate();
    index = 0;
    flags.marked = 1;
    flags.badel = 0;
    flags.reverse = 0;
    flags.illegal = 0;
    flags.illegal_valid = 0;
    flags.badness_valid = 0;
    flags.refflag = 1;
    flags.strongrefflag = false;
    flags.deleted = 0;
    flags.fixed = 0;

    switch (np)
      {
      case 4: typ = TET; break;
      case 5: typ = PYRAMID; break;
      case 6: typ = PRISM; break;
      case 8: typ = HEX; break;
      case 10: typ = TET10; break;
      case 13: typ = PYRAMID13; break;
      case 15: typ = PRISM15; break;
      case 20: typ = HEX20; break;
      default: cerr << "Element::Element: unknown element with " << np << " points" << endl;
      }
    orderx = ordery = orderz = 1;
    is_curved = typ != TET; // false;
  }

  void Element :: SetOrder (const int aorder) 
  { 
    orderx = aorder; 
    ordery = aorder; 
    orderz = aorder;
  }


  void Element :: SetOrder (const int ox, const int oy, const int oz) 
  { 
    orderx = ox; 
    ordery = oy;
    orderz = oz; 
  }


  Element :: Element (ELEMENT_TYPE type)
  {
    SetType (type);

    for (int i = 0; i < ELEMENT_MAXPOINTS; i++)
        pnum[i].Invalidate();
    index = 0;
    flags.marked = 1;
    flags.badel = 0;
    flags.reverse = 0;
    flags.illegal = 0;
    flags.illegal_valid = 0;
    flags.badness_valid = 0;
    flags.refflag = 1;
    flags.strongrefflag = false;
    flags.deleted = 0;
    flags.fixed = 0;
    orderx = ordery = orderz = 1;
    is_curved =  typ != TET; // false;
    // #ifdef PARALLEL
    // partitionNumber = -1;
    // #endif
  }




  /*
  Element & Element :: operator= (const Element & el2)
  {
    typ = el2.typ;
    np = el2.np;
    for (int i = 0; i < ELEMENT_MAXPOINTS; i++)
      pnum[i] = el2.pnum[i];
    index = el2.index;
    flags = el2.flags;
    orderx = el2.orderx;
    ordery = el2.ordery;
    orderz = el2.orderz;
    hp_elnr = el2.hp_elnr;
    flags = el2.flags;
    is_curved = el2.is_curved;
    return *this;
  }
  */


  void Element :: SetNP (int anp)
  {
    np = anp; 
    switch (np)
      {
      case 4: typ = TET; break;
      case 5: typ = PYRAMID; break;
      case 6: typ = PRISM; break;
      case 8: typ = HEX; break;
      case 10: typ = TET10; break;
      case 13: typ = PYRAMID13; break;
      case 15: typ = PRISM15; break;
      case 20: typ = HEX20; break;
        // 
      default: break;
        cerr << "Element::SetNP unknown element with " << np << " points" << endl;
      }
  }



  void Element :: SetType (ELEMENT_TYPE atyp)
  {
    typ = atyp;
    switch (atyp)
      {
      case TET: np = 4; break;
      case PYRAMID: np = 5; break;
      case PRISM: np = 6; break;
      case HEX: np = 8; break;
      case TET10: np = 10; break;
      case PYRAMID13: np = 13; break;
      case PRISM12: np = 12; break;
      case PRISM15: np = 15; break;
      case HEX20: np = 20; break;

      default: break;
        cerr << "Element::SetType unknown type  " << int(typ) << endl;
      }
    is_curved = (np > 4); 
  }



  void Element :: Invert()
  {
    switch (GetNP())
      {
      case 4:
        {
          Swap (PNum(3), PNum(4));
          break;
        }
      case 5:
        {
          Swap (PNum(1), PNum(4));
          Swap (PNum(2), PNum(3));
          break;
        }
      case 6:
        {
          Swap (PNum(1), PNum(4));
          Swap (PNum(2), PNum(5));
          Swap (PNum(3), PNum(6));
          break;
        }
      }
  }


  void Element :: Print (ostream & ost) const
  {
    ost << np << " Points: ";
    for (int i = 1; i <= np; i++)
      ost << pnum[i-1] << " " << endl;
  }

  void Element :: GetBox (const T_POINTS & points, Box3d & box) const
  {
    box.SetPoint (points[PNum(1)]);
    box.AddPoint (points[PNum(2)]);
    box.AddPoint (points[PNum(3)]);
    box.AddPoint (points[PNum(4)]);
  }

  double Element :: Volume (const T_POINTS & points) const
  {
    Vec<3> v1 = points[PNum(2)] - points[PNum(1)];
    Vec<3> v2 = points[PNum(3)] - points[PNum(1)];
    Vec<3> v3 = points[PNum(4)] - points[PNum(1)]; 
  
    return -(Cross (v1, v2) * v3) / 6;	 
  }  


  void Element :: GetFace2 (int i, Element2d & face) const
  {
    static const int tetfaces[][5] = 
      { { 3, 2, 3, 4, 0 },
        { 3, 3, 1, 4, 0 },
        { 3, 1, 2, 4, 0 },
        { 3, 2, 1, 3, 0 } };

    static const int tet10faces[][7] = 
      { { 3, 2, 3, 4, 10, 9, 8 },
        { 3, 3, 1, 4, 7, 10, 6 },
        { 3, 1, 2, 4, 9, 7, 5 },
        { 3, 2, 1, 3, 6, 8, 5 } };

    static const int pyramidfaces[][5] =
      { { 4, 1, 4, 3, 2 },
        { 3, 1, 2, 5, 0 },
        { 3, 2, 3, 5, 0 },
        { 3, 3, 4, 5, 0 },
        { 3, 4, 1, 5, 0 } };
    
    static const int prismfaces[][5] =
      {
        { 3, 1, 3, 2, 0 },
        { 3, 4, 5, 6, 0 },
        { 4, 1, 2, 5, 4 },
        { 4, 2, 3, 6, 5 },
        { 4, 3, 1, 4, 6 }
      };
    
    static const int hexfaces[][5] =
      {
        { 4, 4, 3, 2, 1 },
        { 4, 3, 7, 6, 2 },
        { 4, 7, 8, 5, 6 },
        { 4, 8, 4, 1, 5 },
        { 4, 1, 2, 6, 5 },
        { 4, 3, 4, 8, 7 }
      };
    

    switch (np)
      {
      case 4: // tet
        {
          face.SetType(TRIG);
          for (int j = 1; j <= 3; j++)
            face.PNum(j) = PNum(tetfaces[i-1][j]);
          break;
        }

      case 10: // tet10
        {
          face.SetType(TRIG6);
          for (int j = 1; j <= 6; j++)
            face.PNum(j) = PNum(tet10faces[i-1][j]);
          break;
        }

      case 5: // pyramid
        {
          // face.SetNP(pyramidfaces[i-1][0]);
          face.SetType ( (i == 1) ? QUAD : TRIG);
          for (int j = 1; j <= face.GetNP(); j++)
            face.PNum(j) = PNum(pyramidfaces[i-1][j]);
          break;
        }
      case 6: // prism
        {
          //	face.SetNP(prismfaces[i-1][0]);
          face.SetType ( (i >= 3) ? QUAD : TRIG);
          for (int j = 1; j <= face.GetNP(); j++)
            face.PNum(j) = PNum(prismfaces[i-1][j]);
          break;
        }
      case 8:
        {
          face.SetType(QUAD);
          for (int j = 1; j <= 4; j++)
            face.PNum(j) = PNum(hexfaces[i-1][j]);
          break;
        }
      }
  }



  void Element :: GetTets (NgArray<Element> & locels) const
  {
    GetTetsLocal (locels);
    int i, j;
    for (i = 1; i <= locels.Size(); i++)
      for (j = 1; j <= 4; j++)
        locels.Elem(i).PNum(j) = PNum ( locels.Elem(i).PNum(j) );
  }

  void Element :: GetTetsLocal (NgArray<Element> & locels) const
  {
    int i, j;
    locels.SetSize(0);
    switch (GetType())
      {
      case TET:
        {
          int linels[1][4] = 
            { { 1, 2, 3, 4 },
            };
          for (i = 0; i < 1; i++)
            {
              Element tet(4);
              for (j = 1; j <= 4; j++)
                tet.PNum(j) = linels[i][j-1];
              locels.Append (tet);
            }
          break;
        }
      case TET10:
        {
          int linels[8][4] = 
            { { 1, 5, 6, 7 },
              { 5, 2, 8, 9 },
              { 6, 8, 3, 10 },
              { 7, 9, 10, 4 },
              { 5, 6, 7, 9 },
              { 5, 6, 9, 8 },
              { 6, 7, 9, 10 },
              { 6, 8, 10, 9 } };
          for (i = 0; i < 8; i++)
            {
              Element tet(4);
              for (j = 1; j <= 4; j++)
                tet.PNum(j) = linels[i][j-1];
              locels.Append (tet);
            }
          break;
        }
      case PYRAMID:
        {
          int linels[2][4] = 
            { { 1, 2, 3, 5 },
              { 1, 3, 4, 5 } };
          for (i = 0; i < 2; i++)
            {
              Element tet(4);
              for (j = 1; j <= 4; j++)
                tet.PNum(j) = linels[i][j-1];
              locels.Append (tet);
            }
          break;
        }
      case PRISM:
      case PRISM12:
        {
          int linels[3][4] = 
            { { 1, 2, 3, 4 },
              { 4, 2, 3, 5 },
              { 6, 5, 4, 3 }
            };
          for (i = 0; i < 3; i++)
            {
              Element tet(4);
              for (j = 0; j < 4; j++)
                tet[j] = linels[i][j];
              locels.Append (tet);
            }
          break;
        }
      case HEX:
        {
          int linels[6][4] = 
            { { 1, 7, 2, 3 },
              { 1, 7, 3, 4 },
              { 1, 7, 4, 8 },
              { 1, 7, 8, 5 },
              { 1, 7, 5, 6 },
              { 1, 7, 6, 2 }
            };
          for (i = 0; i < 6; i++)
            {
              Element tet(4);
              for (j = 0; j < 4; j++)
                tet[j] = linels[i][j];
              locels.Append (tet);
            }
          break;
        }
      default:
        {
          cerr << "GetTetsLocal not implemented for el with " << GetNP() << " nodes" << endl;
        }
      }
  }

  bool Element :: operator==(const Element & el2) const
  {
    bool retval = (el2.GetNP() == np);
    for(int i= 0; retval && i<np; i++)
      retval = (el2[i] == (*this)[i]);

    return retval;
  }


#ifdef OLD
  void Element :: GetNodesLocal (NgArray<Point3d> & points) const
  {
    const static double tetpoints[4][3] =
      { { 0, 0, 0 },
        { 1, 0, 0 },
        { 0, 1, 0 },
        { 0, 0, 1 }};
  
    const static double prismpoints[6][3] =
      { { 0, 0, 0 },
        { 1, 0, 0 },
        { 0, 1, 0 },
        { 0, 0, 1 },
        { 1, 0, 1 },
        { 0, 1, 1 } };
  
    const static double pyramidpoints[6][3] =
      { { 0, 0, 0 },
        { 1, 0, 0 },
        { 1, 1, 0 },
        { 0, 1, 0 },
        { 0, 0, 1 } };
  
    const static double tet10points[10][3] =
      { { 0, 0, 0 },
        { 1, 0, 0 },
        { 0, 1, 0 },
        { 0, 0, 1 },
        { 0.5, 0, 0 },
        { 0, 0.5, 0 },
        { 0, 0, 0.5 },
        { 0.5, 0.5, 0 },
        { 0.5, 0, 0.5 },
        { 0, 0.5, 0.5 } };

    const static double hexpoints[8][3] =
      { 
        { 0, 0, 0 },
        { 1, 0, 0 },
        { 1, 1, 0 },
        { 0, 1, 0 },
        { 0, 0, 1 },
        { 1, 0, 1 },
        { 1, 1, 1 },
        { 0, 1, 1 }
      };
  
    int np, i;
    const double (*pp)[3];
    switch (GetType())
      {
      case TET:
        {
          np = 4;
          pp = tetpoints;
          break;
        }
      case PRISM:
      case PRISM12:
        {
          np = 6;
          pp = prismpoints;
          break;
        }
      case TET10:
        {
          np = 10;
          pp = tet10points;
          break;
        }
      case PYRAMID:
        {
          np = 5;
          pp = pyramidpoints;
          break;
        }
      case HEX:
        {
          np = 8;
          pp = hexpoints;
          break;
        }
      default:
        {
          cout << "GetNodesLocal not implemented for element " << GetType() << endl;
          np = 0;
        }
      }
  
    points.SetSize(0);
    for (i = 0; i < np; i++)
      points.Append (Point3d (pp[i][0], pp[i][1], pp[i][2]));
  }
#endif






  void Element :: GetNodesLocalNew (NgArray<Point<3> > & points) const
  {
    const static double tetpoints[4][3] =
      {      
        { 1, 0, 0 },
        { 0, 1, 0 },
        { 0, 0, 1 },
        { 0, 0, 0 }
      };
  
    const static double prismpoints[6][3] =
      {
        { 1, 0, 0 },
        { 0, 1, 0 },
        { 0, 0, 0 },
        { 1, 0, 1 },
        { 0, 1, 1 },
        { 0, 0, 1 }
      };
  
    const static double pyramidpoints[6][3] =
      { { 0, 0, 0 },
        { 1, 0, 0 },
        { 1, 1, 0 },
        { 0, 1, 0 },
        { 0, 0, 1 } };
  
    const static double tet10points[10][3] =
      { { 0, 0, 0 },
        { 1, 0, 0 },
        { 0, 1, 0 },
        { 0, 0, 1 },
        { 0.5, 0, 0 },
        { 0, 0.5, 0 },
        { 0, 0, 0.5 },
        { 0.5, 0.5, 0 },
        { 0.5, 0, 0.5 },
        { 0, 0.5, 0.5 } };

    const static double hexpoints[8][3] =
      { 
        { 0, 0, 0 },
        { 1, 0, 0 },
        { 1, 1, 0 },
        { 0, 1, 0 },
        { 0, 0, 1 },
        { 1, 0, 1 },
        { 1, 1, 1 },
        { 0, 1, 1 }
      };
  

  
    int np, i;
    const double (*pp)[3];
    switch (GetType())
      {
      case TET:
        {
          np = 4;
          pp = tetpoints;
          break;
        }
      case PRISM:
      case PRISM12:
        {
          np = 6;
          pp = prismpoints;
          break;
        }
      case TET10:
        {
          np = 10;
          pp = tet10points;
          break;
        }
      case PYRAMID:
        {
          np = 5;
          pp = pyramidpoints;
          break;
        }
      case HEX:
        {
          np = 8;
          pp = hexpoints;
          break;
        }
      default:
        {
          cout << "GetNodesLocal not implemented for element " << GetType() << endl;
          np = 0;
	  pp = NULL;
        }
      }
  
    points.SetSize(0);
    for (i = 0; i < np; i++)
      points.Append (Point<3> (pp[i][0], pp[i][1], pp[i][2]));
  }

















  void Element :: GetSurfaceTriangles (NgArray<Element2d> & surftrigs) const
  {
    static int tet4trigs[][3] = 
      { { 2, 3, 4 },
        { 3, 1, 4 },
        { 1, 2, 4 },
        { 2, 1, 3 } };

    static int tet10trigs[][3] = 
      { { 2, 8, 9 }, { 3, 10, 8}, { 4, 9, 10 }, { 9, 8, 10 },
        { 3, 6, 10 }, { 1, 7, 6 }, { 4, 10, 7 }, { 6, 7, 10 },
        { 1, 5, 7 }, { 2, 9, 5 }, { 4, 7, 9 }, { 5, 9, 7 },
        { 1, 6, 5 }, { 2, 5, 8 }, { 3, 8, 6 }, { 5, 6, 8 }
      };

    static int pyramidtrigs[][3] =
      {
        { 1, 3, 2 },
        { 1, 4, 3 },
        { 1, 2, 5 },
        { 2, 3, 5 },
        { 3, 4, 5 },
        { 4, 1, 5 }
      };

    static int prismtrigs[][3] =
      {
        { 1, 3, 2 },
        { 4, 5, 6 },
        { 1, 2, 4 },
        { 4, 2, 5 },
        { 2, 3, 5 },
        { 5, 3, 6 },
        { 3, 1, 6 },
        { 6, 1, 4 }
      };
  
    static int hextrigs[][3] = 
      {
        { 1, 3, 2 },
        { 1, 4, 3 }, 
        { 5, 6, 7 },
        { 5, 7, 8 },
        { 1, 2, 6 },
        { 1, 6, 5 },
        { 2, 3, 7 },
        { 2, 7, 6 },
        { 3, 4, 8 },
        { 3, 8, 7 },
        { 4, 1, 8 },
        { 1, 5, 8 }
      };

    int j;

    int nf;
    int (*fp)[3];

    switch (GetType())
      {
      case TET:
        {
          nf = 4;
          fp = tet4trigs;
          break;
        }
      case PYRAMID:
        {
          nf = 6;
          fp = pyramidtrigs;
          break;
        }
      case PRISM:
      case PRISM12:
        {
          nf = 8;
          fp = prismtrigs;
          break;
        }
      case TET10:
        {
          nf = 16;
          fp = tet10trigs;
          break;
        }
      case HEX:
        {
          nf = 12;
          fp = hextrigs;
          break;
        }
      default:
        {
          nf = 0;
          fp = NULL;
        }
      }

  
    surftrigs.SetSize (nf);
    for (j = 0; j < nf; j++)
      {
        surftrigs.Elem(j+1) = Element2d(TRIG);
        surftrigs.Elem(j+1).PNum(1) = fp[j][0];
        surftrigs.Elem(j+1).PNum(2) = fp[j][1];
        surftrigs.Elem(j+1).PNum(3) = fp[j][2];
      }
  }





  NgArray< shared_ptr < IntegrationPointData > > ipdtet;
  NgArray< shared_ptr < IntegrationPointData > > ipdtet10;



  int Element :: GetNIP () const
  {
    int nip;
    switch (typ)
      {
      case TET: nip = 1; break;
      case TET10: nip = 8; break;
      default: nip = 0; break;
      }
    return nip;
  }

  void Element :: 
  GetIntegrationPoint (int ip, Point<3> & p, double & weight) const
  {
    static double eltetqp[1][4] =
      {
        { 0.25, 0.25, 0.25, 1.0/6.0 }
      };

    static double eltet10qp[8][4] =
      {
        { 0.585410196624969, 0.138196601125011, 0.138196601125011, 1.0/24.0 },
        { 0.138196601125011, 0.585410196624969, 0.138196601125011, 1.0/24.0 },
        { 0.138196601125011, 0.138196601125011, 0.585410196624969, 1.0/24.0 },
        { 0.138196601125011, 0.138196601125011, 0.138196601125011, 1.0/24.0 },
        { 1, 0, 0, 1 },
        { 0, 1, 0, 1 },
        { 0, 0, 1, 1 },
        { 0, 0, 0, 1 },
      };
    
    double * pp = NULL;
    switch (typ)
      {
      case TET: pp = &eltetqp[0][0]; break;
      case TET10: pp = &eltet10qp[ip-1][0]; break;
      default:
        throw NgException ("illegal element shape in GetIntegrationPoint");
      }

    p(0) = pp[0];
    p(1) = pp[1];
    p(2) = pp[2];
    weight = pp[3];
  }

  void Element :: 
  GetTransformation (int ip, const T_POINTS & points,
                     DenseMatrix & trans) const
  {
    int np = GetNP();
    DenseMatrix pmat(3, np), dshape(3, np);
    pmat.SetSize (3, np);
    dshape.SetSize (3, np);

    Point<3> p;
    double w;

    GetPointMatrix (points, pmat);
    GetIntegrationPoint (ip, p, w);
    GetDShape (p, dshape);
  
    CalcABt (pmat, dshape, trans);

    /*
      (*testout) << "p = " << p  << endl
      << "pmat = " << pmat << endl
      << "dshape = " << dshape << endl
      << "tans = " << trans << endl;
    */
  }

  void Element :: 
  GetTransformation (int ip, class DenseMatrix & pmat,
                     class DenseMatrix & trans) const
  {
    int np = GetNP();

    if (pmat.Width() != np || pmat.Height() != 3)
      {
        (*testout) << "GetTransofrmation: pmat doesn't fit" << endl;
        return;
      }

    ComputeIntegrationPointData ();
    DenseMatrix * dshapep = 0;
    switch (GetType())
      {
      case TET: dshapep = &ipdtet.Get(ip)->dshape; break;
      case TET10: dshapep = &ipdtet10.Get(ip)->dshape; break;
      default:
        PrintSysError ("Element::GetTransformation, illegal type ", int(typ));
      }
  
    CalcABt (pmat, *dshapep, trans);
  }


  void Element :: GetShape (const Point<3> & hp, Vector & shape) const
  {
    if (shape.Size() != GetNP())
      {
        cerr << "Element::GetShape: Length not fitting" << endl;
        return;
      }

    switch (typ)
      {
      case TET:
        {
          shape(0) = 1 - hp[0] - hp[1] - hp[2]; 
          shape(1) = hp[0];
          shape(2) = hp[1];
          shape(3) = hp[2];
          break;
        }
      case TET10:
        {
          double lam1 = 1 - hp[0] - hp[1] - hp[2];
          double lam2 = hp[0];
          double lam3 = hp[1];
          double lam4 = hp[2];
	
          shape(4) = 4 * lam1 * lam2;
          shape(5) = 4 * lam1 * lam3;
          shape(6) = 4 * lam1 * lam4;
          shape(7) = 4 * lam2 * lam3;
          shape(8) = 4 * lam2 * lam4;
          shape(9) = 4 * lam3 * lam4;
	
          shape(0) = lam1 - 0.5 * (shape(4) + shape(5) + shape(6));
          shape(1) = lam2 - 0.5 * (shape(4) + shape(7) + shape(8));
          shape(2) = lam3 - 0.5 * (shape(5) + shape(7) + shape(9));
          shape(3) = lam4 - 0.5 * (shape(6) + shape(8) + shape(9));
          break;
        }

      case PRISM:
        {
          shape(0) = hp(0) * (1-hp(2));
          shape(1) = hp(1) * (1-hp(2));
          shape(2) = (1-hp(0)-hp(1)) * (1-hp(2));
          shape(3) = hp(0) * hp(2);
          shape(4) = hp(1) * hp(2);
          shape(5) = (1-hp(0)-hp(1)) * hp(2);
          break;
        }
      case HEX:
        {
          shape(0) = (1-hp(0))*(1-hp(1))*(1-hp(2));
          shape(1) = (  hp(0))*(1-hp(1))*(1-hp(2));
          shape(2) = (  hp(0))*(  hp(1))*(1-hp(2));
          shape(3) = (1-hp(0))*(  hp(1))*(1-hp(2));
          shape(4) = (1-hp(0))*(1-hp(1))*(  hp(2));
          shape(5) = (  hp(0))*(1-hp(1))*(  hp(2));
          shape(6) = (  hp(0))*(  hp(1))*(  hp(2));
          shape(7) = (1-hp(0))*(  hp(1))*(  hp(2));
          break;
        }
      default:
        throw NgException("Element :: GetShape not implemented for that element");
      }
  }


  template <typename T>
  void Element :: GetShapeNew (const Point<3,T> & p, TFlatVector<T> shape) const
  {
    /*
      if (shape.Size() < GetNP())
      {
      cerr << "Element::GetShape: Length not fitting" << endl;
      return;
      }
    */

    switch (typ)
      {
      case TET:
        {
          shape(0) = p(0);
          shape(1) = p(1);
          shape(2) = p(2);
          shape(3) = 1-p(0)-p(1)-p(2);
          break;
        }

      case TET10:
        {
          T lam1 = p(0);
          T lam2 = p(1);
          T lam3 = p(2);
          T lam4 = 1-p(0)-p(1)-p(2);
	
          shape(0) = 2 * lam1 * (lam1-0.5);
          shape(1) = 2 * lam2 * (lam2-0.5);
          shape(2) = 2 * lam3 * (lam3-0.5);
          shape(3) = 2 * lam4 * (lam4-0.5);

          shape(4) = 4 * lam1 * lam2;
          shape(5) = 4 * lam1 * lam3;
          shape(6) = 4 * lam1 * lam4;
          shape(7) = 4 * lam2 * lam3;
          shape(8) = 4 * lam2 * lam4;
          shape(9) = 4 * lam3 * lam4;
	
          break;
        }


      case PYRAMID:
        {
          T noz = 1-p(2);
          // if (noz == 0.0) noz = 1e-10;
          noz += T(1e-12);

          T xi  = p(0) / noz;
          T eta = p(1) / noz;
          shape(0) = (1-xi)*(1-eta) * (noz);
          shape(1) = (  xi)*(1-eta) * (noz);
          shape(2) = (  xi)*(  eta) * (noz);
          shape(3) = (1-xi)*(  eta) * (noz);
          shape(4) = p(2);
          break;
        }
      case PYRAMID13:
        {
	  T x = p(0);
	  T y = p(1);
	  T z = p(2);
          z *= 1-1e-12;
          shape[0] = (-z + z*(2*x + z - 1)*(2*y + z - 1)/(-z + 1) + (-2*x - z + 2)*(-2*y - z + 2))*(-0.5*x - 0.5*y - 0.5*z + 0.25);
          shape[1] = (0.5*x - 0.5*y - 0.25)*(-z - z*(2*x + z - 1)*(2*y + z - 1)/(-z + 1) + (2*x + z)*(-2*y - z + 2));
          shape[2] = (-z + z*(2*x + z - 1)*(2*y + z - 1)/(-z + 1) + (2*x + z)*(2*y + z))*(0.5*x + 0.5*y + 0.5*z - 0.75);
          shape[3] = (-0.5*x + 0.5*y - 0.25)*(-z - z*(2*x + z - 1)*(2*y + z - 1)/(-z + 1) + (2*y + z)*(-2*x - z + 2));
          shape[4] = z*(2*z - 1);
          shape[5] = 2*x*(-2*x - 2*z + 2)*(-2*y - 2*z + 2)/(-2*z + 2);
          shape[6] = 4*x*y*(-2*x - 2*z + 2)/(-2*z + 2);
          shape[7] = 2*y*(-2*x - 2*z + 2)*(-2*y - 2*z + 2)/(-2*z + 2);
          shape[8] = 4*x*y*(-2*y - 2*z + 2)/(-2*z + 2);
          shape[9] = z*(-2*x - 2*z + 2)*(-2*y - 2*z + 2)/(-z + 1);
          shape[10] = 2*x*z*(-2*y - 2*z + 2)/(-z + 1);
          shape[11] = 4*x*y*z/(-z + 1);
          shape[12] = 2*y*z*(-2*x - 2*z + 2)/(-z + 1);
          break;
        }
      case PRISM:
        {
          shape(0) = p(0) * (1-p(2));
          shape(1) = p(1) * (1-p(2));
          shape(2) = (1-p(0)-p(1)) * (1-p(2));
          shape(3) = p(0) * p(2);
          shape(4) = p(1) * p(2);
          shape(5) = (1-p(0)-p(1)) * p(2);
          break;
        }
      case PRISM15:
        {
	  T x = p(0);
	  T y = p(1);
	  T z = p(2);
          T lam = 1-x-y;
          T lamz = 1-z;
          shape[0] = (2*x*x-x) * (2*lamz*lamz-lamz);
          shape[1] = (2*y*y-y) * (2*lamz*lamz-lamz);
          shape[2] = (2*lam*lam-lam) * (2*lamz*lamz-lamz);
          shape[3] = (2*x*x-x) * (2*z*z-z);
          shape[4] = (2*y*y-y) * (2*z*z-z);
          shape[5] = (2*lam*lam-lam) * (2*z*z-z);
          shape[6] = 4 * x * y * (2*lamz*lamz-lamz);
          shape[7] = 4 * x * lam * (2*lamz*lamz-lamz);
          shape[8] = 4 * y * lam * (2*lamz*lamz-lamz);
          shape[9] = x * 4 * z * (1-z);
          shape[10] = y * 4 * z * (1-z);
          shape[11] = lam * 4 * z * (1-z);
          shape[12] = 4 * x * y * (2*z*z-z);
          shape[13] = 4 * x * lam * (2*z*z-z);
          shape[14] = 4 * y * lam * (2*z*z-z);
          break;
        }
      case HEX:
        {
          shape(0) = (1-p(0))*(1-p(1))*(1-p(2));
          shape(1) = (  p(0))*(1-p(1))*(1-p(2));
          shape(2) = (  p(0))*(  p(1))*(1-p(2));
          shape(3) = (1-p(0))*(  p(1))*(1-p(2));
          shape(4) = (1-p(0))*(1-p(1))*(  p(2));
          shape(5) = (  p(0))*(1-p(1))*(  p(2));
          shape(6) = (  p(0))*(  p(1))*(  p(2));
          shape(7) = (1-p(0))*(  p(1))*(  p(2));
          break;
        }
      case HEX20:
	{
	  T x = p(0);
	  T y = p(1);
	  T z = p(2);
	  shape[0] = (1-x)*(1-y)*(1-z);
	  shape[1] =    x *(1-y)*(1-z);
	  shape[2] =    x *   y *(1-z);
	  shape[3] = (1-x)*   y *(1-z);
	  shape[4] = (1-x)*(1-y)*(z);
	  shape[5] =    x *(1-y)*(z);
	  shape[6] =    x *   y *(z);
	  shape[7] = (1-x)*   y *(z);

          T sigma[8]={(1-x)+(1-y)+(1-z),x+(1-y)+(1-z),x+y+(1-z),(1-x)+y+(1-z),
                      (1-x)+(1-y)+z,x+(1-y)+z,x+y+z,(1-x)+y+z};

          static const int e[12][2] =
            {
              { 0, 1 }, { 2, 3 }, { 3, 0 }, { 1, 2 },
              { 4, 5 }, { 6, 7 }, { 7, 4 }, { 5, 6 },
              { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 },
            };

          for (int i = 0; i < 12; i++)
            {
              T lame = shape[e[i][0]]+shape[e[i][1]];
              T xi = sigma[e[i][1]]-sigma[e[i][0]];
              shape[8+i] = (1-xi*xi)*lame;
            }
          for (int i = 0; i < 12; i++)
            {
              shape[e[i][0]] -= 0.5 * shape[8+i];
              shape[e[i][1]] -= 0.5 * shape[8+i];
            }
          break;
	}
      default:
        throw NgException("Element :: GetNewShape not implemented for that element");
      }
  }



  void Element :: 
  GetDShape (const Point<3> & hp, DenseMatrix & dshape) const
  {
    int np = GetNP();
    if (dshape.Height() != 3 || dshape.Width() != np)
      {
        cerr << "Element::DShape: Sizes don't fit" << endl;
        return;
      }

    double eps = 1e-6;
    Vector shaper(np), shapel(np);

    for (auto i : Range(3))
      {
        Point<3> pr(hp), pl(hp);
        pr[i] += eps;
        pl[i] -= eps;
      
        GetShape (pr, shaper);
        GetShape (pl, shapel);
        for (int j = 0; j < np; j++)
          dshape(i, j) = (shaper(j) - shapel(j)) / (2 * eps);
      }
  }

  template <typename T>
  void Element :: 
  GetDShapeNew (const Point<3,T> & p, MatrixFixWidth<3,T> & dshape) const
  {
    switch (typ)
      {
      case TET:
        {
          dshape = T(0.0);
          dshape(0,0) = 1;
          dshape(1,1) = 1;
          dshape(2,2) = 1;
          dshape(3,0) = -1;
          dshape(3,1) = -1;
          dshape(3,2) = -1;
          break;
        }
      case PRISM:
        {
          dshape = T(0.0);
          dshape(0,0) = 1-p(2);
          dshape(0,2) = -p(0);
          dshape(1,1) = 1-p(2);
          dshape(1,2) = -p(1);
          dshape(2,0) = -(1-p(2));
          dshape(2,1) = -(1-p(2));
          dshape(2,2) = -(1-p(0)-p(1));

          dshape(3,0) = p(2);
          dshape(3,2) = p(0);
          dshape(4,1) = p(2);
          dshape(4,2) = p(1);
          dshape(5,0) = -p(2);
          dshape(5,1) = -p(2);
          dshape(5,2) = 1-p(0)-p(1);
          break;
        }

      default:
        {
          int np = GetNP();
          double eps = 1e-6;
          NgArrayMem<T,100> mem(2*np);
          TFlatVector<T> shaper(np, &mem[0]);
          TFlatVector<T> shapel(np, &mem[np]);
          // Vector shaper(np), shapel(np);
	
          for (int i = 0; i < 3; i++)
            {
              Point<3,T> pr(p), pl(p);
              pr(i) += eps;
              pl(i) -= eps;
	    
              GetShapeNew (pr, shaper);
              GetShapeNew (pl, shapel);
              for (int j = 0; j < np; j++)
                dshape(j, i) = (shaper(j) - shapel(j)) / (2 * eps);
            }
        }
      }
  }

  template void Element2d :: GetShapeNew (const Point<2,double> & p, TFlatVector<double> shape) const;
  template void Element2d :: GetShapeNew (const Point<2,SIMD<double>> & p, TFlatVector<SIMD<double>> shape) const;

  template void Element2d::GetDShapeNew<double> (const Point<2> &, MatrixFixWidth<2> &) const;
  template void Element2d::GetDShapeNew<SIMD<double>> (const Point<2,SIMD<double>> &, MatrixFixWidth<2,SIMD<double>> &) const;


  template DLL_HEADER void Element :: GetShapeNew (const Point<3,double> & p, TFlatVector<double> shape) const;
  template DLL_HEADER void Element :: GetShapeNew (const Point<3,SIMD<double>> & p, TFlatVector<SIMD<double>> shape) const;
  
  template void Element::GetDShapeNew<double> (const Point<3> &, MatrixFixWidth<3> &) const;
  template void Element::GetDShapeNew<SIMD<double>> (const Point<3,SIMD<double>> &, MatrixFixWidth<3,SIMD<double>> &) const;


  void Element :: 
  GetPointMatrix (const T_POINTS & points,
                  DenseMatrix & pmat) const
  {
    int np = GetNP();
    for (int i = 1; i <= np; i++)
      {
        const auto& p = points[PNum(i)];
        pmat.Elem(1, i) = p[0];
        pmat.Elem(2, i) = p[1];
        pmat.Elem(3, i) = p[2];
      }
  }




  double Element :: CalcJacobianBadness (const T_POINTS & points) const
  {
    int nip = GetNIP();
    DenseMatrix trans(3,3);
    DenseMatrix pmat;
  
    pmat.SetSize (3, GetNP());
    GetPointMatrix (points, pmat);

    double err = 0;
    for (int i = 1; i <= nip; i++)
      {
        GetTransformation (i, pmat, trans);

        // Frobenius norm
        double frob = 0;
        for (int j = 1; j <= 9; j++)
          frob += sqr (trans.Get(j));
        frob = sqrt (frob);
        frob /= 3;

        double det = -trans.Det();
      
        if (det <= 0)
          err += 1e12;
        else
          err += frob * frob * frob / det;
      }

    err /= nip;
    return err;
  }

  double Element :: 
  CalcJacobianBadnessDirDeriv (const T_POINTS & points,
                               int pi, Vec<3> & dir, double & dd) const
  {
    int i, j, k;
    int nip = GetNIP();
    DenseMatrix trans(3,3), dtrans(3,3), hmat(3,3);
    DenseMatrix pmat, vmat;
  
    pmat.SetSize (3, GetNP());
    vmat.SetSize (3, GetNP());

    GetPointMatrix (points, pmat);
  
    for (i = 1; i <= np; i++)
      for (j = 1; j <= 3; j++)
        vmat.Elem(j, i) = 0;
    for (j = 1; j <= 3; j++)
      vmat.Elem(j, pi) = dir(j-1);



    double err = 0;
    dd = 0;

    for (i = 1; i <= nip; i++)
      {
        GetTransformation (i, pmat, trans);
        GetTransformation (i, vmat, dtrans);


        // Frobenius norm
        double frob = 0;
        for (j = 1; j <= 9; j++)
          frob += sqr (trans.Get(j));
        frob = sqrt (frob);
      
        double dfrob = 0;
        for (j = 1; j <= 9; j++)
          dfrob += trans.Get(j) * dtrans.Get(j);
        dfrob = dfrob / frob;
      
        frob /= 3;      
        dfrob /= 3;

      
        double det = trans.Det();
        double ddet = 0;
      
        for (j = 1; j <= 3; j++)
          {
            hmat = trans;
            for (k = 1; k <= 3; k++)
              hmat.Elem(k, j) = dtrans.Get(k, j);
            ddet += hmat.Det();
          }


        det *= -1;
        ddet *= -1;

      
        if (det <= 0)
          err += 1e12;
        else
          {
            err += frob * frob * frob / det;
            dd += (3 * frob * frob * dfrob * det - frob * frob * frob * ddet) / (det * det);
          }
      }

    err /= nip;
    dd /= nip;
    return err;
  }

  double Element :: 
  CalcJacobianBadnessGradient (const T_POINTS & points,
                               int pi, Vec<3> & grad) const
  {
    int nip = GetNIP();
    DenseMatrix trans(3,3), dtrans(3,3), hmat(3,3);
    DenseMatrix pmat, vmat;
  
    pmat.SetSize (3, GetNP());
    vmat.SetSize (3, GetNP());

    GetPointMatrix (points, pmat);
  
    for (int i = 1; i <= np; i++)
      for (int j = 1; j <= 3; j++)
        vmat.Elem(j, i) = 0;
    for (int j = 1; j <= 3; j++)
      vmat.Elem(j, pi) = 1.;


    double err = 0;

    double dfrob[3];

    grad = 0;

    for (int i = 1; i <= nip; i++)
      {
        GetTransformation (i, pmat, trans);
        GetTransformation (i, vmat, dtrans);
 
        // Frobenius norm
        double frob = 0;
        for (int j = 1; j <= 9; j++)
          frob += sqr (trans.Get(j));
        frob = sqrt (frob);

        for(int k = 0; k<3; k++)
          {
            dfrob[k] = 0;
            for (int j = 1; j <= 3; j++)
              dfrob[k] += trans.Get(k+1,j) * dtrans.Get(k+1,j);
            dfrob[k] = dfrob[k] / (3.*frob);
          }

        frob /= 3;      

        double det = trans.Det();
        double ddet[3]; // = 0;
      
        for(int k=1; k<=3; k++)
          {
            int km1 = (k > 1) ? (k-1) : 3;
            int kp1 = (k < 3) ? (k+1) : 1;
            ddet[k-1] = 0;
            for(int j=1; j<=3; j++)
              {
                int jm1 = (j > 1) ? (j-1) : 3;
                int jp1 = (j < 3) ? (j+1) : 1;
	      
                ddet[k-1] += (-1.)* dtrans.Get(k,j) * ( trans.Get(km1,jm1)*trans.Get(kp1,jp1) - 
                                                        trans.Get(km1,jp1)*trans.Get(kp1,jm1) );
              }
          }

      
        det *= -1;
      
        if (det <= 0)
          err += 1e12;
        else
          {
            err += frob * frob * frob / det;
            double fac = (frob * frob)/(det * det);
            for(int j=0; j<3; j++)
              grad(j) += fac * (3 * dfrob[j] * det - frob * ddet[j]);
          }
      }

    err /= nip;
    grad *= 1./nip;
    return err;
  }





  void Element :: ComputeIntegrationPointData () const
  {
    switch (GetType())
      {
      case TET: if (ipdtet.Size()) return; break;
      case TET10: if (ipdtet10.Size()) return; break;
      default:
        PrintSysError ("Element::ComputeIntegrationPoint, illegal type ", int(typ));
      }

    switch (GetType())
      {
      case TET: ipdtet.SetSize(GetNIP()); break;
      case TET10: ipdtet10.SetSize(GetNIP()); break;
      default:
        PrintSysError ("Element::ComputeIntegrationPoint, illegal type2 ", int(typ));
      }


    for (int i = 1; i <= GetNIP(); i++)
      {
        IntegrationPointData * ipd = new IntegrationPointData;
        GetIntegrationPoint (i, ipd->p, ipd->weight);
        ipd->shape.SetSize(GetNP());
        ipd->dshape.SetSize(3, GetNP());

        GetShape (ipd->p, ipd->shape);
        GetDShape (ipd->p, ipd->dshape);

        switch (GetType())
          {
          case TET: ipdtet.Elem(i).reset(ipd); break;
          case TET10: ipdtet10.Elem(i).reset(ipd); break;
          default:
            PrintSysError ("Element::ComputeIntegrationPoint(2), illegal type ", int(typ));
          }
      }
  }







  FaceDescriptor ::  FaceDescriptor()
  { 
    surfnr = domin = domout  = bcprop = 0; 
    domin_singular = domout_singular = 0.;
    // Philippose - 06/07/2009
    // Initialise surface colour
    surfcolour = Vec<4>(0.0,1.0,0.0,1.0);
    tlosurf = -1; 
    // bcname = 0;
    firstelement = -1;
  }

  FaceDescriptor ::  FaceDescriptor(const FaceDescriptor& other)
    : surfnr(other.surfnr), domin(other.domin), domout(other.domout),
      tlosurf(other.tlosurf), bcprop(other.bcprop), 
      surfcolour(other.surfcolour), bcname(other.bcname),
      domin_singular(other.domin_singular), domout_singular(other.domout_singular)
  { 
    firstelement = -1;
  }

  FaceDescriptor :: 
  FaceDescriptor(int surfnri, int domini, int domouti, int tlosurfi)
  { 
    surfnr = surfnri; 
    domin = domini; 
    domout = domouti;
    // Philippose - 06/07/2009
    // Initialise surface colour
    surfcolour = Vec<4>(0.0,1.0,0.0,1.0);
    tlosurf = tlosurfi; 
    bcprop = surfnri;
    domin_singular = domout_singular = 0.;
    // bcname = 0;
    firstelement = -1;
  }

  FaceDescriptor :: FaceDescriptor(const Segment & seg)
  { 
    surfnr = seg.si; 
    domin = seg.domin+1;
    domout = seg.domout+1;
    // Philippose - 06/07/2009
    // Initialise surface colour
    surfcolour = Vec<4>(0.0,1.0,0.0,1.0);
    tlosurf = seg.tlosurf+1;
    bcprop = 0;
    domin_singular = domout_singular = 0.;
    // bcname = 0;
    firstelement = -1;
  }

  int FaceDescriptor ::  SegmentFits (const Segment & seg)
  {
    return
      surfnr == seg.si &&
      domin == seg.domin+1 &&
      domout == seg.domout+1  &&
      tlosurf == seg.tlosurf+1;
  }

  // string FaceDescriptor :: default_bcname = "default";
  /*
  const string & FaceDescriptor :: GetBCName () const
  {
    static string defaultstring = "default";
    if (bcname) return *bcname;
    return defaultstring;
  }
  */

  void FaceDescriptor :: SetBCName (string * bcn)
  {
    if (bcn)
      bcname = *bcn;
    else
      bcname = "default";
  }
  
  void FaceDescriptor :: DoArchive (Archive & ar)
  {
    ar & surfnr & domin & domout & tlosurf & bcprop
      & surfcolour & bcname   
      & domin_singular & domout_singular ;
      // don't need:  firstelement
  }
  

  
  ostream & operator<<(ostream  & s, const FaceDescriptor & fd)
  {
    s << "surfnr = " << fd.SurfNr() 
      << ", domin = " << fd.DomainIn()
      << ", domout = " << fd.DomainOut()
      << ", tlosurf = " << fd.TLOSurface()
      << ", bcprop = " << fd.BCProperty()
      << ", bcname = " << fd.GetBCName()
      << ", domin_sing = " << fd.DomainInSingular()
      << ", domout_sing = " << fd.DomainOutSingular()
      << ", colour = " << fd.SurfColour();
    return s;
  }


  Identifications :: Identifications (Mesh & amesh)
    : mesh(amesh), identifiedpoints(100), identifiedpoints_nr(100)
  {
    // identifiedpoints = new INDEX_2_HASHTABLE<int>(100);
    // identifiedpoints_nr = new INDEX_3_HASHTABLE<int>(100);
    maxidentnr = 0;
  }

  Identifications :: ~Identifications ()
  {
    ;
    // delete identifiedpoints;
    // delete identifiedpoints_nr;
  }

  void Identifications :: Delete ()
  {
    identifiedpoints.DeleteData();
    identifiedpoints_nr.DeleteData();

    /*
    delete identifiedpoints;
    identifiedpoints = new INDEX_2_HASHTABLE<int>(100);
    delete identifiedpoints_nr;
    identifiedpoints_nr = new INDEX_3_HASHTABLE<int>(100);
    */
    maxidentnr = 0;
  }

    void Identifications :: DoArchive (Archive & ar)
    {
      ar & maxidentnr;
      ar & identifiedpoints & identifiedpoints_nr;

      ar & idpoints_table;
      if (ar.Output())
        {
          size_t s = type.Size();
          ar & s;
          for (auto & t : type)
            ar & (unsigned char&)(t);
        }
      else
        {
          size_t s;
          ar & s;
          type.SetSize(s);
          for (auto & t : type)
            ar & (unsigned char&)(t);
        }
    }    

  

  void Identifications :: Add (PointIndex pi1, PointIndex pi2, int identnr)
  {
    //  (*testout) << "Identification::Add, pi1 = " << pi1 << ", pi2 = " << pi2 << ", identnr = " << identnr << endl;
    INDEX_2 pair (pi1, pi2);
    identifiedpoints.Set (pair, identnr);

    INDEX_3 tripl (pi1, pi2, identnr);
    identifiedpoints_nr.Set (tripl, 1);

    if (identnr > maxidentnr) maxidentnr = identnr;
    names.SetSize(maxidentnr);

    if (identnr+1 > idpoints_table.Size())
      idpoints_table.ChangeSize (identnr+1);
    idpoints_table.Add (identnr, pair);
  
    //  timestamp = NextTimeStamp();
  }

  int Identifications :: Get (PointIndex pi1, PointIndex pi2) const
  {
    INDEX_2 pair(pi1, pi2);
    if (identifiedpoints.Used (pair))
      return identifiedpoints.Get(pair);
    else
      return 0;
  }

  bool Identifications :: Get (PointIndex pi1, PointIndex pi2, int nr) const
  {
    INDEX_3 tripl(pi1, pi2, nr);
    if (identifiedpoints_nr.Used (tripl))
      return 1;
    else
      return 0;
  }



  int Identifications :: GetSymmetric (PointIndex pi1, PointIndex pi2) const
  {
    INDEX_2 pair(pi1, pi2);
    if (identifiedpoints.Used (pair))
      return identifiedpoints.Get(pair);

    pair = INDEX_2 (pi2, pi1);
    if (identifiedpoints.Used (pair))
      return identifiedpoints.Get(pair);

    return 0;
  }


  void Identifications :: GetMap (int identnr, NgArray<int,PointIndex::BASE> & identmap, bool symmetric) const
  {
    identmap.SetSize (mesh.GetNP());
    identmap = 0;

    if (identnr)
      for (int i = 0; i < idpoints_table[identnr].Size(); i++)
        {
          INDEX_2 pair = idpoints_table[identnr][i];
          identmap[pair.I1()] = pair.I2();
          if(symmetric)
            identmap[pair.I2()] = pair.I1();
        }

    else
      {
        cout << "getmap, identnr = " << identnr << endl;

        for (int i = 1; i <= identifiedpoints_nr.GetNBags(); i++)
          for (int j = 1; j <= identifiedpoints_nr.GetBagSize(i); j++)
            {
              INDEX_3 i3;
              int dummy;
              identifiedpoints_nr.GetData (i, j, i3, dummy);
	    
              if (i3.I3() == identnr || !identnr)
                {
                  identmap.Elem(i3.I1()) = i3.I2();
                  if(symmetric)
                    identmap.Elem(i3.I2()) = i3.I1();
                }
            }  
      }

  }


  void Identifications :: GetPairs (int identnr, 
                                    NgArray<INDEX_2> & identpairs) const
  {
    identpairs.SetSize(0);
  
    if (identnr == 0)
      for (int i = 1; i <= identifiedpoints.GetNBags(); i++)
        for (int j = 1; j <= identifiedpoints.GetBagSize(i); j++)
          {
            INDEX_2 i2;
            int nr;
            identifiedpoints.GetData (i, j, i2, nr);
            identpairs.Append (i2);
          }  
    else
      for (int i = 1; i <= identifiedpoints_nr.GetNBags(); i++)
        for (int j = 1; j <= identifiedpoints_nr.GetBagSize(i); j++)
          {
            INDEX_3 i3;
            int dummy;
            identifiedpoints_nr.GetData (i, j, i3 , dummy);
	  
            if (i3.I3() == identnr)
              identpairs.Append (INDEX_2(i3.I1(), i3.I2()));
          }  
  }


  void Identifications :: SetMaxPointNr (int maxpnum)
  {
    for (int i = 1; i <= identifiedpoints.GetNBags(); i++)
      for (int j = 1; j <= identifiedpoints.GetBagSize(i); j++)
        {
          INDEX_2 i2;
          int nr;
          identifiedpoints.GetData (i, j, i2, nr);
	
          if (i2.I1() > maxpnum || i2.I2() > maxpnum)
            {
              i2.I1() = i2.I2() = -1;
              identifiedpoints.SetData (i, j, i2, -1);	    
            }
        }
  }


  void Identifications :: Print (ostream & ost) const
  {
    ost << "Identifications:" << endl;
    ost << "pairs: " << endl << identifiedpoints << endl;
    ost << "pairs and nr: " << endl << identifiedpoints_nr << endl;
    ost << "table: " << endl << idpoints_table << endl;
  }


  MeshingParameters :: MeshingParameters ()
  {
    // optimize3d = "cmdmustm";
    //optimize3d = "cmdmstm";
    // optsteps3d = 3;
    // optimize2d = "smsmsmSmSmSm";
    // optsteps2d = 3;
    // opterrpow = 2;
    // blockfill = 1;
    // filldist = 0.1;
    // safety = 5;
    // relinnersafety = 3;
    // uselocalh = 1;
    // grading = 0.3;
    // delaunay = 1;
    // maxh = 1e10;
    // minh = 0;
    // meshsizefilename = NULL;
    // startinsurface = 0;
    // checkoverlap = 1;
    // checkoverlappingboundary = 1;
    // checkchartboundary = 1;
    // curvaturesafety = 2;
    // segmentsperedge = 1;
    // parthread = 0;

    // elsizeweight = 0.2;
    // giveuptol2d = 200;
    // giveuptol = 10;
    // maxoutersteps = 10;
    // starshapeclass = 5;
    // baseelnp = 0;
    // sloppy = 1;

    // badellimit = 175;
    // check_impossible = 0;
    // secondorder = 0;
  }

  void MeshingParameters :: Print (ostream & ost) const
  {
    ost << "Meshing parameters: " << endl
        << "optimize3d = " << optimize3d << endl
        << "optsteps3d = " << optsteps3d << endl
        << " optimize2d = " <<  optimize2d << endl
        << " optsteps2d = " <<  optsteps2d << endl
        << " opterrpow = " <<  opterrpow << endl
        << " blockfill = " <<  blockfill << endl
        << " filldist = " <<  filldist << endl
        << " safety = " <<  safety << endl
        << " relinnersafety = " <<  relinnersafety << endl
        << " uselocalh = " <<  uselocalh << endl
        << " grading = " <<  grading << endl
        << " delaunay = " <<  delaunay << endl
        << " maxh = " <<  maxh << endl
        << " meshsizefilename = " <<  meshsizefilename << endl
        << " startinsurface = " <<  startinsurface << endl
        << " checkoverlap = " <<  checkoverlap << endl
        << " checkchartboundary = " <<  checkchartboundary << endl
        << " curvaturesafety = " <<  curvaturesafety << endl
        << " segmentsperedge = " <<  segmentsperedge << endl
        << " parthread = " <<  parthread << endl
        << " elsizeweight = " <<  elsizeweight << endl
        << " giveuptol2d = " <<  giveuptol2d << endl
        << " giveuptol = " <<  giveuptol << endl
        << " maxoutersteps = " <<  maxoutersteps << endl
        << " starshapeclass = " <<  starshapeclass << endl
        << " baseelnp        = " <<  baseelnp        << endl
        << " sloppy = " <<  sloppy << endl
        << " badellimit = " <<  badellimit << endl
        << " secondorder = " <<  secondorder << endl
        << " elementorder = " <<  elementorder << endl
        << " quad = " <<  quad << endl
        << " inverttets = " <<  inverttets << endl
        << " inverttrigs = " <<  inverttrigs << endl
        << "closeedge enabled = " << closeedgefac.has_value() << endl
        << "closeedgefac = " << closeedgefac.value_or(0.) << endl;
  }

  /*
  void MeshingParameters :: CopyFrom(const MeshingParameters & other)
  {
    //strcpy(optimize3d,other.optimize3d); 
    optimize3d = other.optimize3d;
    optsteps3d = other.optsteps3d;
    //strcpy(optimize2d,other.optimize2d); 
    optimize2d = other.optimize2d;
    optsteps2d = other.optsteps2d;
    opterrpow = other.opterrpow;
    blockfill = other.blockfill;
    filldist = other.filldist;
    safety = other.safety;
    relinnersafety = other.relinnersafety;
    uselocalh = other.uselocalh;
    grading = other.grading;
    delaunay = other.delaunay;
    maxh = other.maxh;
    //strcpy(const_cast<char*>(meshsizefilename), other.meshsizefilename);
    //const_cast<char*>(meshsizefilename) = other.meshsizefilename; //???
    meshsizefilename = other.meshsizefilename;
    startinsurface = other.startinsurface;
    checkoverlap = other.checkoverlap;
    checkoverlappingboundary = other.checkoverlappingboundary;
    checkchartboundary = other.checkchartboundary;
    curvaturesafety = other.curvaturesafety;
    segmentsperedge = other.segmentsperedge;
    parthread = other.parthread;
    elsizeweight = other.elsizeweight;
    giveuptol2d = other.giveuptol2d;
    giveuptol = other.giveuptol;
    maxoutersteps = other.maxoutersteps;
    starshapeclass = other.starshapeclass;
    baseelnp = other.baseelnp;       
    sloppy = other.sloppy;
    badellimit = other.badellimit;
    secondorder = other.secondorder;
    elementorder = other.elementorder;
    quad = other.quad;
    inverttets = other.inverttets;
    inverttrigs = other.inverttrigs;
  }
  */

  DebugParameters :: DebugParameters ()
  {
    slowchecks = 0;
    haltsuccess = 0;
    haltnosuccess = 0;
    haltlargequalclass = 0;
    haltsegment = 0;
    haltsegmentp1 = 0;
    haltsegmentp2 = 0;
  };



}

