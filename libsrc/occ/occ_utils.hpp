#ifndef FILE_OCC_UTILS_INCLUDED
#define FILE_OCC_UTILS_INCLUDED

#include <BRepGProp.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <Standard_Version.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Vertex.hxx>
#include <gp_Trsf.hxx>

#include "meshing.hpp"

#if OCC_VERSION_MAJOR>=7 && OCC_VERSION_MINOR>=4
#define OCC_HAVE_DUMP_JSON
#endif

namespace netgen
{
    typedef Handle(TopoDS_TShape) T_Shape;

    inline Point<3> occ2ng (const gp_Pnt & p)
    {
        return Point<3> (p.X(), p.Y(), p.Z());
    }

    inline Point<2> occ2ng (const gp_Pnt2d & p)
    {
        return Point<2> (p.X(), p.Y());
    }

    inline Vec<3> occ2ng (const gp_Vec & v)
    {
        return Vec<3> (v.X(), v.Y(), v.Z());
    }

    DLL_HEADER Point<3> occ2ng (T_Shape shape);

    inline Point<3> occ2ng (const TopoDS_Shape & s)
    {
        return occ2ng(s.TShape());
    }

    inline Point<3> occ2ng (const TopoDS_Vertex & v)
    {
        return occ2ng (BRep_Tool::Pnt (v));
    }

    DLL_HEADER Transformation<3> occ2ng (const gp_Trsf & t);

    inline gp_Pnt ng2occ (const Point<3> & p)
    {
        return gp_Pnt(p(0), p(1), p(2));
    }

    DLL_HEADER Box<3> GetBoundingBox( const TopoDS_Shape & shape );

    class OCCIdentification
    {
    public:
      T_Shape from;
      T_Shape to;
      Transformation<3> trafo;
      string name;
      Identifications::ID_TYPE type;
      bool opposite_direction;
    };


    class MyExplorer
    {
      class Iterator
      {
        TopExp_Explorer exp;
      public:
        Iterator (TopoDS_Shape ashape, TopAbs_ShapeEnum atoFind, TopAbs_ShapeEnum atoAvoid)
          : exp(ashape, atoFind, atoAvoid) { }
        auto operator*() { return exp.Current(); }
        Iterator & operator++() { exp.Next(); return *this; }
        bool operator!= (nullptr_t nu) { return exp.More(); }
      };

    public:
      TopoDS_Shape shape;
      TopAbs_ShapeEnum toFind;
      TopAbs_ShapeEnum toAvoid;
      MyExplorer (TopoDS_Shape ashape, TopAbs_ShapeEnum atoFind, TopAbs_ShapeEnum atoAvoid = TopAbs_SHAPE)
        : shape(ashape), toFind(atoFind), toAvoid(atoAvoid) { ; }
      Iterator begin() { return Iterator(shape, toFind, toAvoid); }
      auto end() { return nullptr; }
    };

    inline auto Explore (TopoDS_Shape shape, TopAbs_ShapeEnum toFind, TopAbs_ShapeEnum toAvoid = TopAbs_SHAPE)
    {
      return MyExplorer (shape, toFind, toAvoid);
    }


    class IndexMapIterator
    {
      class Iterator
      {
        const TopTools_IndexedMapOfShape & indmap;
        int i;
      public:
        Iterator (const TopTools_IndexedMapOfShape & aindmap, int ai)
          : indmap(aindmap), i(ai) { ; }
        auto operator*() { return tuple(i, indmap(i)); }
        Iterator & operator++() { i++; return *this; }
        bool operator!= (const Iterator & i2) { return i != i2.i; }
      };

    public:
      const TopTools_IndexedMapOfShape & indmap;
      IndexMapIterator (const TopTools_IndexedMapOfShape & aindmap) : indmap(aindmap) { }
      Iterator begin() { return Iterator(indmap, 1); }
      Iterator end() { return Iterator(indmap, indmap.Extent()+1); }
    };

    inline auto Enumerate (const TopTools_IndexedMapOfShape & indmap)
    {
      return IndexMapIterator(indmap);
    }

    struct ShapeLess
    {
      bool operator() (const TopoDS_Shape& s1, const TopoDS_Shape& s2) const
      {
        return s1.TShape() < s2.TShape();
      }
    };

    class ListOfShapes : public std::vector<TopoDS_Shape>
    {
    public:
      DLL_HEADER TopoDS_Shape Max(gp_Vec dir);
      DLL_HEADER TopoDS_Shape Nearest(gp_Pnt pnt);
      DLL_HEADER ListOfShapes SubShapes(TopAbs_ShapeEnum type) const;

      ListOfShapes Solids() const
      {
        return SubShapes(TopAbs_SOLID);
      }
      ListOfShapes Faces() const
      {
        return SubShapes(TopAbs_FACE);
      }
      ListOfShapes Wires() const
      {
        return SubShapes(TopAbs_WIRE);
      }
      ListOfShapes Edges() const
      {
        return SubShapes(TopAbs_EDGE);
      }
      ListOfShapes Vertices() const
      {
        return SubShapes(TopAbs_VERTEX);
      }

      ListOfShapes operator*(const ListOfShapes& other) const
      {
        ListOfShapes common;
        for(const auto& shape : (*this))
          for(const auto& shape_o : other)
            if(shape.IsSame(shape_o))
              common.push_back(shape);
        return common;
      }
    };

    inline ListOfShapes GetSolids(const TopoDS_Shape & shape)
    {
      ListOfShapes sub;
      for (TopExp_Explorer e(shape, TopAbs_SOLID); e.More(); e.Next())
        sub.push_back(e.Current());
      return sub;
    }

    inline ListOfShapes GetFaces(const TopoDS_Shape & shape)
    {
      ListOfShapes sub;
      for (TopExp_Explorer e(shape, TopAbs_FACE); e.More(); e.Next())
        sub.push_back(e.Current());
      return sub;
    }

    inline ListOfShapes GetWires(const TopoDS_Shape & shape)
    {
      ListOfShapes sub;
      for (TopExp_Explorer e(shape, TopAbs_WIRE); e.More(); e.Next())
        sub.push_back(e.Current());
      return sub;
    }

    inline ListOfShapes GetEdges(const TopoDS_Shape & shape)
    {
      ListOfShapes sub;
      for (TopExp_Explorer e(shape, TopAbs_EDGE); e.More(); e.Next())
        sub.push_back(e.Current());
      return sub;
    }

    inline ListOfShapes GetVertices(const TopoDS_Shape & shape)
    {
      ListOfShapes sub;
      for (TopExp_Explorer e(shape, TopAbs_VERTEX); e.More(); e.Next())
        sub.push_back(e.Current());
      return sub;
    }

  class DirectionalInterval
  {
  public:
    gp_Vec dir;
    double minval = -1e99;
    double maxval = 1e99;
    bool openmin = false, openmax = false;

    DirectionalInterval (gp_Vec adir) : dir(adir) { ; }
    DirectionalInterval (const DirectionalInterval & i2)
      : dir(i2.dir), minval(i2.minval), maxval(i2.maxval) { ; }

    DirectionalInterval operator< (double val) const
    {
      DirectionalInterval i2 = *this;
      i2.maxval = val;
      return i2;
    }

    DirectionalInterval operator> (double val) const
    {
      DirectionalInterval i2 = *this;
      i2.minval = val;
      return i2;
    }


    DirectionalInterval Intersect (const DirectionalInterval & i2)
    {
      DirectionalInterval res = *this;
      res.minval = max(res.minval, i2.minval);
      res.maxval = min(res.maxval, i2.maxval);
      return res;
    }

    bool Contains (gp_Pnt p, double eps = 1e-8)
    {
      // cout << "Contains point " << p.X() << "," << p.Y() << "," << p.Z() << " ? " << endl;
      double val = dir.X()*p.X() + dir.Y()*p.Y() + dir.Z() * p.Z();
      // cout << "minval = " << minval << ", val = " << val << " maxval = " << maxval << endl;
      if (openmin) {
        if (val < minval+eps) return false;
      } else {
        if (val < minval-eps) return false;
      }
      if (openmax) {
        if (val > maxval-eps) return false;
      } else {
        if (val > maxval+eps) return false;
      }
      return true;
    }
  };


  inline gp_Pnt Center (TopoDS_Shape shape)
  {
    GProp_GProps props;
    switch (shape.ShapeType())
      {
      case TopAbs_FACE:
        BRepGProp::SurfaceProperties (shape, props); break;
      default:
        BRepGProp::LinearProperties(shape, props);
      }
    return props.CentreOfMass();
  }

}
#endif // FILE_OCC_UTILS_INCLUDED
