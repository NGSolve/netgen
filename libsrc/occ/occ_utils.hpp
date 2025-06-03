#ifndef FILE_OCC_UTILS_INCLUDED
#define FILE_OCC_UTILS_INCLUDED



#include <variant>

// #pragma clang diagnostic push
// #pragma clang diagnostic ignored "-Wdeprecated-declarations"

#include <Standard_Version.hxx>
#include <BRepGProp.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Vertex.hxx>
#include <gp_Trsf.hxx>
#include <gp_GTrsf.hxx>

#define NETGEN_OCC_VERSION_AT_LEAST(MAJOR, MINOR) \
  ((OCC_VERSION_MAJOR > MAJOR) ||                               \
   ((OCC_VERSION_MAJOR == MAJOR) && (OCC_VERSION_MINOR >= MINOR)))
#define NETGEN_OCC_VERSION_AT_LEAST_MAJOR(MAJOR) \
  (NETGEN_OCC_VERSION_AT_LEAST(MAJOR, 0))

// #pragma clang diagnostic pop

#include "meshing.hpp"

#if NETGEN_OCC_VERSION_AT_LEAST(7, 4)
#define OCC_HAVE_DUMP_JSON
#endif

namespace netgen
{
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

    DLL_HEADER Point<3> occ2ng (const TopoDS_Shape & s);

    inline Point<3> occ2ng (const TopoDS_Vertex & v)
    {
        return occ2ng (BRep_Tool::Pnt (v));
    }

    DLL_HEADER Transformation<3> occ2ng (const gp_Trsf & t);
    DLL_HEADER Transformation<3> occ2ng (const gp_GTrsf & t);
    inline Transformation<3> occ2ng (const variant<gp_Trsf, gp_GTrsf> & t)
    {
      if(auto t1 = get_if<gp_Trsf>(&t))
        return occ2ng(*t1);
      return occ2ng(get<gp_GTrsf>(t));
    }

    inline gp_Pnt ng2occ (const Point<3> & p)
    {
        return gp_Pnt(p(0), p(1), p(2));
    }

    DLL_HEADER Box<3> GetBoundingBox( const TopoDS_Shape & shape );

    struct OCCIdentification
    {
      TopoDS_Shape from;
      TopoDS_Shape to;
      optional<Transformation<3>> trafo = nullopt;
      string name;
      Identifications::ID_TYPE type;
      bool opposite_direction = false;
    };

    Standard_Integer BuildTriangulation( const TopoDS_Shape & shape );


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

      ListOfShapes GetHighestDimShapes() const
      {
        for (auto type : {TopAbs_SOLID, TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX})
        {
          auto ret = SubShapes(type);
          if (ret.size() > 0)
            return ret;
        }
        return ListOfShapes();
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

    inline ListOfShapes GetHighestDimShapes(const TopoDS_Shape & shape)
    {
      auto ret = GetSolids(shape); if(ret.size() > 0) return ret;
      ret = GetFaces(shape); if(ret.size() > 0) return ret;
      ret = GetEdges(shape); if(ret.size() > 0) return ret;
      ret = GetVertices(shape); if(ret.size() > 0) return ret;
      return ListOfShapes();
    }


  class DirectionalInterval
  {
  public:
    gp_Vec dir;
    double minval = -1e99;
    double maxval = 1e99;
    bool openmin = false, openmax = false;

    DirectionalInterval (gp_Vec adir) : dir(adir) { ; }
    DirectionalInterval (const DirectionalInterval & i2) = default;

    DirectionalInterval operator< (double val) const
    {
      DirectionalInterval i2 = *this;
      i2.maxval = val;
      i2.openmax = true;
      return i2;
    }

    DirectionalInterval operator> (double val) const
    {
      DirectionalInterval i2 = *this;
      i2.minval = val;
      i2.openmin = true;
      return i2;
    }

    DirectionalInterval operator<= (double val) const
    {
      DirectionalInterval i2 = *this;
      i2.maxval = val;
      i2.openmax = false;
      return i2;
    }

    DirectionalInterval operator>= (double val) const
    {
      DirectionalInterval i2 = *this;
      i2.minval = val;
      i2.openmin = false;
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

  inline auto Properties (TopoDS_Shape shape)
  {
    GProp_GProps props;
    double tol;
    switch (shape.ShapeType())
    {
      case TopAbs_SOLID:
      case TopAbs_COMPOUND:
      case TopAbs_COMPSOLID:
        tol = 1e-2 * BRep_Tool::MaxTolerance(shape, TopAbs_FACE);
        BRepGProp::VolumeProperties (shape, props, tol); break;
      case TopAbs_FACE:
      case TopAbs_SHELL:
        tol = 1e-2 * BRep_Tool::MaxTolerance(shape, TopAbs_FACE);
        BRepGProp::SurfaceProperties (shape, props, tol); break;
      case TopAbs_WIRE:
      case TopAbs_EDGE:
        tol = 1e-2 * BRep_Tool::MaxTolerance(shape, TopAbs_EDGE);
        BRepGProp::LinearProperties(shape, props, tol); break;
      default:
        BRepGProp::LinearProperties(shape, props);
    }
    return props;
  }

  inline gp_Pnt Center (TopoDS_Shape shape)
  {
    return Properties(shape).CentreOfMass();
  }

  inline double Mass (TopoDS_Shape shape)
  {
    return Properties(shape).Mass();
  }
  
}
#endif // FILE_OCC_UTILS_INCLUDED
