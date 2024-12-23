#ifndef FILE_BASEGEOM
#define FILE_BASEGEOM

/**************************************************************************/
/* File:   basegeom.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   23. Aug. 09                                                    */
/**************************************************************************/

#include <gprim/geomobjects.hpp>
#include <gprim/transform3d.hpp>

#include "meshtype.hpp"
#include "meshclass.hpp"

struct Tcl_Interp;

namespace netgen
{
  class Refinement;

  struct ShapeProperties
  {
    optional<string> name;
    optional<Vec<4>> col;
    double maxh = 1e99;
    double hpref = 0;  // number of hp refinement levels (will be multiplied by factor later)
    int layer = 1;
    optional<bool> quad_dominated;
    optional<Array<double>> partition;
    void Merge(const ShapeProperties & prop2)
    {
      if (!name && prop2.name) name = prop2.name;
      if (!col && prop2.col) col = prop2.col;
      if (!partition && prop2.partition) partition = prop2.partition;
      maxh = min2(maxh, prop2.maxh);
      hpref = max2(hpref, prop2.hpref);
      if(!quad_dominated.has_value()) quad_dominated = prop2.quad_dominated;
      layer = max(layer, prop2.layer);
    }

    string GetName() const { return name ? *name : "default"; }
    Vec<4> GetColor() { return col ? *col : Vec<4>{0., 1., 0., 1.}; }

    void DoArchive(Archive& ar)
    {
        ar & name & col & maxh & hpref & layer;
    }
  };

  class GeometryShape;

  struct ShapeIdentification
  {
    GeometryShape * from;
    GeometryShape * to;
    optional<Transformation<3>> trafo;
    Identifications::ID_TYPE type;
    string name = "";
  };

  class DLL_HEADER GeometryShape
  {
  public:
    int nr = -1;
    int layer = 1;
    ShapeProperties properties;
    Array<ShapeIdentification> identifications;
    GeometryShape * primary;
    optional<Transformation<3>> primary_to_me = nullopt;

    virtual ~GeometryShape() {}
    virtual bool IsMappedShape( const GeometryShape & other, const Transformation<3> & trafo, double tolerance ) const;
  };


  class DLL_HEADER GeometryVertex : public GeometryShape
  {
  public:
    virtual Point<3> GetPoint() const = 0;
    virtual bool IsMappedShape( const GeometryShape & other, const Transformation<3> & trafo, double tolerance ) const override;
  };

  class DLL_HEADER GeometryEdge : public GeometryShape
  {
  protected:
      GeometryVertex *start, *end;
  public:
    // Neighboring domains in 2d
    // In 3d unused, EXCEPT for free floating edges in a domain,
    // then both are pointing to the containing domain
    int domin=-1, domout=-1;

    GeometryEdge( GeometryVertex &start_, GeometryVertex &end_ )
        : start(&start_), end(&end_)
    {}

    virtual const GeometryVertex& GetStartVertex() const { return *start; }
    virtual const GeometryVertex& GetEndVertex() const { return *end; }
    virtual GeometryVertex& GetStartVertex() { return *start; }
    virtual GeometryVertex& GetEndVertex() { return *end; }
    virtual double GetLength() const = 0;
    virtual Point<3> GetCenter() const = 0;
    virtual Point<3> GetPoint(double t) const = 0;
    // Calculate parameter step respecting edges sag value
    virtual double CalcStep(double t, double sag) const = 0;
    virtual bool IsDegenerated(double tol = 1e-10) const {
      return GetLength() < tol;
    }
    virtual void ProjectPoint(Point<3>& p, EdgePointGeomInfo* gi) const = 0;
    virtual void PointBetween(const Point<3>& p1,
                              const Point<3>& p2,
                              double secpoint,
                              const EdgePointGeomInfo& gi1,
                              const EdgePointGeomInfo& gi2,
                              Point<3>& newp,
                              EdgePointGeomInfo& newgi) const
    {
      newp = p1 + secpoint * (p2-p1);
      newgi = gi1;
      ProjectPoint(newp, &newgi);
    }
    virtual Vec<3> GetTangent(double t) const = 0;
    virtual bool IsMappedShape( const GeometryShape & other, const Transformation<3> & trafo, double tolerance ) const override;
    virtual void Divide(const MeshingParameters & mparam, const Mesh & mesh, Array<Point<3>> & points, Array<double> & params);
  };

  class DLL_HEADER GeometryFace : public GeometryShape
  {
  public:
    Array<GeometryEdge*> edges;
    int domin=-1, domout=-1;

    virtual Point<3> GetCenter() const = 0;
    virtual size_t GetNBoundaries() const = 0;
    virtual Array<Segment> GetBoundary(const Mesh& mesh) const = 0;

    virtual PointGeomInfo Project(Point<3>& p) const = 0;
    // Project point using geo info. Fast if point is close to
    // parametrization in geo info.
    virtual bool ProjectPointGI(Point<3>& p, PointGeomInfo& gi) const =0;
    virtual bool CalcPointGeomInfo(const Point<3>& p, PointGeomInfo& gi) const
    {
      auto pnew = p;
      gi = Project(pnew);
      return (p-pnew).Length() < 1e-10 * GetBoundingBox().Diam() ;
    }
    virtual Point<3> GetPoint(const PointGeomInfo& gi) const = 0;
    virtual void CalcEdgePointGI(const GeometryEdge& edge,
                                 double t,
                                 EdgePointGeomInfo& egi) const = 0;
    virtual Box<3> GetBoundingBox() const = 0;

    // Get curvature in point from local coordinates in PointGeomInfo
    virtual double GetCurvature(const PointGeomInfo& gi) const = 0;

    virtual void RestrictH(Mesh& mesh, const MeshingParameters& mparam) const = 0;
    virtual Vec<3> GetNormal(const Point<3>& p, const PointGeomInfo* gi = nullptr) const = 0;

    virtual void PointBetween(const Point<3>& p1,
                              const Point<3>& p2,
                              double secpoint,
                              const PointGeomInfo& gi1,
                              const PointGeomInfo& gi2,
                              Point<3>& newp,
                              PointGeomInfo& newgi) const
    {
      newp = p1 + secpoint * (p2-p1);
      newgi.trignum = gi1.trignum;
      /*
      newgi.u = 0.5 * (gi1.u + gi1.u);
      newgi.v = 0.5 * (gi1.v + gi2.v);
      */
      newgi.u = gi1.u + secpoint*(gi2.u - gi1.u);
      newgi.v = gi1.v + secpoint*(gi2.v - gi1.v);
      if(!ProjectPointGI(newp, newgi))
        newgi = Project(newp);
    }

    virtual bool IsMappedShape( const GeometryShape & other, const Transformation<3> & trafo, double tolerance ) const override;
    virtual bool IsConnectingCloseSurfaces() const;

  protected:
    void RestrictHTrig(Mesh& mesh,
                       const PointGeomInfo& gi0,
                       const PointGeomInfo& gi1,
                       const PointGeomInfo& gi2,
                       const MeshingParameters& mparam,
                       int depth = 0, double h = 0.) const;
  };

  class DLL_HEADER GeometrySolid : public GeometryShape
  {
  public:
    Array<GeometryEdge*> free_edges; // edges with no adjacent face
  };

  class DLL_HEADER NetgenGeometry
  {
    unique_ptr<Refinement> ref;
  protected:
    Array<unique_ptr<GeometryVertex>> vertices;
    Array<unique_ptr<GeometryEdge>> edges;
    Array<unique_ptr<GeometryFace>> faces;
    Array<unique_ptr<GeometrySolid>> solids;
    Array<std::pair<Point<3>, double>> restricted_h;
    Box<3> bounding_box;
    int dimension = 3;

  public:

    NetgenGeometry()
    {
      ref = make_unique<Refinement>(*this);
    }
    virtual ~NetgenGeometry () { ; }

    size_t GetNVertices() const { return vertices.Size(); }
    size_t GetNEdges() const { return edges.Size(); }
    size_t GetNFaces() const { return faces.Size(); }
    size_t GetNSolids() const { return solids.Size(); }

    const GeometrySolid & GetSolid(int i) const { return *solids[i]; }
    const GeometryFace & GetFace(int i) const { return *faces[i]; }
    const GeometryEdge & GetEdge(int i) const { return *edges[i]; }
    const GeometryVertex & GetVertex(int i) const { return *vertices[i]; }

    auto Solids() const { return FlatArray{solids}; }
    auto Faces() const { return FlatArray{faces}; }
    auto Edges() const { return FlatArray{edges}; }
    auto Vertices() const { return FlatArray{vertices}; }

    virtual Array<const GeometryVertex*> GetFaceVertices(const GeometryFace& face) const { return Array<const GeometryVertex*>{}; }

    void Clear();

    virtual int GenerateMesh (shared_ptr<Mesh> & mesh, MeshingParameters & mparam);

    void RestrictH(const Point<3>& pnt, double maxh)
    {
      restricted_h.Append({pnt, maxh});
    }

    virtual const Refinement & GetRefinement () const
    {
      return *ref;
    }

    virtual void DoArchive(Archive&)
  { throw NgException("DoArchive not implemented for " + Demangle(typeid(*this).name())); }

    virtual Mesh::GEOM_TYPE GetGeomType() const { return Mesh::NO_GEOM; }
    virtual void ProcessIdentifications();
    virtual void Analyse(Mesh& mesh,
                         const MeshingParameters& mparam) const;
    virtual void FindEdges(Mesh& mesh, const MeshingParameters& mparam) const;
    virtual void MeshSurface(Mesh& mesh, const MeshingParameters& mparam) const;
    virtual bool MeshFace(Mesh& mesh, const MeshingParameters& mparam,
                     int nr, FlatArray<int, PointIndex> glob2loc) const;
    virtual void MapSurfaceMesh( Mesh & mesh, const GeometryFace & dst, std::map<tuple<PointIndex, int>, PointIndex> & mapto) const;
    virtual void OptimizeSurface(Mesh& mesh, const MeshingParameters& mparam) const;

    virtual void FinalizeMesh(Mesh& mesh) const;

    virtual PointGeomInfo ProjectPoint (int surfind, Point<3> & p) const
    {
      if(surfind <= faces.Size() && surfind > 0)
        return faces[surfind-1]->Project(p);
      return PointGeomInfo();
    }

  virtual void ProjectPointEdge (int surfind, int surfind2, Point<3> & p, EdgePointGeomInfo* gi = nullptr) const
  {
    if(gi && gi->edgenr < edges.Size() && gi->edgenr >= 0)
      edges[gi->edgenr]->ProjectPoint(p, gi);
  }

    virtual bool CalcPointGeomInfo(int surfind, PointGeomInfo& gi, const Point<3> & p3) const
    {
      return faces[surfind-1]->CalcPointGeomInfo(p3, gi);
    }
    virtual bool ProjectPointGI (int surfind, Point<3> & p, PointGeomInfo & gi) const
    {
      if(surfind > 0 && surfind <= faces.Size())
        return faces[surfind-1]->ProjectPointGI(p, gi);
      return false;
    }

    virtual Vec<3> GetNormal(int surfind, const Point<3> & p, const PointGeomInfo* gi = nullptr) const
    {
      if(surfind > 0 && surfind <= faces.Size())
        return faces[surfind-1]->GetNormal(p, gi);
      else
        return {0., 0., 0.};
    }

    virtual void PointBetween (const Point<3> & p1,
                               const Point<3> & p2, double secpoint,
                               int surfi,
                               const PointGeomInfo & gi1,
                               const PointGeomInfo & gi2,
                               Point<3> & newp,
                               PointGeomInfo & newgi) const
    {
      if(faces.Size() >= surfi && surfi > 0)
        {
          faces[surfi-1]->PointBetween(p1, p2, secpoint, gi1, gi2, newp, newgi);
          return;
        }
      newp = p1 + secpoint * (p2-p1);
    }

    virtual void PointBetweenEdge(const Point<3> & p1,
                                  const Point<3> & p2, double secpoint,
                                  int surfi1, int surfi2,
                                  const EdgePointGeomInfo & ap1,
                                  const EdgePointGeomInfo & ap2,
                                  Point<3> & newp,
                                  EdgePointGeomInfo & newgi) const
    {
      if(ap1.edgenr < edges.Size() && ap1.edgenr >= 0)
        {
          edges[ap1.edgenr]->PointBetween(p1, p2, secpoint,
                                          ap1, ap2, newp, newgi);
          return;
        }
      newp = p1+secpoint*(p2-p1);
    }

    virtual Vec<3> GetTangent(const Point<3> & p, int surfi1,
                              int surfi2,
                              const EdgePointGeomInfo & egi) const
    {
      throw Exception("Base geometry get tangent called");
    }

    virtual void Save (const filesystem::path & filename) const;
    virtual void SaveToMeshFile (ostream & /* ost */) const { ; }
  };





  class DLL_HEADER GeometryRegister
  {
  public:
    virtual ~GeometryRegister();
    virtual NetgenGeometry * Load (const filesystem::path & filename) const = 0;
    virtual NetgenGeometry * LoadFromMeshFile (istream & /* ist */, string) const { return NULL; }
    virtual class VisualScene * GetVisualScene (const NetgenGeometry * /* geom */) const
    { return NULL; }
    virtual void SetParameters (Tcl_Interp * /* interp */) { ; }
  };

  class DLL_HEADER GeometryRegisterArray : public NgArray<GeometryRegister*>
  {
  public:
    virtual ~GeometryRegisterArray()
    {
      for (int i = 0; i < Size(); i++)
        delete (*this)[i];
    }

    virtual shared_ptr<NetgenGeometry> LoadFromMeshFile (istream & ist) const;
  };

  // extern DLL_HEADER NgArray<GeometryRegister*> geometryregister; 
  extern DLL_HEADER GeometryRegisterArray geometryregister; 
}



#endif
