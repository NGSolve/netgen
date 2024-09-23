#ifndef FILE_OCCGEOM
#define FILE_OCCGEOM

/* *************************************************************************/
/* File:   occgeom.hpp                                                     */
/* Author: Robert Gaisbauer                                                */
/* Date:   26. May  03                                                     */
/* *************************************************************************/

#ifdef OCCGEOMETRY

#include <set>

#include <meshing.hpp>
#include "occ_utils.hpp"
#include "occmeshsurf.hpp"

#include <BOPAlgo_BuilderShape.hxx>
#include <BRepTools_ReShape.hxx>
#include <BRepBuilderAPI_MakeShape.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <Quantity_ColorRGBA.hxx>
#include <STEPCAFControl_Reader.hxx>
#include <StepBasic_MeasureValueMember.hxx>
#include <StepRepr_CompoundRepresentationItem.hxx>
#include <StepRepr_IntegerRepresentationItem.hxx>
#include <StepRepr_ValueRepresentationItem.hxx>
#include <TCollection_HAsciiString.hxx>
#include <TDocStd_Document.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <Transfer_FinderProcess.hxx>

#if OCC_VERSION_MAJOR>=7 && OCC_VERSION_MINOR>=4
#define OCC_HAVE_HISTORY
#endif

namespace netgen
{

  // extern DLL_HEADER MeshingParameters mparam;

#define PROJECTION_TOLERANCE 1e-10

#define ENTITYISVISIBLE 1
#define ENTITYISHIGHLIGHTED 2
#define ENTITYISDRAWABLE 4

#define OCCGEOMETRYVISUALIZATIONNOCHANGE   0
#define OCCGEOMETRYVISUALIZATIONFULLCHANGE 1  // Compute transformation matrices and redraw
#define OCCGEOMETRYVISUALIZATIONHALFCHANGE 2  // Redraw

  bool IsMappedShape(const Transformation<3> & trafo, const TopoDS_Shape & me, const TopoDS_Shape & you);

  class EntityVisualizationCode
  {
    int code;

  public:

    EntityVisualizationCode()
    {  code = ENTITYISVISIBLE + !ENTITYISHIGHLIGHTED + ENTITYISDRAWABLE;}

    int IsVisible ()
    {  return code & ENTITYISVISIBLE;}

    int IsHighlighted ()
    {  return code & ENTITYISHIGHLIGHTED;}

    int IsDrawable ()
    {  return code & ENTITYISDRAWABLE;}

    void Show ()
    {  code |= ENTITYISVISIBLE;}

    void Hide ()
    {  code &= ~ENTITYISVISIBLE;}

    void Highlight ()
    {  code |= ENTITYISHIGHLIGHTED;}

    void Lowlight ()
    {  code &= ~ENTITYISHIGHLIGHTED;}

    void SetDrawable ()
    {  code |= ENTITYISDRAWABLE;}

    void SetNotDrawable ()
    {  code &= ~ENTITYISDRAWABLE;}
  };



  class Line
  {
  public:
    Point<3> p0, p1;
    int layer = 1;
    double Dist (Line l);
    double Length () { return (p1-p0).Length(); }
  };
  


  inline double Det3 (double a00, double a01, double a02,
                      double a10, double a11, double a12,
                      double a20, double a21, double a22)
  {
    return a00*a11*a22 + a01*a12*a20 + a10*a21*a02 - a20*a11*a02 - a10*a01*a22 - a21*a12*a00;
  }
  
  class DLL_HEADER OCCParameters
  {
  public:

    /// Factor for meshing close edges, moved to meshingparameters
    // double resthcloseedgefac = 2.;

    /// Enable / Disable detection of close edges
    // int resthcloseedgeenable = true;

    /// Minimum edge length to be used for dividing edges to mesh points
    // double resthminedgelen = 0.001;
    double resthminedgelen = 1e-4;

    /// Enable / Disable use of the minimum edge length (by default use 1e-4)
    int resthminedgelenenable = false;

    /*!
      Dump all the OpenCascade specific meshing parameters 
      to console
    */
    void Print (ostream & ost) const;
  };


  class DLL_HEADER OCCGeometry : public NetgenGeometry
  {
    Point<3> center;
    OCCParameters occparam;
  public:
    static TopTools_IndexedMapOfShape global_shape_property_indices;
    static std::vector<ShapeProperties> global_shape_properties;
    static TopTools_IndexedMapOfShape global_identification_indices;
    static std::vector<std::vector<OCCIdentification>> global_identifications;

    static ShapeProperties& GetProperties(const TopoDS_Shape& shape)
    {
      auto index = OCCGeometry::global_shape_property_indices.FindIndex(shape);
      if(index > 0)
        return OCCGeometry::global_shape_properties
          [index-1];
      OCCGeometry::global_shape_property_indices.Add(shape);
      OCCGeometry::global_shape_properties.push_back({});
      return OCCGeometry::global_shape_properties.back();
    }
    static bool HaveProperties(const TopoDS_Shape& shape)
    {
      return OCCGeometry::global_shape_property_indices.FindIndex(shape) > 0;
    }
    static std::vector<OCCIdentification>& GetIdentifications(const TopoDS_Shape& shape)
    {
      auto index = OCCGeometry::global_identification_indices.FindIndex(shape);
      if(index > 0)
        return OCCGeometry::global_identifications[index-1];
      OCCGeometry::global_identification_indices.Add(shape);
      OCCGeometry::global_identifications.push_back({});
      return OCCGeometry::global_identifications.back();
    }
    static bool HaveIdentifications(const TopoDS_Shape& shape)
    {
      return OCCGeometry::global_identification_indices.FindIndex(shape) > 0;
    }

    TopoDS_Shape shape;
    TopTools_IndexedMapOfShape fmap, emap, vmap, somap, shmap, wmap;
    NgArray<bool> fsingular, esingular, vsingular;
    Box<3> boundingbox;

    mutable int changed;
    mutable NgArray<int> facemeshstatus;

    // Philippose - 15/01/2009
    // Maximum mesh size for a given face
    // (Used to explicitly define mesh size limits on individual faces)
    NgArray<double> face_maxh;
     
    // Philippose - 14/01/2010
    // Boolean array to detect whether a face has been explicitly modified 
    // by the user or not
    NgArray<bool> face_maxh_modified;
     
    // Philippose - 15/01/2009
    // Indicates which faces have been selected by the user in geometry mode
    // (Currently handles only selection of one face at a time, but an array would
    //  help to extend this to multiple faces)
    NgArray<bool> face_sel_status;
     
    NgArray<EntityVisualizationCode> fvispar, evispar, vvispar;
     
    double tolerance;
    bool fixsmalledges;
    bool fixspotstripfaces;
    bool sewfaces;
    bool makesolids;
    bool splitpartitions;

    OCCGeometry()
    {
      somap.Clear();
      shmap.Clear();
      fmap.Clear();
      wmap.Clear();
      emap.Clear();
      vmap.Clear();
    }

    OCCGeometry(const TopoDS_Shape& _shape, int aoccdim = 3, bool copy = false);

    Mesh::GEOM_TYPE GetGeomType() const override
    { return Mesh::GEOM_OCC; }

    void SetDimension(int dim)
    {
        dimension = dim;
        BuildFMap();
    }

    void SetOCCParameters(const OCCParameters& par)
    { occparam = par; }

    using NetgenGeometry::GetVertex;
    using NetgenGeometry::GetEdge;
    using NetgenGeometry::GetFace;
    using NetgenGeometry::GetSolid;

    GeometryShape & GetShape(const TopoDS_Shape & shape)
    {
        return const_cast<GeometryShape&>(as_const(*this).GetShape(shape));
    }
    GeometryVertex & GetVertex(const TopoDS_Shape & shape)
    {
        return const_cast<GeometryVertex&>(as_const(*this).GetVertex(shape));
    }

    GeometryEdge & GetEdge(const TopoDS_Shape & shape)
    {
        return const_cast<GeometryEdge&>(as_const(*this).GetEdge(shape));
    }

    GeometryFace & GetFace(const TopoDS_Shape & shape)
    {
        return const_cast<GeometryFace&>(as_const(*this).GetFace(shape));
    }

    GeometrySolid & GetSolid(const TopoDS_Shape & shape)
    {
        return const_cast<GeometrySolid&>(as_const(*this).GetSolid(shape));
    }

    const GeometryShape & GetShape(const TopoDS_Shape & shape) const;
    const GeometryVertex & GetVertex(const TopoDS_Shape & shape) const;
    const GeometryEdge & GetEdge(const TopoDS_Shape & shape) const;
    const GeometryFace & GetFace(const TopoDS_Shape & shape) const;
    const GeometrySolid & GetSolid(const TopoDS_Shape & shape) const;

    void Analyse(Mesh& mesh,
                 const MeshingParameters& mparam) const override;
    bool MeshFace(Mesh& mesh, const MeshingParameters& mparam,
                     int nr, FlatArray<int, PointIndex> glob2loc) const override;
    // void OptimizeSurface(Mesh& mesh, const MeshingParameters& mparam) const override {}
 
    void Save (const filesystem::path & filename) const override;
    void SaveToMeshFile (ostream & /* ost */) const override;
     
    void DoArchive(Archive& ar) override;

    void BuildFMap();

    auto GetShape() const { return shape; }
    Box<3> GetBoundingBox() const
    { return boundingbox; }

    int NrSolids() const
    { return somap.Extent(); }

    // Philippose - 17/01/2009
    // Total number of faces in the geometry
    int NrFaces() const
    { return fmap.Extent(); }

    void SetCenter()
    { center = boundingbox.Center(); }

    Point<3> Center() const
    { return center; }

    OCCSurface GetSurface (int surfi)
    {
      cout << "OCCGeometry::GetSurface using PLANESPACE" << endl;
      return OCCSurface (TopoDS::Face(fmap(surfi)), PLANESPACE);
    }

    void CalcBoundingBox ();
    void BuildVisualizationMesh (double deflection);
    
    void RecursiveTopologyTree (const TopoDS_Shape & sh,
                                stringstream & str,
                                TopAbs_ShapeEnum l,
                                bool free,
                                const char * lname);

    void GetTopologyTree (stringstream & str);

    void PrintNrShapes ();

    void CheckIrregularEntities (stringstream & str);

    void SewFaces();

    void MakeSolid();

    Array<const GeometryVertex*> GetFaceVertices(const GeometryFace& face) const override;

    void FixFaceOrientation();
    void HealGeometry();
    void GlueGeometry();

    // Philippose - 15/01/2009
    // Sets the maximum mesh size for a given face
    // (Note: Local mesh size limited by the global max mesh size)
    void SetFaceMaxH(int facenr, double faceh, const MeshingParameters & mparam)
    {
      if((facenr> 0) && (facenr <= fmap.Extent()))
        {
          face_maxh[facenr-1] = min(mparam.maxh,faceh);
            
          // Philippose - 14/01/2010
          // If the face maxh is greater than or equal to the 
          // current global maximum, then identify the face as 
          // not explicitly controlled by the user any more
          if(faceh >= mparam.maxh)
            {
              face_maxh_modified[facenr-1] = 0;
            }
          else
            {
              face_maxh_modified[facenr-1] = 1;
            }
        }
    }

    void SetFaceMaxH(size_t facenr, double faceh)
    {
      if(facenr >= fmap.Extent())
        throw RangeException("OCCGeometry faces", facenr, 0, fmap.Extent());
      face_maxh[facenr] = faceh;
      face_maxh_modified[facenr] = true;
    }

    // Philippose - 15/01/2009
    // Returns the local mesh size of a given face
    double GetFaceMaxH(int facenr)
    {
      if((facenr> 0) && (facenr <= fmap.Extent()))
        {
          return face_maxh[facenr-1];
        }
      else
        {
          return 0.0;
        }
    }
      
    // Philippose - 14/01/2010
    // Returns the flag whether the given face 
    // has a mesh size controlled by the user or not
    bool GetFaceMaxhModified(int facenr)
    {
      return face_maxh_modified[facenr-1];
    }
      
    // Philippose - 17/01/2009
    // Returns the index of the currently selected face
    int SelectedFace()
    {
      for(int i = 1; i <= fmap.Extent(); i++)
        {
          if(face_sel_status[i-1])
            {
              return i;
            }
        }

      return 0;
    }

    // Philippose - 17/01/2009
    // Sets the currently selected face
    void SetSelectedFace(int facenr)
    {
      face_sel_status = 0;

      if((facenr >= 1) && (facenr <= fmap.Extent()))
        {
          face_sel_status[facenr-1] = 1;
        }
    }

    void LowLightAll()
    {
      for (int i = 1; i <= fmap.Extent(); i++)
        fvispar[i-1].Lowlight();
      for (int i = 1; i <= emap.Extent(); i++)
        evispar[i-1].Lowlight();
      for (int i = 1; i <= vmap.Extent(); i++)
        vvispar[i-1].Lowlight();
    }

    void GetUnmeshedFaceInfo (stringstream & str);
    void GetNotDrawableFaces (stringstream & str);
    bool ErrorInSurfaceMeshing ();

    //      void WriteOCC_STL(char * filename);

  private:
    //bool FastProject (int surfi, Point<3> & ap, double& u, double& v) const;
  };

  DLL_HEADER void Identify(const ListOfShapes & me, const ListOfShapes & you, string name, Identifications::ID_TYPE type, Transformation<3> trafo);
  DLL_HEADER void Identify(const TopoDS_Shape & me, const TopoDS_Shape & you, string name, Identifications::ID_TYPE type, std::optional<std::variant<gp_Trsf, gp_GTrsf>> opt_trafo);
   

  void PrintContents (OCCGeometry * geom);

  DLL_HEADER OCCGeometry * LoadOCC_IGES (const filesystem::path & filename);
  DLL_HEADER OCCGeometry * LoadOCC_STEP (const filesystem::path & filename);
  DLL_HEADER OCCGeometry * LoadOCC_BREP (const filesystem::path & filename);

  // Philippose - 31.09.2009
  // External access to the mesh generation functions within the OCC
  // subsystem (Not sure if this is the best way to implement this....!!)
  DLL_HEADER extern void OCCSetLocalMeshSize(const OCCGeometry & geom, Mesh & mesh, const MeshingParameters & mparam,
                                             const OCCParameters& occparam);

  DLL_HEADER extern bool OCCMeshFace (const OCCGeometry & geom, Mesh & mesh, FlatArray<int, PointIndex> glob2loc,
                       const MeshingParameters & mparam, int nr, int projecttype, bool delete_on_failure);

  inline auto GetModified(BRepBuilderAPI_MakeShape & builder, TopoDS_Shape shape) { return builder.Modified(shape); }
  inline auto GetModified(BRepTools_History & history, TopoDS_Shape shape) { return history.Modified(shape); }
  inline auto GetModified(BOPAlgo_BuilderShape & builder, TopoDS_Shape shape) { return builder.Modified(shape); }
  inline ArrayMem<TopoDS_Shape, 1> GetModified(BRepBuilderAPI_Sewing& builder, TopoDS_Shape shape) { return {builder.Modified(shape)}; }
  inline auto GetModified(BRepTools_ReShape& reshape, TopoDS_Shape shape) {
    auto history = reshape.History();
    return history->Modified(shape);
  }

  template <class TBuilder>
  void PropagateIdentifications (TBuilder & builder, TopoDS_Shape shape, std::optional<Transformation<3>> trafo = nullopt)
  {
    TopTools_IndexedMapOfShape mod_indices;
    Array<TopTools_IndexedMapOfShape> modifications;
    TopTools_MapOfShape shape_handled;

    Transformation<3> trafo_inv;
    if(trafo)
        trafo_inv = trafo->CalcInverse();
  
    for (auto typ : { TopAbs_SOLID, TopAbs_FACE,  TopAbs_EDGE, TopAbs_VERTEX })
      for (TopExp_Explorer e(shape, typ); e.More(); e.Next())
      {
          auto s = e.Current();
          mod_indices.Add(s);
          modifications.Append(TopTools_IndexedMapOfShape());
          modifications.Last().Add(s);
      }
  
    for (auto typ : { TopAbs_SOLID, TopAbs_FACE,  TopAbs_EDGE, TopAbs_VERTEX })
      for (TopExp_Explorer e(shape, typ); e.More(); e.Next())
        {
          auto s = e.Current();
            for (auto mods : GetModified(builder, s))
              {
                auto index = mod_indices.FindIndex(s)-1;
                modifications[index].Add(mods);
              }
        }
  
    for (auto typ : { TopAbs_SOLID, TopAbs_FACE,  TopAbs_EDGE, TopAbs_VERTEX })
      for (TopExp_Explorer e(shape, typ); e.More(); e.Next())
      {
          auto s = e.Current();
  
          if(shape_handled.Contains(s))
              continue;
          shape_handled.Add(s);

          if(!OCCGeometry::HaveIdentifications(s))
            continue;
          auto& identifications = OCCGeometry::GetIdentifications(s);

          // auto& shape_mapped = modifications[mod_indices.FindIndex(s)-1];
  
          for(auto ident : identifications)
          {
              auto i1 = mod_indices.FindIndex(ident.to);
              auto i2 = mod_indices.FindIndex(ident.from);
              if(i1 == 0 || i2 == 0) // not in geometry
                continue;
              auto& mods_to = modifications[i1-1];
              auto& mods_from = modifications[i2-1];
              if(mods_to.Extent()==1 && mods_from.Extent() ==1)
                continue;
  
              auto from = ident.from;
              auto to = ident.to;
  
              for(auto it = mods_from.cbegin(); it != mods_from.cend(); it++)
                for(auto it2 = mods_to.cbegin(); it2 != mods_to.cend(); it2++)
                  {
                      auto& from_mapped = it.Iterator().Value();
                      auto& to_mapped = it2.Iterator().Value();
                      if(from.IsSame(from_mapped) && to.IsSame(to_mapped))
                        continue;
  
                      if(!ident.trafo) continue;
                      Transformation<3> trafo_mapped = *ident.trafo;

                      if(trafo)
                      {
                          Transformation<3> trafo_temp;
                          trafo_temp.Combine(*ident.trafo, trafo_inv);
                          trafo_mapped.Combine(*trafo, trafo_temp);
                      }
  
                      if(!IsMappedShape(trafo_mapped, from_mapped, to_mapped))
                          continue;
  
                      OCCIdentification id_new = ident;
                      id_new.to = to_mapped;
                      id_new.from = from_mapped;
                      id_new.trafo = trafo_mapped;
                      auto id_owner = from.IsSame(s) ? from_mapped : to_mapped;
                      OCCGeometry::GetIdentifications(id_owner).push_back(id_new);
                  }
          }
      }
  }
  
  template <class TBuilder>
  void PropagateProperties (TBuilder & builder, TopoDS_Shape shape, std::optional<Transformation<3>> trafo = nullopt)
  {
    bool have_identifications = false;
  
    for (auto typ : { TopAbs_SOLID, TopAbs_FACE,  TopAbs_EDGE, TopAbs_VERTEX })
      for (TopExp_Explorer e(shape, typ); e.More(); e.Next())
        {
          auto s = e.Current();
          have_identifications |= OCCGeometry::HaveIdentifications(s);
          if(!OCCGeometry::HaveProperties(s))
            continue;
          auto prop = OCCGeometry::GetProperties(s);
          for (auto mods : GetModified(builder, s))
            OCCGeometry::GetProperties(mods).Merge(prop);
        }
    if(have_identifications)
        PropagateIdentifications(builder, shape, trafo);
  }

  namespace step_utils
  {
      inline Handle(TCollection_HAsciiString) MakeName (string s)
      {
          return new TCollection_HAsciiString(s.c_str());
      };

      inline Handle(StepRepr_RepresentationItem) MakeInt (int n, string name = "")
      {
          Handle(StepRepr_IntegerRepresentationItem) int_obj = new StepRepr_IntegerRepresentationItem;
          int_obj->Init(MakeName(name), n);
          return int_obj;
      }

      inline int ReadInt (Handle(StepRepr_RepresentationItem) item)
      {
          return Handle(StepRepr_IntegerRepresentationItem)::DownCast(item)->Value();
      }

      inline Handle(StepRepr_RepresentationItem) MakeReal (double val, string name = "")
      {
            Handle(StepBasic_MeasureValueMember) value_member = new StepBasic_MeasureValueMember;
            value_member->SetReal(val);
            Handle(StepRepr_ValueRepresentationItem) value_repr = new StepRepr_ValueRepresentationItem;
            value_repr->Init(MakeName(name), value_member);
            return value_repr;
      }

      inline double ReadReal (Handle(StepRepr_RepresentationItem) item)
      {
          return Handle(StepRepr_ValueRepresentationItem)::DownCast(item)
              ->ValueComponentMember()->Real();
      }


      inline Handle(StepRepr_RepresentationItem) MakeCompound( FlatArray<Handle(StepRepr_RepresentationItem)> items, string name = "" )
      {
            Handle(StepRepr_HArray1OfRepresentationItem) array_repr = new StepRepr_HArray1OfRepresentationItem(1,items.Size());

            for(auto i : Range(items))
                array_repr->SetValue(i+1, items[i]);

            Handle(StepRepr_CompoundRepresentationItem) comp = new StepRepr_CompoundRepresentationItem;
            comp->Init( MakeName(name), array_repr );
            return comp;
      }

      void WriteIdentifications(const Handle(Interface_InterfaceModel) model, const TopoDS_Shape & shape, const Handle(Transfer_FinderProcess) finder);
      void ReadIdentifications(Handle(StepRepr_RepresentationItem) item, Handle(Transfer_TransientProcess) transProc);

      inline Quantity_ColorRGBA MakeColor(const Vec<4> & c)
      {
          return Quantity_ColorRGBA (c[0], c[1], c[2], c[3]);
      }

      inline Vec<4> ReadColor (const Quantity_ColorRGBA & c)
      {
          auto rgb = c.GetRGB();
          return {rgb.Red(), rgb.Green(), rgb.Blue(), c.Alpha()};
      }


      void LoadProperties(const TopoDS_Shape & shape,
                          const STEPCAFControl_Reader & reader,
                          const Handle(TDocStd_Document) step_doc);
      void WriteProperties(const Handle(Interface_InterfaceModel) model, const Handle(Transfer_FinderProcess) finder, const TopoDS_Shape & shape);

      void WriteSTEP(const TopoDS_Shape & shape, const filesystem::path & filename);

      inline void WriteSTEP(const OCCGeometry & geo, const filesystem::path & filename)
      {
          WriteSTEP(geo.GetShape(), filename);
      }

      // deep copy, also ensures consistent shape ordering (face numbers etc.)
      TopoDS_Shape WriteAndRead(const TopoDS_Shape shape);
  } // namespace step_utils
}

#endif

#endif
