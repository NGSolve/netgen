
#include <BOPAlgo_Builder.hxx>
#include <BRepAlgoAPI_Splitter.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <Geom_Plane.hxx>
#include <Geom_Surface.hxx>
#include <Precision.hxx>
#include <ShapeBuild_ReShape.hxx>
#include <ShapeFix_Shape.hxx>
#include <Standard_Handle.hxx>
#include <Standard_TypeDef.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <algorithm>
#include <gp.hxx>
#include <gp_Ax3.hxx>
#include <gp_Pln.hxx>
#include <gp_Trsf.hxx>
#include <limits>

#include "core/array.hpp"
#include "occ_utils.hpp"
#include "occgeom.hpp"

namespace netgen
{
  TopoDS_Face MakePlaneFaceCoveringShape
  (const TopoDS_Shape& shape,
   const gp_Pln& pln,
   double rel_margin = 0.05,
   double abs_margin = 0.0)
  {
    Bnd_Box b;
    BRepBndLib::Add(shape, b);
    b.SetGap(0.0);

    Standard_Real xmin, ymin, zmin, xmax, ymax, zmax;
    b.Get(xmin, ymin, zmin, xmax, ymax, zmax);

    gp_Trsf w2l;
    w2l.SetTransformation(pln.Position());

    auto upd = [&] (double x, double y, double z, double& umin, double& umax, double& vmin, double& vmax) {
      gp_Pnt p(x, y, z);
      p.Transform(w2l);
      umin = std::min(umin, p.X());
      umax = std::max(umax, p.X());
      vmin = std::min(vmin, p.Y());
      vmax = std::max(vmax, p.Y());
    };

    double umin = std::numeric_limits<double>::infinity();
    double umax = -std::numeric_limits<double>::infinity();
    double vmin = std::numeric_limits<double>::infinity();
    double vmax = -std::numeric_limits<double>::infinity();

    upd(xmin, ymin, zmin, umin, umax, vmin, vmax);
    upd(xmin, ymin, zmax, umin, umax, vmin, vmax);
    upd(xmin, ymax, zmin, umin, umax, vmin, vmax);
    upd(xmin, ymax, zmax, umin, umax, vmin, vmax);
    upd(xmax, ymin, zmin, umin, umax, vmin, vmax);
    upd(xmax, ymin, zmax, umin, umax, vmin, vmax);
    upd(xmax, ymax, zmin, umin, umax, vmin, vmax);
    upd(xmax, ymax, zmax, umin, umax, vmin, vmax);

    double du = umax - umin;
    double dv = vmax - vmin;
    double pad = std::max(abs_margin, rel_margin * std::max(du, dv));

    umin -= pad;
    umax += pad;
    vmin -= pad;
    vmax += pad;

    double eps = std::max(Precision::Confusion(), 1e-9);
    if (umax - umin < eps)
      {
        umin -= 1.0;
        umax += 1.0;
      }
    if (vmax - vmin < eps)
      {
        vmin -= 1.0;
        vmax += 1.0;
      }
    return BRepBuilderAPI_MakeFace(pln, umin, umax, vmin, vmax).Face();
  }

  TopoDS_Shape CrossSection (const TopoDS_Shape& shape,
                             const gp_Ax3& axes)
  {
    gp_Pln pln(axes);
    TopoDS_Face planeFace = MakePlaneFaceCoveringShape(shape, axes);

    auto bb = GetBoundingBox(shape);
    auto diam = bb.Diam();

    // Heal / unify the input
    //  - especially if shape comes from STEP/IGES
    ShapeFix_Shape fix(shape);
    fix.Perform();
    TopoDS_Shape fixed = fix.Shape();

    BRepAlgoAPI_Splitter splitter;
    splitter.SetNonDestructive(true);
    splitter.SetFuzzyValue(1e-7 * diam);

    TopTools_ListOfShape args;
    args.Append(shape);

    TopTools_ListOfShape tools;
    tools.Append(planeFace);

    splitter.SetArguments(args);
    splitter.SetTools(tools);
    splitter.Build();
    if (!splitter.IsDone())
      return TopoDS_Shape();

    // for each input face/edge propagate properties to generated sub-shapes
    for (auto typ : {TopAbs_FACE, TopAbs_EDGE})
      for (TopExp_Explorer e(shape, typ); e.More(); e.Next())
        {
          if (!OCCGeometry::HaveProperties(e.Current()))
            continue;
          auto prop = OCCGeometry::GetProperties(e.Current());
          for (auto mods : splitter.Generated(e.Current()))
            OCCGeometry::GetProperties(mods).Merge(prop);
        }
    // for each input solid, find the cut faces and propagate properties
    for (auto solid : GetSolids(shape))
      {
        if (!OCCGeometry::HaveProperties(solid))
          continue;
        auto prop = OCCGeometry::GetProperties(solid);
        for (auto mods : splitter.Modified(solid))
          {
            if (mods.ShapeType() != TopAbs_SOLID)
              continue;
            for (TopExp_Explorer ef(mods, TopAbs_FACE); ef.More(); ef.Next())
              OCCGeometry::GetProperties(ef.Current()).Merge(prop);
          }
      }
    TopoDS_Shape res = splitter.Shape();
    const gp_Ax1 ax = axes.Axis();
    const gp_Dir n = ax.Direction();
    const gp_Pnt p0 = axes.Location();
    Array<TopoDS_Shape> out_shapes;
    for (auto f : GetFaces(res))
      {
        Handle(Geom_Surface) s = BRep_Tool::Surface(TopoDS::Face(f));
        Handle(Geom_Plane) gp = Handle(Geom_Plane)::DownCast(s);
        if (gp.IsNull())
          continue;

        // Check coplanarity (orientation can flip)
        gp_Pln fpln = gp->Pln();
        auto my_n = fpln.Position().Axis().Direction();
        if (!my_n.IsParallel(n, 1e-10))
          continue;

        // Check distance of face plane to our plane
        Standard_Real dist = fpln.Distance(p0);
        if (dist > 1e-7 * diam)
          continue;
        if (my_n.Dot(n) > 0)
          f.Reverse();
        out_shapes.Append(f);
      }
    BRepBuilderAPI_Sewing sew(1e-7 * diam);
    for (auto s : out_shapes)
      sew.Add(TopoDS::Face(s));
    sew.Perform();
    for (auto s : out_shapes)
      PropagateProperties(sew, s);
    auto out = sew.SewedShape();
    gp_Ax3 xy(gp::XOY());

    gp_Trsf T;
    T.SetTransformation(axes, xy);
    BRepBuilderAPI_Transform btransf(out, T, Standard_True);
    PropagateProperties(btransf, out);
    return btransf.Shape();
  }
} // namespace netgen
