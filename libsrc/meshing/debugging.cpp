#include <meshing.hpp>

namespace netgen
{
    unique_ptr<Mesh> GetOpenElements( const Mesh & m, int dom = 0 )
    {
        static Timer t("GetOpenElements"); RegionTimer rt(t);
        auto mesh = make_unique<Mesh>();
        *mesh = m;

        Array<bool, PointIndex> interesting_points(mesh->GetNP());
        interesting_points = false;

        mesh->FindOpenElements(dom);
        NgArray<Element2d> openelements;
        openelements = mesh->OpenElements();

        for (auto & el : openelements)
            for (auto i : el.PNums())
                interesting_points[i] = true;

        for (auto & el : mesh->VolumeElements())
        {
            int num_interesting_points = 0;

            for (auto pi : el.PNums())
                if(interesting_points[pi])
                    num_interesting_points++;

            if(num_interesting_points==0)
                el.Delete();
            el.SetIndex(num_interesting_points);
        }

        mesh->SetMaterial(1, "1_point");
        mesh->SetMaterial(2, "2_points");
        mesh->SetMaterial(3, "3_points");
        mesh->SetMaterial(4, "4_points");

        mesh->ClearSurfaceElements();

        for (auto & el : openelements)
            mesh->AddSurfaceElement( el );

        mesh->Compress();
        return mesh;
    }

    unique_ptr<Mesh> FilterMesh( const Mesh & m, FlatArray<PointIndex> points, FlatArray<SurfaceElementIndex> sels, FlatArray<ElementIndex> els )
    {
        static Timer t("GetOpenElements"); RegionTimer rt(t);
        auto mesh_ptr = make_unique<Mesh>();
        auto & mesh = *mesh_ptr;
        mesh = m;

        Array<bool, PointIndex> keep_point(mesh.GetNP());
        Array<bool, SurfaceElementIndex> keep_sel(mesh.GetNSE());
        Array<bool, ElementIndex> keep_el(mesh.GetNE());
        mesh.LineSegments().DeleteAll();

        keep_point = false;
        for(auto pi : points)
            keep_point[pi] = true;

        auto set_keep = [&] (auto & input, auto & keep_array, auto & els)
        {
            keep_array = false;
            for(auto ind : input)
                keep_array[ind] = true;

            for(auto ind : Range(els))
            {
                bool & keep = keep_array[ind];
                if(keep) continue;

                for(auto pi : mesh[ind].PNums())
                    keep |= keep_point[pi];

                if(!keep)
                    mesh[ind].Delete();
            }

            for(auto i = 0; i<els.Size(); i++)
                if(els[i].IsDeleted())
                {
                    els.DeleteElement(i);
                    i--;
                }
        };

        set_keep(sels, keep_sel, mesh.SurfaceElements());
        set_keep(els, keep_el, mesh.VolumeElements());
        //mesh.Compress();

        return mesh_ptr;
    }



}
