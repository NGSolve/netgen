#include <meshclass.hpp>


namespace netgen
{
    inline unique_ptr<Mesh> GetOpenElements( const Mesh & m, int dom = 0 )
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
        mesh->Compress();

        mesh->ClearSurfaceElements();

        for (auto & el : openelements)
            mesh->AddSurfaceElement( el );

        return mesh;
    }


}
