#include <mystdlib.h>
#include "geomsearch.hpp"
#include "adfront3.hpp"


namespace netgen
{
  GeomSearch3d :: GeomSearch3d() 
  {
    size[0] = 0; size[1] = 0; size[2] = 0; 
  };

  GeomSearch3d :: ~GeomSearch3d()
  {
    //delete old Hashtable:
    if (size[0] != 0)
      {
        for (int i = 0; i < size[0]*size[1]*size[2]; i++)
          delete hashtable[i];
      } 
  }

  void GeomSearch3d :: Init (Array <FrontPoint3,Front3PointIndex> *pointsi, Array <FrontFace> *facesi)
  {
    points = pointsi;
    faces = facesi;
    size[0] = 0; size[1] = 0; size[2] = 0; 
    reset = 1;
    hashcount = 1;
  }

  void GeomSearch3d :: ElemMaxExt(Point<3>& minp, Point<3>& maxp, const FrontElement2d& elem)
  {
    maxp(0)=(*points)[elem.PNum(1)].P()(0);
    maxp(1)=(*points)[elem.PNum(1)].P()(1);
    maxp(2)=(*points)[elem.PNum(1)].P()(2);
    minp(0)=(*points)[elem.PNum(1)].P()(0);
    minp(1)=(*points)[elem.PNum(1)].P()(1);
    minp(2)=(*points)[elem.PNum(1)].P()(2);
  
    for (int i=2; i <= 3; i++)
      {
        maxp(0)=max2((*points)[elem.PNum(i)].P()(0),maxp(0));
        maxp(1)=max2((*points)[elem.PNum(i)].P()(1),maxp(1));
        maxp(2)=max2((*points)[elem.PNum(i)].P()(2),maxp(2));
        minp(0)=min2((*points)[elem.PNum(i)].P()(0),minp(0));
        minp(1)=min2((*points)[elem.PNum(i)].P()(1),minp(1));
        minp(2)=min2((*points)[elem.PNum(i)].P()(2),minp(2));
      }
  }

  void GeomSearch3d :: MinCoords(const Point<3>& p1, Point<3>& p2)
  {
    p2(0)=min2(p1(0),p2(0));
    p2(1)=min2(p1(1),p2(1));
    p2(2)=min2(p1(2),p2(2));
  }

  void GeomSearch3d :: MaxCoords(const Point<3>& p1, Point<3>& p2)
  {
    p2(0)=max2(p1(0),p2(0));
    p2(1)=max2(p1(1),p2(1));
    p2(2)=max2(p1(2),p2(2));
  }

  void GeomSearch3d :: Create()
  {
    INDEX i,j,k;
    if (reset)
      {
        const double hashelemsizefactor = 4;
        reset = 0;
        /*
          minext=Point<3>(MAXDOUBLE, MAXDOUBLE, MAXDOUBLE);
          maxext=Point<3>(MINDOUBLE, MINDOUBLE, MINDOUBLE);
        */
        ElemMaxExt(minext, maxext, faces->operator[](0).Face());
        Point<3> maxp, minp;
        Vec<3> midext(0,0,0);
      
        //get max Extension of Frontfaces
        for (i = 1; i <= faces->Size(); i++)
          {
            ElemMaxExt(minp, maxp, faces->operator[](i-1).Face());
            MinCoords(minp, minext);
            MaxCoords(maxp, maxext);
            midext+=maxp-minp;
          }


        maxextreal = maxext;
        maxext = maxext + 1e-4 * (maxext - minext);

        midext*=1./faces->Size();
        Vec<3> boxext = maxext - minext;
      
        //delete old Hashtable:
        if (size[0] != 0)
          {
            for (i = 1; i <= size[0]*size[1]*size[2]; i++)
              {
                delete hashtable[i-1];
              }
          } 
      
        size[0] = int (boxext(0)/midext(0)/hashelemsizefactor+1);
        size[1] = int (boxext(1)/midext(1)/hashelemsizefactor+1);
        size[2] = int (boxext(2)/midext(2)/hashelemsizefactor+1);

        int nfaces = faces->Size();
        size[0] = min(size[0], nfaces);
        size[1] = min(size[1], nfaces);
        size[2] = min(size[2], nfaces);

        // PrintMessage (5, "hashsizes = ", size[0], ", ", size[1], ", ", size[2]);
      
        elemsize(0)=boxext(0)/size[0];
        elemsize(1)=boxext(1)/size[1];
        elemsize(2)=boxext(2)/size[2];

        //create Hasharrays:
        hashtable.SetSize(size[0]*size[1]*size[2]);
        for (i = 1; i <= size[0]; i++)
          {
            for (j = 1; j <= size[1]; j++)
              {
                for (k = 1; k <= size[2]; k++)
                  {
                    INDEX ind=i+(j-1)*size[0]+(k-1)*size[1]*size[0];
                    hashtable[ind-1] = new Array <int> ();
                  }
              }
          }
      }
    else
      {
        //Clear all Hash-Arrays
        for (i = 1; i <= size[0]; i++)
          {
            for (j = 1; j <= size[1]; j++)
              {
                for (k = 1; k <= size[2]; k++)
                  {
                    INDEX ind=i+(j-1)*size[0]+(k-1)*size[1]*size[0];
                    hashtable[ind-1]->SetSize(0);
                  }
              }
          }       
      }
  
    //Faces in Hashtable einfuegen:
    for (i = 1; i <= faces->Size(); i++)
      {
        AddElem(faces->operator[](i-1).Face(),i);
      }
  
  }

  void GeomSearch3d :: AddElem(const FrontElement2d& elem, INDEX elemnum)
  {
    Point<3> minp, maxp;
    ElemMaxExt(minp, maxp, elem);
    int sx = int ((minp(0)-minext(0))/elemsize(0)+1.);
    int ex = int ((maxp(0)-minext(0))/elemsize(0)+1.);
    int sy = int ((minp(1)-minext(1))/elemsize(1)+1.);
    int ey = int ((maxp(1)-minext(1))/elemsize(1)+1.);
    int sz = int ((minp(2)-minext(2))/elemsize(2)+1.);
    int ez = int ((maxp(2)-minext(2))/elemsize(2)+1.);
  
    for (int ix = sx; ix <= ex; ix++)
      for (int iy = sy; iy <= ey; iy++)
        for (int iz = sz; iz <= ez; iz++)
          {
            INDEX ind=ix+(iy-1)*size[0]+(iz-1)*size[1]*size[0];
            if (ind < 1 || ind > size[0] * size[1] * size[2])
              {
                cerr << "Illegal hash-position";
                cerr << "Position: " << ix << "," << iy << "," << iz << endl;
                    throw NgException ("Illegal position in Geomsearch");
              }
            hashtable[ind-1]->Append(elemnum);                
          }
  }

  void GeomSearch3d :: GetLocals(Array<FrontElement2d> & locfaces,  Array<INDEX> & findex,
                                 INDEX fstind, const Point<3>& p0, double xh)
  {
    hashcount++;
  
    Point<3> minp, maxp, midp; 

    minp=p0-Vec<3>(xh,xh,xh); //lay cube over sphere
    maxp=p0+Vec<3>(xh,xh,xh);

    MaxCoords(minext,minp); //cube may not be out of hash-region
    MinCoords(maxextreal,maxp);


    Front3PointIndex cluster = faces->operator[](fstind-1).Cluster();
  
    int sx = int((minp(0)-minext(0))/elemsize(0)+1.);
    int ex = int((maxp(0)-minext(0))/elemsize(0)+1.);
    int sy = int((minp(1)-minext(1))/elemsize(1)+1.);
    int ey = int((maxp(1)-minext(1))/elemsize(1)+1.);
    int sz = int((minp(2)-minext(2))/elemsize(2)+1.);
    int ez = int((maxp(2)-minext(2))/elemsize(2)+1.);
    int ix,iy,iz,i,k;

    [[maybe_unused]] int cnt1 = 0;  // test, how efficient hashtable is
    [[maybe_unused]] int cnt2 = 0;
    [[maybe_unused]] int cnt3 = 0;
  
    for (ix = sx; ix <= ex; ix++)
      {
        for (iy = sy; iy <= ey; iy++)
          {
            for (iz = sz; iz <= ez; iz++)
              {
                INDEX ind=ix+(iy-1)*size[0]+(iz-1)*size[1]*size[0];
              
                //go through all elements in one hash area
                const Array <int> & area = *hashtable[ind-1];
                for (k = 1; k <= area.Size(); k++)
                  {
                    cnt2++;
                    i = area[k-1];
                    if (faces->operator[](i-1).Cluster() == cluster && 
                        faces->operator[](i-1).Valid() &&
                        faces->operator[](i-1).HashValue() != hashcount && 
                        i != fstind)
                      {
                        cnt1++;
                        const FrontElement2d & face = faces->operator[](i-1).Face();
                      
                        const Point<3> & p1 = (*points)[face.PNum(1)].P();
                        const Point<3> & p2 = (*points)[face.PNum(2)].P();
                        const Point<3> & p3 = (*points)[face.PNum(3)].P();
                      
                        midp = Center (p1, p2, p3);
                      
                        // if (Dist2 (midp, p0) <= xh*xh)  
                        if((Dist2 (p1, p0) <= xh*xh) ||
                           (Dist2 (p2, p0) <= xh*xh) ||
                           (Dist2 (p3, p0) <= xh*xh) ||
                           (Dist2 (midp, p0) <= xh*xh) )  // by Jochen Wild
                          {
                            cnt3++;
                            locfaces.Append(faces->operator[](i-1).Face());
                            findex.Append(i);
                            faces->operator[](i-1).SetHashValue(hashcount);
                          }
                      }
                  }
              }
          }
      }
    /*
      if (faces->Size() != 0 && hashcount % 200 == 0)
      {
      (*mycout) << "n.o.f= " << faces->Size();
      (*mycout) << ", n.o.lf= " << locfaces.Size();
      (*mycout) << ", hashf= " << (double)cnt2/(double)faces->Size();
      (*mycout) << " (" << (double)cnt1/(double)faces->Size();
      (*mycout) << ", " << (double)cnt3/(double)faces->Size() << ")" << endl;
      }
    */

  }

}
