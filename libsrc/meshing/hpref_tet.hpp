






enum VNUM { V1, V2, V3, V4,
            E12, E13, E14,
            E21, E23, E24,
            E31, E32, E34,
            E41, E42, E43,
            F123, F124, F134,
            F213, F214, F234,
            F312, F314, F324,
            F412, F413, F423
};

class El
{
public:
  HPREF_ELEMENT_TYPE type;
  std::vector<VNUM> vertices;

  El (HPREF_ELEMENT_TYPE atype,
      std::vector<VNUM> avertices)
    : type(atype), vertices(avertices) { } 
};

extern std::map<HPREF_ELEMENT_TYPE, HPRef_Struct*> & GetHPRegistry();

template <HPREF_ELEMENT_TYPE GEOM>
class HPRefStruct : public HPRef_Struct
{
  typedef int int3[3];
  typedef int int4[4];
  typedef int int8[8];  
  std::vector<std::array<int,3>> refedges;
  std::vector<std::array<int,4>> reffaces;
  std::vector<HPREF_ELEMENT_TYPE> neweltypes_vec;
  std::vector<std::array<int,8>> newelverts;

public:
  HPRefStruct(HPREF_ELEMENT_TYPE type,
              std::vector<El> list)
  {
    GetHPRegistry()[type] = this;
    
    geom = GEOM;
    std::map<VNUM, int> mapnums;
    int ii = 0;
    for (auto v : { V1, V2, V3, V4})
      mapnums[v] = ++ii;
    
    for (auto el : list)
      for (auto v : el.vertices)
        if (mapnums.count(v)==0)
          mapnums[v] = ++ii;

    int elist[][3] =
      { { 1, 2, E12 },
        { 1, 3, E13 },
        { 1, 4, E14 },
        { 2, 1, E21 },
        { 2, 3, E23 },
        { 2, 4, E24 },
        { 3, 1, E31 },
        { 3, 2, E32 },
        { 3, 4, E34 },
        { 4, 1, E41 },
        { 4, 2, E42 },
        { 4, 3, E43 }
      };
    int flist[][4] =
      { { 1, 2, 3, F123 },
        { 1, 2, 4, F124 },
        { 1, 3, 4, F134 },
        { 2, 1, 3, F213 },
        { 2, 1, 4, F214 },
        { 2, 3, 4, F234 },
        { 3, 1, 2, F312 },
        { 3, 1, 4, F314 },
        { 3, 2, 4, F324 },
        { 4, 1, 2, F412 },
        { 4, 1, 3, F413 },
        { 4, 2, 3, F423 }
      };
      
    /*
      // too advanced ...
    for (auto [i1,i2,inew] : elist)
      if (mapnums.count(VNUM(inew)))
        refedges.push_back( { i1, i2, mapnums[VNUM(inew)] });
    */
    for (int i = 0; i < size(elist); i++)
      {
        int i1 = elist[i][0];
        int i2 = elist[i][1];
        int inew = elist[i][2];
        if (mapnums.count(VNUM(inew)))
          refedges.push_back( { i1, i2, mapnums[VNUM(inew)] });
      }
    
    refedges.push_back( { 0, 0, 0 } );
    splitedges = (int3*) &refedges[0][0];

    /*
      // too advanced ...
    for (auto [i1,i2,i3,inew] : flist)
      if (mapnums.count(VNUM(inew)))
        reffaces.push_back( { i1, i2, i3, mapnums[VNUM(inew)] });
    */
    for (int i = 0; i < size(flist); i++)
      {
        int i1 = flist[i][0];
        int i2 = flist[i][1];
        int i3 = flist[i][2];
        int inew = flist[i][3];
        if (mapnums.count(VNUM(inew)))
          reffaces.push_back( { i1, i2, i3, mapnums[VNUM(inew)] });
      }



    reffaces.push_back( { 0, 0, 0 } );
    splitfaces = (int4*) &reffaces[0][0];
    


    splitelements = nullptr;
    
    for (auto el : list)
      {
        neweltypes_vec.push_back (el.type);
        std::array<int,8> verts;
        for (int j = 0; j < std::min(verts.size(), el.vertices.size()); j++)
          verts[j] = mapnums[VNUM(el.vertices[j])];
        newelverts.push_back(verts);
      }
    
    neweltypes_vec.push_back (HP_NONE);

    neweltypes = &neweltypes_vec[0];
    newels = (int8*) &newelverts[0][0];    

    /*
    int ind = 0;
    cout << "rule, split edges:" << endl;
    while (splitedges[ind][0])
      {
        cout << splitedges[ind][0] << "-" << splitedges[ind][1] << ": " << splitedges[ind][2] << endl;
        ind++;
      }
    
    ind = 0;
    cout << "rule, split faces:" << endl;
    while (splitfaces[ind][0])
      {
        cout << splitfaces[ind][0] << "-" << splitfaces[ind][1]
             << "-" << splitfaces[ind][2] << ": " << splitfaces[ind][3] << endl;
        ind++;
      }

    ind = 0;
    cout << "rule, new els:" << endl;
    while (neweltypes[ind] != HP_NONE)
      {
        cout << "new type " << neweltypes[ind] << ", verts: ";
        for (int j = 0; j < 8; j++)
          cout << newels[ind][j] << " ";
        ind++;
      }
    */
  }
};





 
// HP_NONETET
int refnonetet_splitedges[][3] =
{
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE refnonetet_newelstypes[] =
{
  HP_TET,  
  HP_NONE,
};
int refnonetet_newels[][8] =
{
  { 1, 1, 1, 1 },
};
HPRef_Struct refnonetet =
{
  HP_TET,
  refnonetet_splitedges, 
  0, 0,
  refnonetet_newelstypes, 
  refnonetet_newels
};





// HP_TET
int reftet_splitedges[][3] =
{
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftet_newelstypes[] =
{
  HP_TET,
  HP_NONE,
};
int reftet_newels[][8] =
{
  { 1, 2, 3, 4 },
};
HPRef_Struct reftet =
{
  HP_TET,
  reftet_splitedges, 
  0, 0,
  reftet_newelstypes, 
  reftet_newels
};



/* *********** Tet - Refinement - 0 edges *************** */

// HP_TET_0E_1V
int reftet_0e_1v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftet_0e_1v_newelstypes[] =
{
  HP_TET_0E_1V,
  HP_PRISM,
  HP_NONE,
};
int reftet_0e_1v_newels[][8] =
{
  { 1, 5, 6, 7 },
  { 5, 6, 7, 2, 3, 4 }
};
HPRef_Struct reftet_0e_1v =
{
  HP_TET,
  reftet_0e_1v_splitedges, 
  0, 0,
  reftet_0e_1v_newelstypes, 
  reftet_0e_1v_newels
};

/*
  // new syntax ???
HPRef_Struct2 str =
  {
    HP_TET_0E_1V, HP_TET,
    El(HP_TET_0E_1V, { V1, V12, V13, V14 })
  };
*/  




// HP_TET_0E_2V
int reftet_0e_2v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftet_0e_2v_newelstypes[] =
{
  HP_TET_0E_1V,
  HP_TET_0E_1V,
  HP_PRISM,
  HP_PRISM,
  HP_NONE,
};
int reftet_0e_2v_newels[][8] =
{
  { 1, 5, 6, 7 },
  { 2, 10, 9, 8 },
  { 5, 6, 7, 8, 9, 10 },
  { 4, 10, 7, 3, 9, 6 },
};
HPRef_Struct reftet_0e_2v =
{
  HP_TET,
  reftet_0e_2v_splitedges, 
  0, 0,
  reftet_0e_2v_newelstypes, 
  reftet_0e_2v_newels
};





// HP_TET_0E_3V
int reftet_0e_3v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 0, 0, 0 }
};
int reftet_0e_3v_splitfaces[][4] =
  {
    { 1, 2, 3, 14 },
    { 2, 3, 1, 15 },
    { 3, 1, 2, 16 },
    { 0, 0, 0, 0 },
  };
HPREF_ELEMENT_TYPE reftet_0e_3v_newelstypes[] =
{
  HP_PYRAMID_0E_1V,
  HP_PYRAMID_0E_1V,
  HP_PYRAMID_0E_1V,
  HP_PRISM,
  HP_PRISM,
  HP_PRISM,
  HP_PRISM,
  HP_TET,
  HP_NONE,
};
int reftet_0e_3v_newels[][8] =
{
  { 1, 5, 14, 6, 7 },
  { 2, 9, 15, 8, 10 },
  { 3, 11, 16, 12, 13 },
  { 5, 14, 7, 8, 15, 10 },
  { 9, 15, 10, 12, 16, 13 },
  { 6, 7, 14, 11, 13, 16 },
  { 14, 15, 16, 7, 10, 13 },
  { 7, 10, 13, 4 }
};
HPRef_Struct reftet_0e_3v =
{
  HP_TET,
  reftet_0e_3v_splitedges, 
  reftet_0e_3v_splitfaces, 
  0,
  reftet_0e_3v_newelstypes, 
  reftet_0e_3v_newels
};





// HP_TET_0E_4V
int reftet_0e_4v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_0e_4v_splitfaces[][4] =
  {
    { 1, 2, 3, 17 },
    { 1, 2, 4, 18 },
    { 1, 3, 4, 19 },

    { 2, 1, 3, 20 },
    { 2, 1, 4, 21 },
    { 2, 3, 4, 22 },

    { 3, 1, 2, 23 },
    { 3, 1, 4, 24 },
    { 3, 2, 4, 25 },

    { 4, 1, 2, 26 },
    { 4, 1, 3, 27 },
    { 4, 2, 3, 28 },
    { 0, 0, 0, 0 },
  };
int reftet_0e_4v_splitelements[][5] =
  {
    { 1, 2, 3, 4, 29 },
    { 2, 3, 4, 1, 30 },
    { 3, 4, 1, 2, 31 },
    { 4, 1, 2, 3, 32 },
    { 0 },
  };
HPREF_ELEMENT_TYPE reftet_0e_4v_newelstypes[] =
{
  HP_HEX_0E_1V,
  HP_HEX_0E_1V,
  HP_HEX_0E_1V,
  HP_HEX_0E_1V,
  HP_PRISM, HP_PRISM, 
  HP_PRISM, HP_PRISM, 
  HP_PRISM, HP_PRISM, 
  HP_PRISM, HP_PRISM, 
  HP_PRISM, HP_PRISM, 
  HP_PRISM, HP_PRISM, 
  HP_PRISM,
  HP_PRISM,
  HP_PRISM,
  HP_PRISM,
  HP_TET,
  HP_NONE,
};
int reftet_0e_4v_newels[][8] =
{
  { 1, 5, 17, 6, 7, 18, 29, 19 },
  { 2, 9, 20, 8, 10, 22, 30, 21 },
  { 3, 11, 23, 12, 13, 24, 31, 25 },
  { 4, 15, 26, 14, 16, 28, 32, 27 },
  { 5, 17, 18, 8, 20, 21 },
  { 18, 17, 29, 21, 20, 30 },
  { 6, 19, 17,  11, 24, 23 },
  { 17, 19, 29,  23, 24, 31 },
  { 7, 18, 19, 14, 26, 27 },
  { 19, 18, 29, 27, 26, 32 },
  { 9, 20, 22, 12, 23, 25 },
  { 22, 20, 30, 25, 23, 31 },
  { 10, 22, 21, 15, 28, 26 },
  { 21, 22, 30, 26, 28, 32 },
  { 13, 24, 25, 16, 27, 28 },
  { 25, 24, 31, 28, 27, 32 },
  { 17, 20, 23, 29, 30, 31 },
  { 18, 26, 21, 29, 32, 30 },
  { 19, 24, 27, 29, 31, 32 },
  { 22, 28, 25, 30, 32, 31 },
  { 29, 30, 31, 32 },
};
HPRef_Struct reftet_0e_4v =
{
  HP_TET,
  reftet_0e_4v_splitedges, 
  reftet_0e_4v_splitfaces, 
  reftet_0e_4v_splitelements, 
  reftet_0e_4v_newelstypes, 
  reftet_0e_4v_newels
};

















/* *********** Tet - Refinement - 1 edge *************** */



// HP_TET_1E_0V
int reftet_1e_0v_splitedges[][3] =
{
  { 1, 3, 5 },
  { 1, 4, 6 },
  { 2, 3, 7 },
  { 2, 4, 8 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftet_1e_0v_newelstypes[] =
{
  HP_PRISM_SINGEDGE,
  HP_PRISM,
  HP_NONE,
};
int reftet_1e_0v_newels[][8] =
{
  { 1, 5, 6, 2, 7, 8 },
  { 7, 3, 5, 8, 4, 6 }
};
HPRef_Struct reftet_1e_0v =
{
  HP_TET,
  reftet_1e_0v_splitedges, 
  0, 0,
  reftet_1e_0v_newelstypes, 
  reftet_1e_0v_newels
};





// HP_TET_1E_1VA
int reftet_1e_1va_splitedges[][3] =
{
  { 1, 3, 5 },
  { 1, 4, 6 },
  { 2, 3, 7 },
  { 2, 4, 8 },
  { 1, 2, 9 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftet_1e_1va_newelstypes[] =
{
  HP_TET_1E_1VA,
  HP_PRISM_SINGEDGE,
  HP_PRISM,
  HP_NONE,
};
int reftet_1e_1va_newels[][8] =
{
  { 1, 9, 5, 6 },
  { 9, 5, 6, 2, 7, 8 },
  { 7, 3, 5, 8, 4, 6 }
};
HPRef_Struct reftet_1e_1va =
{
  HP_TET,
  reftet_1e_1va_splitedges, 
  0, 0,
  reftet_1e_1va_newelstypes, 
  reftet_1e_1va_newels
};






// HP_TET_1E_1VB
int reftet_1e_1vb_splitedges[][3] =
{
  { 1, 3, 5 },
  { 1, 4, 6 },
  { 2, 3, 7 },
  { 2, 4, 8 },
  { 4, 1, 9 },
  { 4, 2, 10 },
  { 4, 3, 11 },
  { 0, 0, 0 }
};
int reftet_1e_1vb_splitelements[][5] =
{
  { 4, 1, 2, 3, 12 },
  { 0 }
};

HPREF_ELEMENT_TYPE reftet_1e_1vb_newelstypes[] =
{
  HP_PRISM_SINGEDGE,
  HP_TET_0E_1V,
  HP_PYRAMID,
  HP_TET,
  HP_PYRAMID, 
  HP_TET,
  HP_PYRAMID,
  HP_TET,
  HP_PYRAMID,
  HP_TET,
  HP_NONE,
};
int reftet_1e_1vb_newels[][8] =
{
  { 1, 5, 6, 2, 7, 8 },
  { 4, 11, 10, 9 },
  { 7, 8, 10, 11, 12 },
  { 3, 7, 11, 12 },
  { 5, 11, 9, 6, 12 },
  { 5, 3, 11, 12 },
  { 6, 9, 10, 8, 12 },
  { 5, 7, 3, 12 },
  { 5, 6, 8, 7, 12 },
  { 9, 11, 10, 12 }
};
HPRef_Struct reftet_1e_1vb =
{
  HP_TET,
  reftet_1e_1vb_splitedges, 
  0,
  reftet_1e_1vb_splitelements, 
  reftet_1e_1vb_newelstypes, 
  reftet_1e_1vb_newels
};








// HP_TET_1E_2VA
int reftet_1e_2va_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftet_1e_2va_newelstypes[] =
{
  HP_TET_1E_1VA,
  HP_TET_1E_1VA,
  HP_PRISM_SINGEDGE,
  HP_PRISM,
  HP_NONE,
};
int reftet_1e_2va_newels[][8] =
{
  { 1, 5, 6, 7 },
  { 2, 8, 10, 9 },
  { 5, 6, 7, 8, 9, 10 },
  { 4, 10, 7, 3, 9, 6 },
};
HPRef_Struct reftet_1e_2va =
{
  HP_TET,
  reftet_1e_2va_splitedges, 
  0, 0,
  reftet_1e_2va_newelstypes, 
  reftet_1e_2va_newels
};







// HP_TET_1E_2VB
int reftet_1e_2vb_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 3, 8 },
  { 2, 4, 9 },
  { 3, 1, 10 },
  { 3, 2, 11 },
  { 3, 4, 12 },
  { 0, 0, 0 }
};
int reftet_1e_2vb_splitelements[][5] =
{
  { 3, 4, 1, 2, 13 },
  { 0 }
};

HPREF_ELEMENT_TYPE reftet_1e_2vb_newelstypes[] =
{
  HP_TET_1E_1VA,
  HP_PRISM_SINGEDGE,
  HP_TET_0E_1V,
  HP_PYRAMID,
  HP_TET,
  HP_PYRAMID, 
  HP_TET,
  HP_PYRAMID,
  HP_TET,
  HP_PYRAMID,
  HP_TET,
  HP_NONE,
};
int reftet_1e_2vb_newels[][8] =
{
  { 1, 5, 6, 7 },
  { 5, 6, 7, 2, 8, 9 },
  { 3, 10, 11, 12 },

  { 8, 9, 12, 11, 13 },
  { 4, 12, 9, 13 },
  { 6, 10, 12, 7, 13 },
  { 4, 7, 12, 13 },
  { 6, 8, 11, 10, 13 },
  { 4, 9, 7, 13 },
  { 6, 7, 9, 8, 13 },
  { 10, 11, 12, 13 },
};
HPRef_Struct reftet_1e_2vb =
{
  HP_TET,
  reftet_1e_2vb_splitedges, 
  0,
  reftet_1e_2vb_splitelements, 
  reftet_1e_2vb_newelstypes, 
  reftet_1e_2vb_newels
};






// HP_TET_1E_2VC
int reftet_1e_2vc_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 3, 8 },
  { 2, 4, 9 },
  { 4, 1, 10 },
  { 4, 2, 11 },
  { 4, 3, 12 },
  { 0, 0, 0 }
};
int reftet_1e_2vc_splitelements[][5] =
{
  { 4, 1, 2, 3, 13 },
  { 0 }
};

HPREF_ELEMENT_TYPE reftet_1e_2vc_newelstypes[] =
{
  HP_TET_1E_1VA,
  HP_PRISM_SINGEDGE,
  HP_TET_0E_1V,
  HP_PYRAMID,
  HP_TET,
  HP_PYRAMID, 
  HP_TET,
  HP_PYRAMID,
  HP_TET,
  HP_PYRAMID,
  HP_TET,
  HP_NONE,
};
int reftet_1e_2vc_newels[][8] =
{
  { 1, 5, 6, 7 },
  { 5, 6, 7, 2, 8, 9 },
  { 4, 11, 10, 12 },
  { 8, 9, 11, 12, 13 },
  { 3, 8, 12, 13 },
  { 7, 6, 12, 10, 13 },
  { 3, 12, 6, 13 },
  { 9, 7, 10, 11, 13 },
  { 3, 6, 8, 13 },
  { 6, 7, 9, 8, 13 },
  { 10, 12, 11, 13 }
};
HPRef_Struct reftet_1e_2vc =
{
  HP_TET,
  reftet_1e_2vc_splitedges, 
  0,
  reftet_1e_2vc_splitelements, 
  reftet_1e_2vc_newelstypes, 
  reftet_1e_2vc_newels
};








/*

// HP_TET_1E_2VD
int reftet_1e_2vd_splitedges[][3] =
{
  { 1, 3, 5 },
  { 1, 4, 6 },
  { 2, 3, 7 },
  { 2, 4, 8 },
  { 3, 1, 9 },
  { 3, 2, 10 },
  { 3, 4, 11 },
  { 4, 1, 12 },
  { 4, 2, 13 },
  { 4, 3, 14 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftet_1e_2vd_newelstypes[] =
{
  HP_PRISM_SINGEDGE,
  HP_TET_0E_1V,
  HP_TET_0E_1V,
  HP_PRISM,
  HP_HEX,
  HP_NONE,
};
int reftet_1e_2vd_newels[][8] =
{
  { 1, 5, 6, 2, 7, 8 },
  { 4, 13, 12, 14 },
  { 3, 10, 11, 9 },
  { 14, 13, 12, 11, 10, 9 },
  { 6, 12, 13, 8, 5, 9, 10, 7 },
};
HPRef_Struct reftet_1e_2vd =
{
  HP_TET,
  reftet_1e_2vd_splitedges, 
  0, 0,
  reftet_1e_2vd_newelstypes, 
  reftet_1e_2vd_newels
};

*/




//  HP_TET_1E_2VD,  // 1 v on edge
int reftet_1e_2vd_splitedges[][3] =
{
  // { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  // { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_1e_2vd_splitfaces[][4] =
  {
    { 1, 3, 4, 19 },
    { 2, 3, 4, 22 },
    { 3, 1, 4, 24 },
    { 3, 2, 4, 25 },
    { 4, 1, 3, 27 },
    { 4, 2, 3, 28 },
    { 0, 0, 0, 0 }
  };
HPREF_ELEMENT_TYPE reftet_1e_2vd_newelstypes[] =
  {
    HP_PRISM_SINGEDGE,
    HP_TET_0E_1V,
    HP_TET_0E_1V,
    HP_PRISM,
    HP_HEX,
    HP_PYRAMID,
    HP_HEX,
    HP_PYRAMID,
    HP_PRISM,
    HP_PRISM,
    HP_NONE,
  };
int reftet_1e_2vd_newels[][8] =
{
  { 1, 6, 7, 2, 9, 10 },
  { 3, 11, 12, 13 },
  { 4, 16, 15, 14 },
  { 7, 6, 19, 10, 9, 22 },
  { 7, 19, 27, 14, 10, 22, 28, 15 },
  { 14, 15, 28, 27, 16 },
  { 9, 6, 19, 22, 12, 11, 24, 25 },
  { 12, 11, 24, 25, 13 },
  { 19, 24, 27, 22, 25, 28 },
  { 16, 28, 27, 13, 25, 24 }
};
HPRef_Struct reftet_1e_2vd =
{
  HP_TET,
  reftet_1e_2vd_splitedges, 
  reftet_1e_2vd_splitfaces, 
  0,
  reftet_1e_2vd_newelstypes, 
  reftet_1e_2vd_newels
};















// HP_TET_1E_3VA
int reftet_1e_3va_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 0, 0, 0 }
};
int reftet_1e_3va_splitelements[][5] =
{
  { 1, 2, 3, 4, 14 },
  { 0 }
};

HPREF_ELEMENT_TYPE reftet_1e_3va_newelstypes[] =
{
  HP_PRISM_SINGEDGE,
  HP_TET_1E_1VA,
  HP_TET_1E_1VA,
  HP_TET_0E_1V,

  HP_PYRAMID,
  HP_TET,
  HP_PYRAMID, 
  HP_TET,
  HP_PYRAMID,
  HP_TET,
  HP_PYRAMID,
  HP_TET,
  HP_NONE,
};
int reftet_1e_3va_newels[][8] =
{
  { 5, 6, 7, 8, 9, 10 },
  { 1, 5, 6, 7 },
  { 2, 8, 10, 9 },
  { 3, 11, 12, 13 },

  { 6, 7, 10, 9, 14 },
  { 4, 10, 7, 14 },
  { 9, 10, 13, 12, 14 },
  { 4, 13, 10, 14 },
  { 6, 11, 13, 7, 14 },
  { 4, 7, 13, 14 },
  { 6, 11, 12, 9, 14 },
  { 11, 13, 12, 14 },
};

HPRef_Struct reftet_1e_3va =
{
  HP_TET,
  reftet_1e_3va_splitedges, 
  0,
  reftet_1e_3va_splitelements, 
  reftet_1e_3va_newelstypes, 
  reftet_1e_3va_newels
};






















//  HP_TET_1E_3VB,  // 1 v on edge
int reftet_1e_3vb_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  // { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_1e_3vb_splitfaces[][4] =
  {
    { 1, 3, 4, 19 },
    { 2, 3, 4, 22 },
    { 3, 1, 4, 24 },
    { 3, 2, 4, 25 },
    { 4, 1, 3, 27 },
    { 4, 2, 3, 28 },
    { 0, 0, 0, 0 }
  };
HPREF_ELEMENT_TYPE reftet_1e_3vb_newelstypes[] =
  {
    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_TET_0E_1V,
    HP_TET_0E_1V,
    HP_PRISM,
    HP_HEX,
    HP_PYRAMID,
    HP_HEX,
    HP_PYRAMID,
    HP_PRISM,
    HP_PRISM,
    HP_NONE,
  };
int reftet_1e_3vb_newels[][8] =
{
  { 1, 5, 6, 7 },
  { 5, 6, 7, 2, 9, 10 },
  { 3, 11, 12, 13 },
  { 4, 16, 15, 14 },
  { 7, 6, 19, 10, 9, 22 },
  { 7, 19, 27, 14, 10, 22, 28, 15 },
  { 14, 15, 28, 27, 16 },
  { 9, 6, 19, 22, 12, 11, 24, 25 },
  { 12, 11, 24, 25, 13 },
  { 19, 24, 27, 22, 25, 28 },
  { 16, 28, 27, 13, 25, 24 }
};
HPRef_Struct reftet_1e_3vb =
{
  HP_TET,
  reftet_1e_3vb_splitedges, 
  reftet_1e_3vb_splitfaces, 
  0,
  reftet_1e_3vb_newelstypes, 
  reftet_1e_3vb_newels
};






/*
// HP_TET_1E_4V
int reftet_1e_4v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_1e_4v_splitfaces[][4] =
  {
    { 1, 2, 3, 17 },
    { 1, 2, 4, 18 },
    { 1, 3, 4, 19 },

    { 2, 1, 3, 20 },
    { 2, 1, 4, 21 },
    { 2, 3, 4, 22 },

    { 3, 1, 2, 23 },
    { 3, 1, 4, 24 },
    { 3, 2, 4, 25 },

    { 4, 1, 2, 26 },
    { 4, 1, 3, 27 },
    { 4, 2, 3, 28 },
    { 0, 0, 0, 0 },
  };
int reftet_1e_4v_splitelements[][5] =
  {
    { 1, 2, 3, 4, 29 },
    { 2, 3, 4, 1, 30 },
    { 3, 4, 1, 2, 31 },
    { 4, 1, 2, 3, 32 },
    { 0 },
  };
HPREF_ELEMENT_TYPE reftet_1e_4v_newelstypes[] =
{
  HP_HEX_1E_1V,
  HP_HEX_1E_1V,
  HP_HEX_0E_1V,
  HP_HEX_0E_1V,
  HP_PRISM_SINGEDGE, HP_PRISM, 
  HP_PRISM, HP_PRISM, 
  HP_PRISM, HP_PRISM, 
  HP_PRISM, HP_PRISM, 
  HP_PRISM, HP_PRISM, 
  HP_PRISM, HP_PRISM, 
  HP_PRISM,
  HP_PRISM,
  HP_PRISM,
  HP_PRISM,
  HP_TET,
  HP_NONE,
};
int reftet_1e_4v_newels[][8] =
{
  { 1, 5, 17, 6, 7, 18, 29, 19 },
  //  { 2, 9, 20, 8, 10, 22, 30, 21 },
  { 2, 8, 21, 10, 9, 20, 30, 22 },
  { 3, 11, 23, 12, 13, 24, 31, 25 },
  { 4, 15, 26, 14, 16, 28, 32, 27 },
  { 5, 17, 18, 8, 20, 21 },
  { 18, 17, 29, 21, 20, 30 },
  { 6, 19, 17,  11, 24, 23 },
  { 17, 19, 29,  23, 24, 31 },
  { 7, 18, 19, 14, 26, 27 },
  { 19, 18, 29, 27, 26, 32 },
  { 9, 20, 22, 12, 23, 25 },
  { 22, 20, 30, 25, 23, 31 },
  { 10, 22, 21, 15, 28, 26 },
  { 21, 22, 30, 26, 28, 32 },
  { 13, 24, 25, 16, 27, 28 },
  { 25, 24, 31, 28, 27, 32 },
  { 17, 20, 23, 29, 30, 31 },
  { 18, 26, 21, 29, 32, 30 },
  { 19, 24, 27, 29, 31, 32 },
  { 22, 28, 25, 30, 32, 31 },

  { 29, 30, 31, 32 },
};
HPRef_Struct reftet_1e_4v =
{
  HP_TET,
  reftet_1e_4v_splitedges, 
  reftet_1e_4v_splitfaces, 
  reftet_1e_4v_splitelements, 
  reftet_1e_4v_newelstypes, 
  reftet_1e_4v_newels
};
*/




// HP_TET_1E_4V
int reftet_1e_4v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_1e_4v_splitfaces[][4] =
  {
    { 1, 3, 4, 17 },
    { 2, 3, 4, 18 },

    { 3, 1, 4, 19 },
    { 3, 2, 4, 20 },

    { 4, 1, 3, 21 },
    { 4, 2, 3, 22 },
    { 0, 0, 0, 0 },
  };

HPREF_ELEMENT_TYPE reftet_1e_4v_newelstypes[] =
{
  HP_TET_1E_1VA,
  HP_TET_1E_1VA,
  //  HP_TET_1E_1VA,
  //  HP_TET_1E_1VA,
  HP_PRISM_SINGEDGE,
  HP_PRISM,
  HP_HEX, 
  HP_HEX, 
  HP_PRISM,
  HP_PRISM,

  HP_PYRAMID,
  HP_TET_0E_1V,

  HP_PYRAMID,
  HP_TET_0E_1V,

  HP_NONE,
};

int reftet_1e_4v_newels[][8] =
{
  { 1, 5, 6, 7 },
  { 2, 8, 10, 9 },

  { 5, 6, 7, 8, 9, 10 },
  { 7, 6, 17, 10, 9, 18 },

  { 7, 10, 18, 17, 14, 15, 22, 21 },
  { 9, 6, 17, 18, 12, 11, 19, 20 },

  { 17, 19, 21, 18, 20, 22 },
  { 16, 22, 21, 13, 20, 19 },

  { 14, 15, 22, 21, 16 },
  { 4, 14, 16, 15 },
  { 12, 11, 19, 20, 13 },
  { 3, 11, 12, 13 },



  { 1, 5, 17, 6, 7, 18, 29, 19 },
  //  { 2, 9, 20, 8, 10, 22, 30, 21 },
  { 2, 8, 21, 10, 9, 20, 30, 22 },
  { 3, 11, 23, 12, 13, 24, 31, 25 },
  { 4, 15, 26, 14, 16, 28, 32, 27 },
  { 5, 17, 18, 8, 20, 21 },
  { 18, 17, 29, 21, 20, 30 },
  { 6, 19, 17,  11, 24, 23 },
  { 17, 19, 29,  23, 24, 31 },
  { 7, 18, 19, 14, 26, 27 },
  { 19, 18, 29, 27, 26, 32 },
  { 9, 20, 22, 12, 23, 25 },
  { 22, 20, 30, 25, 23, 31 },
  { 10, 22, 21, 15, 28, 26 },
  { 21, 22, 30, 26, 28, 32 },
  { 13, 24, 25, 16, 27, 28 },
  { 25, 24, 31, 28, 27, 32 },
  { 17, 20, 23, 29, 30, 31 },
  { 18, 26, 21, 29, 32, 30 },
  { 19, 24, 27, 29, 31, 32 },
  { 22, 28, 25, 30, 32, 31 },

  { 29, 30, 31, 32 },
};
HPRef_Struct reftet_1e_4v =
{
  HP_TET,
  reftet_1e_4v_splitedges, 
  reftet_1e_4v_splitfaces, 
  0, 
  reftet_1e_4v_newelstypes, 
  reftet_1e_4v_newels
};













//  HP_TET_2EA_0V,  // 2 edges connected
int reftet_2ea_0v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 0, 0, 0 }
};
int reftet_2ea_0v_splitfaces[][4] =
  {
    { 1, 2, 3, 17 },
    { 0, 0, 0, 0 }
  };
HPREF_ELEMENT_TYPE reftet_2ea_0v_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_PRISM,
    HP_TET,
    HP_NONE,
  };
int reftet_2ea_0v_newels[][8] =
{
  { 1, 5, 17, 6, 7 },
  { 5, 17, 7, 2, 9, 10 },
  { 6, 7, 17, 3, 13, 12 },
  { 17, 9, 12, 7, 10, 13 },
  { 7, 10, 13, 4 },
};
HPRef_Struct reftet_2ea_0v =
{
  HP_TET,
  reftet_2ea_0v_splitedges, 
  reftet_2ea_0v_splitfaces, 
  0,
  reftet_2ea_0v_newelstypes, 
  reftet_2ea_0v_newels
};






//  HP_TET_2EA_1VA,  // 2 edges connected
int reftet_2ea_1va_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 0, 0, 0 }
};
int reftet_2ea_1va_splitfaces[][4] =
  {
    { 1, 2, 3, 17 },
    { 0, 0, 0, 0 }
  };
HPREF_ELEMENT_TYPE reftet_2ea_1va_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    HP_PRISM_SINGEDGE,
    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM,
    HP_TET,
    HP_NONE,
  };
int reftet_2ea_1va_newels[][8] =
{
  { 1, 5, 17, 6, 7 },
  { 5, 17, 7, 8, 9, 10 },
  { 2, 8, 10, 9 },
  { 6, 7, 17, 3, 13, 12 },
  { 17, 9, 12, 7, 10, 13 },
  { 7, 10, 13, 4 },
};
HPRef_Struct reftet_2ea_1va =
{
  HP_TET,
  reftet_2ea_1va_splitedges, 
  reftet_2ea_1va_splitfaces, 
  0,
  reftet_2ea_1va_newelstypes, 
  reftet_2ea_1va_newels
};








//  HP_TET_2EA_1VB, 
int reftet_2ea_1vb_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 0, 0, 0 }
};
int reftet_2ea_1vb_splitfaces[][4] =
  {
    { 1, 2, 3, 17 },
    { 0, 0, 0, 0 }
  };
HPREF_ELEMENT_TYPE reftet_2ea_1vb_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_PRISM,
    HP_TET,
    HP_NONE,
  };
int reftet_2ea_1vb_newels[][8] =
{
  { 1, 5, 17, 6, 7 },
  { 3, 11, 12, 13 },
  { 5, 17, 7, 2, 9, 10 },
  { 6, 7, 17, 11, 13, 12 },
  { 17, 9, 12, 7, 10, 13 },
  { 7, 10, 13, 4 },
};
HPRef_Struct reftet_2ea_1vb =
{
  HP_TET,
  reftet_2ea_1vb_splitedges, 
  reftet_2ea_1vb_splitfaces, 
  0,
  reftet_2ea_1vb_newelstypes, 
  reftet_2ea_1vb_newels
};






//  HP_TET_2EA_1VC,  // 2 edges connected
int reftet_2ea_1vc_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  //  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_2ea_1vc_splitfaces[][4] =
  {
    { 1, 2, 3, 17 },
    { 2, 3, 4, 18 },
    { 3, 4, 2, 19 },
    { 4, 2, 3, 20 },
    { 0, 0, 0, 0 }
  };
int reftet_2ea_1vc_splitelements[][5] =
  {
    { 1, 2, 3, 4, 21 },
    { 0, 0, 0, 0 }
  };
HPREF_ELEMENT_TYPE reftet_2ea_1vc_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    //    HP_TET_1E_1VA,
    HP_TET_0E_1V,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,

    HP_TET, HP_TET, HP_TET, HP_TET, 
    HP_PYRAMID, HP_PYRAMID, HP_PYRAMID, 
    HP_PYRAMID, HP_PYRAMID, HP_TET,
    HP_PYRAMID, HP_PYRAMID, HP_TET,
    //     HP_PRISM,
    //    HP_PRISM,
    HP_NONE,
  };
int reftet_2ea_1vc_newels[][8] =
{
  { 1, 5, 17, 6, 7 },
  // { 3, 11, 12, 13 },
  { 4, 15, 14, 16 }, 
  { 5, 17, 7, 2, 9, 10 },
  { 6, 7, 17, 3, 13, 12 },
 
  { 9, 10, 18, 21 },
  { 13, 12, 19, 21 },
  { 15, 16, 20, 21 },
  { 18, 20, 19, 21 },
  { 10, 15, 20, 18, 21 },
  { 13, 19, 20, 16, 21 },
  { 9, 18, 19, 12, 21 },
  
  { 7, 13, 16, 14, 21 },
  { 7, 14, 15, 10, 21 },
  { 9, 12, 17, 21 },
  { 7, 10, 9, 17, 21 },
  { 7, 17, 12, 13, 21 },
  { 14, 16, 15, 21 },
  //  { 17, 9, 12, 7, 10, 13 },
  //  { 7, 10, 13, 14, 15, 16 },
};
HPRef_Struct reftet_2ea_1vc =
{
  HP_TET,
  reftet_2ea_1vc_splitedges, 
  reftet_2ea_1vc_splitfaces, 
  reftet_2ea_1vc_splitelements, 
  reftet_2ea_1vc_newelstypes, 
  reftet_2ea_1vc_newels
};











 
//  HP_TET_2EA_2VA, 
int reftet_2ea_2va_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 0, 0, 0 }
};
int reftet_2ea_2va_splitfaces[][4] =
  {
    { 1, 2, 3, 17 },
    { 0, 0, 0, 0 }
  };
HPREF_ELEMENT_TYPE reftet_2ea_2va_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    HP_TET_1E_1VA,
    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_PRISM,
    HP_TET,
    HP_NONE,
  };
int reftet_2ea_2va_newels[][8] =
{
  { 1, 5, 17, 6, 7 },
  { 3, 11, 12, 13 },
  { 2, 8, 10, 9 },
  { 5, 17, 7, 8, 9, 10 },
  { 6, 7, 17, 11, 13, 12 },
  { 17, 9, 12, 7, 10, 13 },
  { 7, 10, 13, 4 },
};
HPRef_Struct reftet_2ea_2va =
{
  HP_TET,
  reftet_2ea_2va_splitedges, 
  reftet_2ea_2va_splitfaces, 
  0,
  reftet_2ea_2va_newelstypes, 
  reftet_2ea_2va_newels
};











//  HP_TET_2EA_2VB,  // 2 edges connected
int reftet_2ea_2vb_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  //  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_2ea_2vb_splitfaces[][4] =
  {
    { 1, 2, 3, 17 },
    { 2, 3, 4, 18 },
    { 3, 4, 2, 19 },
    { 4, 2, 3, 20 },
    { 0, 0, 0, 0 }
  };
int reftet_2ea_2vb_splitelements[][5] =
  {
    { 1, 2, 3, 4, 21 },
    { 0, 0, 0, 0 }
  };
HPREF_ELEMENT_TYPE reftet_2ea_2vb_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    HP_TET_1E_1VA,
    //  HP_TET_1E_1VA,
    HP_TET_0E_1V,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,

    HP_TET, HP_TET, HP_TET, HP_TET, 
    HP_PYRAMID, HP_PYRAMID, HP_PYRAMID, 
    HP_PYRAMID, HP_PYRAMID, HP_TET,
    HP_PYRAMID, HP_PYRAMID, HP_TET,
    //     HP_PRISM,
    //    HP_PRISM,
    HP_NONE,
  };
int reftet_2ea_2vb_newels[][8] =
{
  { 1, 5, 17, 6, 7 },
  { 2, 8, 10, 9 },
  //  { 3, 11, 12, 13 },
  { 4, 15, 14, 16 }, 
  { 5, 17, 7, 8, 9, 10 },
  { 6, 7, 17, 3, 13, 12 },
 
  { 9, 10, 18, 21 },
  { 13, 12, 19, 21 },
  { 15, 16, 20, 21 },
  { 18, 20, 19, 21 },
  { 10, 15, 20, 18, 21 },
  { 13, 19, 20, 16, 21 },
  { 9, 18, 19, 12, 21 },
  
  { 7, 13, 16, 14, 21 },
  { 7, 14, 15, 10, 21 },
  { 9, 12, 17, 21 },
  { 7, 10, 9, 17, 21 },
  { 7, 17, 12, 13, 21 },
  { 14, 16, 15, 21 },
  //  { 17, 9, 12, 7, 10, 13 },
  //  { 7, 10, 13, 14, 15, 16 },
};
HPRef_Struct reftet_2ea_2vb =
{
  HP_TET,
  reftet_2ea_2vb_splitedges, 
  reftet_2ea_2vb_splitfaces, 
  reftet_2ea_2vb_splitelements, 
  reftet_2ea_2vb_newelstypes, 
  reftet_2ea_2vb_newels
};







 


//  HP_TET_2EA_2VC,  // 2 edges connected
int reftet_2ea_2vc_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  //  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_2ea_2vc_splitfaces[][4] =
  {
    { 1, 2, 3, 17 },
    { 2, 3, 4, 18 },
    { 3, 4, 2, 19 },
    { 4, 2, 3, 20 },
    { 0, 0, 0, 0 }
  };
int reftet_2ea_2vc_splitelements[][5] =
  {
    { 1, 2, 3, 4, 21 },
    { 0, 0, 0, 0 }
  };
HPREF_ELEMENT_TYPE reftet_2ea_2vc_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    HP_TET_1E_1VA,
    //    HP_TET_1E_1VA,
    HP_TET_0E_1V,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,

    HP_TET, HP_TET, HP_TET, HP_TET, 
    HP_PYRAMID, HP_PYRAMID, HP_PYRAMID, 
    HP_PYRAMID, HP_PYRAMID, HP_TET,
    HP_PYRAMID, HP_PYRAMID, HP_TET,
    //     HP_PRISM,
    //    HP_PRISM,
    HP_NONE,
  };
int reftet_2ea_2vc_newels[][8] =
{
  { 1, 5, 17, 6, 7 },
  //  { 2, 8, 10, 9 },
  { 3, 11, 12, 13 },
  { 4, 15, 14, 16 }, 
  { 5, 17, 7, 2, 9, 10 },
  { 6, 7, 17, 11, 13, 12 },
 
  { 9, 10, 18, 21 },
  { 13, 12, 19, 21 },
  { 15, 16, 20, 21 },
  { 18, 20, 19, 21 },
  { 10, 15, 20, 18, 21 },
  { 13, 19, 20, 16, 21 },
  { 9, 18, 19, 12, 21 },
  
  { 7, 13, 16, 14, 21 },
  { 7, 14, 15, 10, 21 },
  { 9, 12, 17, 21 },
  { 7, 10, 9, 17, 21 },
  { 7, 17, 12, 13, 21 },
  { 14, 16, 15, 21 },
  //  { 17, 9, 12, 7, 10, 13 },
  //  { 7, 10, 13, 14, 15, 16 },
};
HPRef_Struct reftet_2ea_2vc =
{
  HP_TET,
  reftet_2ea_2vc_splitedges, 
  reftet_2ea_2vc_splitfaces, 
  reftet_2ea_2vc_splitelements, 
  reftet_2ea_2vc_newelstypes, 
  reftet_2ea_2vc_newels
};








//  HP_TET_2EA_3V,  // 2 edges connected
int reftet_2ea_3v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_2ea_3v_splitfaces[][4] =
  {
    { 1, 2, 3, 17 },
    { 2, 3, 4, 18 },
    { 3, 4, 2, 19 },
    { 4, 2, 3, 20 },
    { 0, 0, 0, 0 }
  };
int reftet_2ea_3v_splitelements[][5] =
  {
    { 1, 2, 3, 4, 21 },
    { 0, 0, 0, 0 }
  };
HPREF_ELEMENT_TYPE reftet_2ea_3v_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    HP_TET_1E_1VA,
    HP_TET_1E_1VA,
    HP_TET_0E_1V,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,

    HP_TET, HP_TET, HP_TET, HP_TET, 
    HP_PYRAMID, HP_PYRAMID, HP_PYRAMID, 
    HP_PYRAMID, HP_PYRAMID, HP_TET,
    HP_PYRAMID, HP_PYRAMID, HP_TET,
    //     HP_PRISM,
    //    HP_PRISM,
    HP_NONE,
  };
int reftet_2ea_3v_newels[][8] =
{
  { 1, 5, 17, 6, 7 },
  { 2, 8, 10, 9 },
  { 3, 11, 12, 13 },
  { 4, 15, 14, 16 }, 
  { 5, 17, 7, 8, 9, 10 },
  { 6, 7, 17, 11, 13, 12 },
 
  { 9, 10, 18, 21 },
  { 13, 12, 19, 21 },
  { 15, 16, 20, 21 },
  { 18, 20, 19, 21 },
  { 10, 15, 20, 18, 21 },
  { 13, 19, 20, 16, 21 },
  { 9, 18, 19, 12, 21 },
  
  { 7, 13, 16, 14, 21 },
  { 7, 14, 15, 10, 21 },
  { 9, 12, 17, 21 },
  { 7, 10, 9, 17, 21 },
  { 7, 17, 12, 13, 21 },
  { 14, 16, 15, 21 },
  //  { 17, 9, 12, 7, 10, 13 },
  //  { 7, 10, 13, 14, 15, 16 },
};
HPRef_Struct reftet_2ea_3v =
{
  HP_TET,
  reftet_2ea_3v_splitedges, 
  reftet_2ea_3v_splitfaces, 
  reftet_2ea_3v_splitelements, 
  reftet_2ea_3v_newelstypes, 
  reftet_2ea_3v_newels
};







//  HP_TET_2EB_0V,  // 2 opposite edges
int reftet_2eb_0v_splitedges[][3] =
{
  { 1, 3, 5 },
  { 1, 4, 6 },
  { 2, 3, 7 },
  { 2, 4, 8 },
  { 3, 1, 9 },
  { 3, 2, 10 },
  { 4, 1, 11 },
  { 4, 2, 12 },
  { 0, 0, 0 }
};

HPREF_ELEMENT_TYPE reftet_2eb_0v_newelstypes[] =
  {
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_HEX,
    HP_NONE,
  };
int reftet_2eb_0v_newels[][8] =
{
  { 1, 5, 6, 2, 7, 8 },
  { 3, 9, 10, 4, 11, 12 },
  { 6, 11, 12, 8, 5, 9, 10, 7 },
};
HPRef_Struct reftet_2eb_0v =
{
  HP_TET,
  reftet_2eb_0v_splitedges, 
  0, 0,
  reftet_2eb_0v_newelstypes, 
  reftet_2eb_0v_newels
};


//  HP_TET_2EB_1V,    // V1


//  HP_TET_2EB_1V,  // 2 opposite edges, V1
int reftet_2eb_1v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};

HPREF_ELEMENT_TYPE reftet_2eb_1v_newelstypes[] =
  {
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_TET_1E_1VA,
    //    HP_TET_1E_1VA,
    //    HP_TET_1E_1VA,
    //    HP_TET_1E_1VA,
    HP_HEX,
    HP_NONE,
  };
int reftet_2eb_1v_newels[][8] =
{
  { 5, 6, 7, 2, 9, 10 },
  { 4, 15, 14, 3, 12, 11 },
  { 1, 5, 6, 7 },
  //  { 2, 8, 10, 9 },
  //  { 3, 13, 11, 12 },
  //  { 4, 16, 15, 14 },
  { 7, 14, 15, 10, 6, 11, 12, 9 }
};
HPRef_Struct reftet_2eb_1v =
{
  HP_TET,
  reftet_2eb_1v_splitedges, 
  0, 0,
  reftet_2eb_1v_newelstypes, 
  reftet_2eb_1v_newels
};



//  HP_TET_2EB_2VA,  // 2 opposite edges, V1,2
int reftet_2eb_2va_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};

HPREF_ELEMENT_TYPE reftet_2eb_2va_newelstypes[] =
  {
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_TET_1E_1VA,
    HP_TET_1E_1VA,
    //    HP_TET_1E_1VA,
    //    HP_TET_1E_1VA,
    HP_HEX,
    HP_NONE,
  };
int reftet_2eb_2va_newels[][8] =
{
  { 5, 6, 7, 8, 9, 10 },
  { 4, 15, 14, 3, 12, 11 },
  { 1, 5, 6, 7 },
  { 2, 8, 10, 9 },
  //  { 3, 13, 11, 12 },
  //  { 4, 16, 15, 14 },
  { 7, 14, 15, 10, 6, 11, 12, 9 }
};
HPRef_Struct reftet_2eb_2va =
{
  HP_TET,
  reftet_2eb_2va_splitedges, 
  0, 0,
  reftet_2eb_2va_newelstypes, 
  reftet_2eb_2va_newels
};


//  HP_TET_2EB_2VB,   // V1,3
int reftet_2eb_2vb_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};

HPREF_ELEMENT_TYPE reftet_2eb_2vb_newelstypes[] =
  {
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_TET_1E_1VA,
    // HP_TET_1E_1VA,
    HP_TET_1E_1VA,
    // HP_TET_1E_1VA,
    HP_HEX,
    HP_NONE,
  };
int reftet_2eb_2vb_newels[][8] =
{
  { 5, 6, 7, 2, 9, 10 },
  { 4, 15, 14, 13, 12, 11 },
  { 1, 5, 6, 7 },
  // { 2, 8, 10, 9 },
  { 3, 13, 11, 12 },
  // { 4, 16, 15, 14 },
  { 7, 14, 15, 10, 6, 11, 12, 9 }
};
HPRef_Struct reftet_2eb_2vb =
{
  HP_TET,
  reftet_2eb_2vb_splitedges, 
  0, 0,
  reftet_2eb_2vb_newelstypes, 
  reftet_2eb_2vb_newels
};




//  HP_TET_2EB_2VC,   // V1,4
int reftet_2eb_2vc_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};

HPREF_ELEMENT_TYPE reftet_2eb_2vc_newelstypes[] =
  {
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_TET_1E_1VA,
    // HP_TET_1E_1VA,
    // HP_TET_1E_1VA,
    HP_TET_1E_1VA,
    HP_HEX,
    HP_NONE,
  };
int reftet_2eb_2vc_newels[][8] =
{
  { 5, 6, 7, 2, 9, 10 },
  { 16, 15, 14, 3, 12, 11 },
  { 1, 5, 6, 7 },
  // { 2, 8, 10, 9 },
  // { 3, 13, 11, 12 },
  { 4, 16, 15, 14 },
  { 7, 14, 15, 10, 6, 11, 12, 9 }
};
HPRef_Struct reftet_2eb_2vc =
{
  HP_TET,
  reftet_2eb_2vc_splitedges, 
  0, 0,
  reftet_2eb_2vc_newelstypes, 
  reftet_2eb_2vc_newels
};






//  HP_TET_2EB_3V,    // V1,2,3
int reftet_2eb_3v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};

HPREF_ELEMENT_TYPE reftet_2eb_3v_newelstypes[] =
  {
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_TET_1E_1VA,
    HP_TET_1E_1VA,
    HP_TET_1E_1VA,
    // HP_TET_1E_1VA,
    HP_HEX,
    HP_NONE,
  };
int reftet_2eb_3v_newels[][8] =
{
  { 5, 6, 7, 8, 9, 10 },
  { 4, 15, 14, 13, 12, 11 },
  { 1, 5, 6, 7 },
  { 2, 8, 10, 9 },
  { 3, 13, 11, 12 },
  // { 4, 16, 15, 14 },
  { 7, 14, 15, 10, 6, 11, 12, 9 }
};
HPRef_Struct reftet_2eb_3v =
{
  HP_TET,
  reftet_2eb_3v_splitedges, 
  0, 0,
  reftet_2eb_3v_newelstypes, 
  reftet_2eb_3v_newels
};






//  HP_TET_2EB_4V,  // 2 opposite edges
int reftet_2eb_4v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};

HPREF_ELEMENT_TYPE reftet_2eb_4v_newelstypes[] =
  {
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_TET_1E_1VA,
    HP_TET_1E_1VA,
    HP_TET_1E_1VA,
    HP_TET_1E_1VA,
    HP_HEX,
    HP_NONE,
  };
int reftet_2eb_4v_newels[][8] =
{
  { 5, 6, 7, 8, 9, 10 },
  { 16, 15, 14, 13, 12, 11 },
  { 1, 5, 6, 7 },
  { 2, 8, 10, 9 },
  { 3, 13, 11, 12 },
  { 4, 16, 15, 14 },
  { 7, 14, 15, 10, 6, 11, 12, 9 }
};
HPRef_Struct reftet_2eb_4v =
{
  HP_TET,
  reftet_2eb_4v_splitedges, 
  0, 0,
  reftet_2eb_4v_newelstypes, 
  reftet_2eb_4v_newels
};

















//  HP_TET_3EA_0V,  
int reftet_3ea_0v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 3, 8 },
  { 2, 4, 9 },
  { 3, 2, 10 },
  { 3, 4, 11 },
  { 4, 2, 12 },
  { 4, 3, 13 },
  { 0, 0, 0 }
};
int reftet_3ea_0v_splitfaces[][4] =
  {
    { 1, 2, 3, 14 },
    { 1, 2, 4, 15 },
    { 1, 3, 4, 16 },
    { 2, 3, 4, 17 },
    { 3, 4, 2, 18 },
    { 4, 2, 3, 19 },
    { 0, 0, 0, 0 }
  };
int reftet_3ea_0v_splitelements[][5] =
  {
    { 1, 2, 3, 4, 20 },
    { 0 },
  };

HPREF_ELEMENT_TYPE reftet_3ea_0v_newelstypes[] =
  {
    HP_HEX_3E_0V,
    HP_HEX_1E_0V,
    HP_HEX_1E_0V,
    HP_HEX_1E_0V,
    HP_PRISM,
    HP_PRISM,
    HP_PRISM,
    HP_TET,
    HP_NONE,
  };
int reftet_3ea_0v_newels[][8] =
{
  { 1, 5, 14, 6, 7, 15, 20, 16 },
  { 5, 2, 8, 14, 15, 9, 17, 20 },
  { 3, 6, 14, 10, 11, 16, 20, 18 },
  { 7, 4, 12, 15, 16, 13, 19, 20 },
  { 11, 13, 16, 18, 19, 20 },
  { 15, 12, 9, 20, 19, 17 },
  { 8, 10, 14, 17, 18, 20 },
  { 20, 17, 18, 19 },
};
HPRef_Struct reftet_3ea_0v =
{
  HP_TET,
  reftet_3ea_0v_splitedges, 
  reftet_3ea_0v_splitfaces, 
  reftet_3ea_0v_splitelements, 
  reftet_3ea_0v_newelstypes, 
  reftet_3ea_0v_newels
};










//  HP_TET_3EA_1V,  
int reftet_3ea_1v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 3, 8 },
  { 2, 4, 9 },
  { 3, 2, 10 },
  { 3, 4, 11 },
  { 4, 2, 12 },
  { 4, 3, 13 },
  { 2, 1, 21 },
  { 3, 1, 22 },
  { 4, 1, 23 },
  { 0, 0, 0 }
};
int reftet_3ea_1v_splitfaces[][4] =
  {
    { 1, 2, 3, 14 },
    { 1, 2, 4, 15 },
    { 1, 3, 4, 16 },
    { 2, 3, 4, 17 },
    { 3, 4, 2, 18 },
    { 4, 2, 3, 19 },
    { 0, 0, 0, 0 }
  };
int reftet_3ea_1v_splitelements[][5] =
  {
    { 1, 2, 3, 4, 20 },
    { 0 },
  };

HPREF_ELEMENT_TYPE reftet_3ea_1v_newelstypes[] =
  {
    HP_HEX_3E_0V,
    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM,
    //    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM,
    //    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM,

    HP_PRISM,
    HP_PRISM,
    HP_PRISM,
    HP_TET,
    HP_NONE,
  };
int reftet_3ea_1v_newels[][8] =
{
  { 1, 5, 14, 6, 7, 15, 20, 16 },

  { 2, 21, 9, 8 },
  { 5, 14, 15, 21, 8, 9 },
  { 15, 14, 20, 9, 8, 17 },
  //  { 3, 22, 10, 11 },
  //  { 6, 16, 14, 22, 11, 10 },
  { 6, 16, 14, 3, 11, 10 },
  { 14, 16, 20, 10, 11, 18 },
  //  { 4, 23, 13, 12 },
  //  { 7, 15, 16, 23, 12, 13 },
  { 7, 15, 16, 4, 12, 13 },
  { 16, 15, 20, 13, 12, 19 },

  { 11, 13, 16, 18, 19, 20 },
  { 15, 12, 9, 20, 19, 17 },
  { 8, 10, 14, 17, 18, 20 },
  { 20, 17, 18, 19 },
};
HPRef_Struct reftet_3ea_1v =
{
  HP_TET,
  reftet_3ea_1v_splitedges, 
  reftet_3ea_1v_splitfaces, 
  reftet_3ea_1v_splitelements, 
  reftet_3ea_1v_newelstypes, 
  reftet_3ea_1v_newels
};










//  HP_TET_3EA_2V,  
int reftet_3ea_2v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 3, 8 },
  { 2, 4, 9 },
  { 3, 2, 10 },
  { 3, 4, 11 },
  { 4, 2, 12 },
  { 4, 3, 13 },
  { 2, 1, 21 },
  { 3, 1, 22 },
  { 4, 1, 23 },
  { 0, 0, 0 }
};
int reftet_3ea_2v_splitfaces[][4] =
  {
    { 1, 2, 3, 14 },
    { 1, 2, 4, 15 },
    { 1, 3, 4, 16 },
    { 2, 3, 4, 17 },
    { 3, 4, 2, 18 },
    { 4, 2, 3, 19 },
    { 0, 0, 0, 0 }
  };
int reftet_3ea_2v_splitelements[][5] =
  {
    { 1, 2, 3, 4, 20 },
    { 0 },
  };

HPREF_ELEMENT_TYPE reftet_3ea_2v_newelstypes[] =
  {
    HP_HEX_3E_0V,
    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM,
    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM,
    //    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM,

    HP_PRISM,
    HP_PRISM,
    HP_PRISM,
    HP_TET,
    HP_NONE,
  };
int reftet_3ea_2v_newels[][8] =
{
  { 1, 5, 14, 6, 7, 15, 20, 16 },

  { 2, 21, 9, 8 },
  { 5, 14, 15, 21, 8, 9 },
  { 15, 14, 20, 9, 8, 17 },
  { 3, 22, 10, 11 },
  { 6, 16, 14, 22, 11, 10 },
  { 14, 16, 20, 10, 11, 18 },
  //  { 4, 23, 13, 12 },
  { 7, 15, 16, 4, 12, 13 },
  { 16, 15, 20, 13, 12, 19 },

  { 11, 13, 16, 18, 19, 20 },
  { 15, 12, 9, 20, 19, 17 },
  { 8, 10, 14, 17, 18, 20 },
  { 20, 17, 18, 19 },
};
HPRef_Struct reftet_3ea_2v =
{
  HP_TET,
  reftet_3ea_2v_splitedges, 
  reftet_3ea_2v_splitfaces, 
  reftet_3ea_2v_splitelements, 
  reftet_3ea_2v_newelstypes, 
  reftet_3ea_2v_newels
};








//  HP_TET_3EA_3V,  
int reftet_3ea_3v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 3, 8 },
  { 2, 4, 9 },
  { 3, 2, 10 },
  { 3, 4, 11 },
  { 4, 2, 12 },
  { 4, 3, 13 },
  { 2, 1, 21 },
  { 3, 1, 22 },
  { 4, 1, 23 },
  { 0, 0, 0 }
};
int reftet_3ea_3v_splitfaces[][4] =
  {
    { 1, 2, 3, 14 },
    { 1, 2, 4, 15 },
    { 1, 3, 4, 16 },
    { 2, 3, 4, 17 },
    { 3, 4, 2, 18 },
    { 4, 2, 3, 19 },
    { 0, 0, 0, 0 }
  };
int reftet_3ea_3v_splitelements[][5] =
  {
    { 1, 2, 3, 4, 20 },
    { 0 },
  };

HPREF_ELEMENT_TYPE reftet_3ea_3v_newelstypes[] =
  {
    HP_HEX_3E_0V,
    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM,
    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM,
    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM,

    HP_PRISM,
    HP_PRISM,
    HP_PRISM,
    HP_TET,
    HP_NONE,
  };
int reftet_3ea_3v_newels[][8] =
{
  { 1, 5, 14, 6, 7, 15, 20, 16 },

  { 2, 21, 9, 8 },
  { 5, 14, 15, 21, 8, 9 },
  { 15, 14, 20, 9, 8, 17 },
  { 3, 22, 10, 11 },
  { 6, 16, 14, 22, 11, 10 },
  { 14, 16, 20, 10, 11, 18 },
  { 4, 23, 13, 12 },
  { 7, 15, 16, 23, 12, 13 },
  { 16, 15, 20, 13, 12, 19 },

  { 11, 13, 16, 18, 19, 20 },
  { 15, 12, 9, 20, 19, 17 },
  { 8, 10, 14, 17, 18, 20 },
  { 20, 17, 18, 19 },
};
HPRef_Struct reftet_3ea_3v =
{
  HP_TET,
  reftet_3ea_3v_splitedges, 
  reftet_3ea_3v_splitfaces, 
  reftet_3ea_3v_splitelements, 
  reftet_3ea_3v_newelstypes, 
  reftet_3ea_3v_newels
};







//  HP_TET_3EV_0V,  
int reftet_3eb_0v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  //  { 3, 2, 12 },
  { 3, 4, 13 },
  //  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_3eb_0v_splitfaces[][4] =
  {
    { 1, 2, 4, 17 },
    { 2, 1, 3, 18 },
    { 0, 0, 0, 0 }
  };
int reftet_3eb_0v_splitelements[][5] =
  {
    { 1, 2, 3, 4, 20 },
    { 0 },
  };

HPREF_ELEMENT_TYPE reftet_3eb_0v_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    HP_PYRAMID_EDGES,
    //    HP_TET_1E_1VA,
    //    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    
    HP_PYRAMID,
    HP_PYRAMID,
    HP_TET,
    HP_TET,
    HP_PYRAMID,
    HP_PYRAMID,
    HP_PYRAMID,
    HP_NONE,
  };
int reftet_3eb_0v_newels[][8] =
{
  { 1, 7, 17, 5, 6 },
  { 2, 9, 18, 8, 10 },
  //  { 3, 12, 13, 11 },
  //  { 4, 14, 16, 15 },
  { 5, 6, 17, 8, 18, 10 },
  { 7, 17, 6, 4, 15, 16 },
  { 9, 18, 10, 3, 11, 13 },
  
  { 10, 15, 16, 13, 20 },
  { 6, 11, 13, 16, 20 },
  { 10, 17, 15, 20 },
  { 6, 18, 11, 20 },
  { 6, 17, 10, 18, 20 },
  { 6, 16, 15, 17, 20 },
  { 18, 10, 13, 11, 20 },
};
HPRef_Struct reftet_3eb_0v =
{
  HP_TET,
  reftet_3eb_0v_splitedges, 
  reftet_3eb_0v_splitfaces, 
  reftet_3eb_0v_splitelements, 
  reftet_3eb_0v_newelstypes, 
  reftet_3eb_0v_newels
};









//  HP_TET_3EV_1V,  
int reftet_3eb_1v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  //  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_3eb_1v_splitfaces[][4] =
  {
    { 1, 2, 4, 17 },
    { 2, 1, 3, 18 },
    { 0, 0, 0, 0 }
  };
int reftet_3eb_1v_splitelements[][5] =
  {
    { 1, 2, 3, 4, 20 },
    { 0 },
  };

HPREF_ELEMENT_TYPE reftet_3eb_1v_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    HP_PYRAMID_EDGES,
    HP_TET_1E_1VA,
    //    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    
    HP_PYRAMID,
    HP_PYRAMID,
    HP_TET,
    HP_TET,
    HP_PYRAMID,
    HP_PYRAMID,
    HP_PYRAMID,
    HP_NONE,
  };
int reftet_3eb_1v_newels[][8] =
{
  { 1, 7, 17, 5, 6 },
  { 2, 9, 18, 8, 10 },
  { 3, 12, 13, 11 },
  //  { 4, 14, 16, 15 },
  { 5, 6, 17, 8, 18, 10 },
  { 7, 17, 6, 4, 15, 16 },
  { 9, 18, 10, 12, 11, 13 },
  
  { 10, 15, 16, 13, 20 },
  { 6, 11, 13, 16, 20 },
  { 10, 17, 15, 20 },
  { 6, 18, 11, 20 },
  { 6, 17, 10, 18, 20 },
  { 6, 16, 15, 17, 20 },
  { 18, 10, 13, 11, 20 },
};
HPRef_Struct reftet_3eb_1v =
{
  HP_TET,
  reftet_3eb_1v_splitedges, 
  reftet_3eb_1v_splitfaces, 
  reftet_3eb_1v_splitelements, 
  reftet_3eb_1v_newelstypes, 
  reftet_3eb_1v_newels
};








//  HP_TET_3EV_2V,  
int reftet_3eb_2v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_3eb_2v_splitfaces[][4] =
  {
    { 1, 2, 4, 17 },
    { 2, 1, 3, 18 },
    { 0, 0, 0, 0 }
  };
int reftet_3eb_2v_splitelements[][5] =
  {
    { 1, 2, 3, 4, 20 },
    { 0 },
  };

HPREF_ELEMENT_TYPE reftet_3eb_2v_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    HP_PYRAMID_EDGES,
    HP_TET_1E_1VA,
    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    
    HP_PYRAMID,
    HP_PYRAMID,
    HP_TET,
    HP_TET,
    HP_PYRAMID,
    HP_PYRAMID,
    HP_PYRAMID,
    HP_NONE,
  };
int reftet_3eb_2v_newels[][8] =
{
  { 1, 7, 17, 5, 6 },
  { 2, 9, 18, 8, 10 },
  { 3, 12, 13, 11 },
  { 4, 14, 16, 15 },
  { 5, 6, 17, 8, 18, 10 },
  { 7, 17, 6, 14, 15, 16 },
  { 9, 18, 10, 12, 11, 13 },
  
  { 10, 15, 16, 13, 20 },
  { 6, 11, 13, 16, 20 },
  { 10, 17, 15, 20 },
  { 6, 18, 11, 20 },
  { 6, 17, 10, 18, 20 },
  { 6, 16, 15, 17, 20 },
  { 18, 10, 13, 11, 20 },
};
HPRef_Struct reftet_3eb_2v =
{
  HP_TET,
  reftet_3eb_2v_splitedges, 
  reftet_3eb_2v_splitfaces, 
  reftet_3eb_2v_splitelements, 
  reftet_3eb_2v_newelstypes, 
  reftet_3eb_2v_newels
};













//  HP_TET_3EC_0V,  
int reftet_3ec_0v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  //  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  //  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_3ec_0v_splitfaces[][4] =
  {
    { 1, 2, 3, 17 },
    { 2, 1, 4, 18 },
    { 0, 0, 0, 0 }
  };
int reftet_3ec_0v_splitelements[][5] =
  {
    { 1, 2, 3, 4, 20 },
    { 0 },
  };

HPREF_ELEMENT_TYPE reftet_3ec_0v_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    HP_PYRAMID_EDGES,
    //    HP_TET_1E_1VA,
    //    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    
    HP_PYRAMID,
    HP_PYRAMID,
    HP_TET,
    HP_TET,
    HP_PYRAMID,
    HP_PYRAMID,
    HP_PYRAMID,
    HP_NONE,
  };
int reftet_3ec_0v_newels[][8] =
{
  { 1, 5, 17, 6, 7 },
  { 2, 8, 18, 10, 9 },
  //  { 3, 11, 12, 13 },
  //  { 4, 15, 14, 16 },
  { 5, 17, 7, 8, 9, 18 },
  { 6, 7, 17, 3, 13, 12 },
  { 10, 9, 18, 4, 16, 14 },
  
  { 9, 16, 13, 12, 20 },
  { 7, 13, 16, 14, 20 },
  { 7, 14, 18, 20 },
  { 9, 12, 17, 20 },
  { 17, 7, 18, 9, 20 },
  { 7, 17, 12, 13, 20 },
  { 9, 18, 14, 16, 20 },
};
HPRef_Struct reftet_3ec_0v =
{
  HP_TET,
  reftet_3ec_0v_splitedges, 
  reftet_3ec_0v_splitfaces, 
  reftet_3ec_0v_splitelements, 
  reftet_3ec_0v_newelstypes, 
  reftet_3ec_0v_newels
};






 


//  HP_TET_3EC_1V,  
int reftet_3ec_1v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  // { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_3ec_1v_splitfaces[][4] =
  {
    { 1, 2, 3, 17 },
    { 2, 1, 4, 18 },
    { 0, 0, 0, 0 }
  };
int reftet_3ec_1v_splitelements[][5] =
  {
    { 1, 2, 3, 4, 20 },
    { 0 },
  };

HPREF_ELEMENT_TYPE reftet_3ec_1v_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    HP_PYRAMID_EDGES,
    HP_TET_1E_1VA,
    //    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    
    HP_PYRAMID,
    HP_PYRAMID,
    HP_TET,
    HP_TET,
    HP_PYRAMID,
    HP_PYRAMID,
    HP_PYRAMID,
    HP_NONE,
  };
int reftet_3ec_1v_newels[][8] =
{
  { 1, 5, 17, 6, 7 },
  { 2, 8, 18, 10, 9 },
  { 3, 11, 12, 13 },
  //  { 4, 15, 14, 16 },
  { 5, 17, 7, 8, 9, 18 },
  { 6, 7, 17, 11, 13, 12 },
  { 10, 9, 18, 4, 16, 14 },
  
  { 9, 16, 13, 12, 20 },
  { 7, 13, 16, 14, 20 },
  { 7, 14, 18, 20 },
  { 9, 12, 17, 20 },
  { 17, 7, 18, 9, 20 },
  { 7, 17, 12, 13, 20 },
  { 9, 18, 14, 16, 20 },
};
HPRef_Struct reftet_3ec_1v =
{
  HP_TET,
  reftet_3ec_1v_splitedges, 
  reftet_3ec_1v_splitfaces, 
  reftet_3ec_1v_splitelements, 
  reftet_3ec_1v_newelstypes, 
  reftet_3ec_1v_newels
};








//  HP_TET_3EC_2V,  
int reftet_3ec_2v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 4, 1, 14 },
  { 4, 2, 15 },
  { 4, 3, 16 },
  { 0, 0, 0 }
};
int reftet_3ec_2v_splitfaces[][4] =
  {
    { 1, 2, 3, 17 },
    { 2, 1, 4, 18 },
    { 0, 0, 0, 0 }
  };
int reftet_3ec_2v_splitelements[][5] =
  {
    { 1, 2, 3, 4, 20 },
    { 0 },
  };

HPREF_ELEMENT_TYPE reftet_3ec_2v_newelstypes[] =
  {
    HP_PYRAMID_EDGES,
    HP_PYRAMID_EDGES,
    HP_TET_1E_1VA,
    HP_TET_1E_1VA,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    
    HP_PYRAMID,
    HP_PYRAMID,
    HP_TET,
    HP_TET,
    HP_PYRAMID,
    HP_PYRAMID,
    HP_PYRAMID,
    HP_NONE,
  };
int reftet_3ec_2v_newels[][8] =
{
  { 1, 5, 17, 6, 7 },
  { 2, 8, 18, 10, 9 },
  { 3, 11, 12, 13 },
  { 4, 15, 14, 16 },
  { 5, 17, 7, 8, 9, 18 },
  { 6, 7, 17, 11, 13, 12 },
  { 10, 9, 18, 15, 16, 14 },
  
  { 9, 16, 13, 12, 20 },
  { 7, 13, 16, 14, 20 },
  { 7, 14, 18, 20 },
  { 9, 12, 17, 20 },
  { 17, 7, 18, 9, 20 },
  { 7, 17, 12, 13, 20 },
  { 9, 18, 14, 16, 20 },
};
HPRef_Struct reftet_3ec_2v =
{
  HP_TET,
  reftet_3ec_2v_splitedges, 
  reftet_3ec_2v_splitfaces, 
  reftet_3ec_2v_splitelements, 
  reftet_3ec_2v_newelstypes, 
  reftet_3ec_2v_newels
};





//  HP_TET_3ED_3V,  
int reftet_3ed_3v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 2, 3, 9 },
  { 2, 4, 10 },
  { 3, 1, 11 },
  { 3, 2, 12 },
  { 3, 4, 13 },
  { 0, 0, 0 }
};
int reftet_3ed_3v_splitfaces[][4] =
  {
    { 1, 2, 3, 14 },
    { 2, 3, 1, 15 },
    { 3, 1, 2, 16 },
    { 0, 0, 0, 0 }
  };
int reftet_3ed_3v_splitelements[][5] =
  {
    { 0 },
  };

HPREF_ELEMENT_TYPE reftet_3ed_3v_newelstypes[] =
  {
    HP_TET,
    HP_PRISM,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_PRISM_SINGEDGE,
    HP_TET_1E_1VA,
    HP_TET_1E_1VA,
    HP_TET_1E_1VA,    
    HP_TET_1E_1VA,
    HP_TET_1E_1VA,    
    HP_TET_1E_1VA,    
    HP_NONE,
  };
int reftet_3ed_3v_newels[][8] =
{
  { 7, 10, 13, 4 },
  { 14, 15, 16, 7, 10, 13 },
  { 5, 14, 7, 8, 15, 10 },
  { 9, 15, 10, 12, 16, 13 },  
  { 11, 16, 13, 6, 14, 7 },
  { 1, 5, 14, 7 },
  { 1, 6, 7, 14 },
  { 2, 8, 10, 15 },
  { 2, 9, 15, 10 },
  { 3, 12, 13, 16 },
  { 3, 11, 16, 13 }
};

HPRef_Struct reftet_3ed_3v =
{
  HP_TET,
  reftet_3ed_3v_splitedges, 
  reftet_3ed_3v_splitfaces, 
  reftet_3ed_3v_splitelements, 
  reftet_3ed_3v_newelstypes, 
  reftet_3ed_3v_newels
};











/* ************************ 1 singular face ******************** */


// HP_TET_1F_0E_0V
int reftet_1f_0e_0v_splitedges[][3] =
{
  { 2, 1, 5 },
  { 3, 1, 6 },
  { 4, 1, 7 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftet_1f_0e_0v_newelstypes[] =
{
  HP_PRISM_1FA_0E_0V,
  HP_TET,
  HP_NONE,
};
int reftet_1f_0e_0v_newels[][8] =
{
  { 3, 2, 4, 6, 5, 7 },
  { 5, 7, 6, 1 }
};
HPRef_Struct reftet_1f_0e_0v =
{
  HP_TET,
  reftet_1f_0e_0v_splitedges, 
  0, 0,
  reftet_1f_0e_0v_newelstypes, 
  reftet_1f_0e_0v_newels
};




/*
// HP_TET_1F_0E_1VA    ... singular vertex in face
int reftet_1f_0e_1va_splitedges[][3] =
{
  { 2, 1, 5 },
  { 2, 3, 6 },
  { 2, 4, 7 },
  { 3, 1, 8 },
  { 4, 1, 9 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftet_1f_0e_1va_newelstypes[] =
{
  // HP_HEX_1F_0E_0V,
  HP_HEX7_1FA,
  HP_TET_1F_0E_1VA,
  HP_TET,
  HP_NONE,
};
int reftet_1f_0e_1va_newels[][8] =
{
  // { 3, 6, 7, 4, 8, 5, 5, 9 },
  { 4, 3, 6, 7, 9, 8, 5 },  
  { 5, 2, 6, 7 },
  { 5, 9, 8, 1 },
};
HPRef_Struct reftet_1f_0e_1va =
{
  HP_TET,
  reftet_1f_0e_1va_splitedges, 
  0, 0,
  reftet_1f_0e_1va_newelstypes, 
  reftet_1f_0e_1va_newels
};
*/


HPRefStruct<HP_TET> reftet_1f_0e_1va
  {
    HP_TET_1F_0E_1VA,
    {
      El(HP_HEX7_1FA, { V4, V3, E23, E24, E41, E31, E21 }),
      El(HP_TET_1F_0E_1VA, { E21, V2, E23, E24 }),
      El(HP_TET, {  E21, E41, E31, V1 })
    }
  };



// HP_TET_1F_0E_1VB    ... singular vertex not in face
int reftet_1f_0e_1vb_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 1, 4, 7 },
  { 2, 1, 8 },
  { 3, 1, 9 },
  { 4, 1, 10 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftet_1f_0e_1vb_newelstypes[] =
{
  HP_PRISM_1FA_0E_0V,
  HP_PRISM,
  HP_TET_0E_1V,
  HP_NONE,
};
int reftet_1f_0e_1vb_newels[][8] =
{
  { 2, 4, 3, 8, 10, 9 },
  { 8, 10, 9, 5, 7, 6 }, 
  { 1, 5, 6, 7 },
};
HPRef_Struct reftet_1f_0e_1vb =
{
  HP_TET,
  reftet_1f_0e_1vb_splitedges, 
  0, 0,
  reftet_1f_0e_1vb_newelstypes, 
  reftet_1f_0e_1vb_newels
};







// HP_TET_1F_0E_2V    ... face 234, sing verts v2,v3
int reftet_1f_0e_2v_splitedges[][3] =
{
  { 2, 1, 5 },
  { 2, 3, 6 },
  { 2, 4, 7 },
  { 3, 1, 8 },
  { 3, 2, 9 },
  { 3, 4, 10 },
  { 4, 1, 11 },
  { 0, 0, 0 }
};
HPREF_ELEMENT_TYPE reftet_1f_0e_2v_newelstypes[] =
{
  HP_TET_0E_1V,
  HP_PRISM_1FA_0E_0V,
  HP_PRISM_1FB_0E_0V,
  HP_TET_1F_0E_1VA,
  HP_TET_1F_0E_1VA,
  HP_NONE,
};
int reftet_1f_0e_2v_newels[][8] =
{
  { 1, 5, 8, 11 },
  { 4, 10, 7, 11, 8, 5 },
  { 9, 10, 8, 6, 7, 5 },
  { 5, 2, 6, 7 },
  { 8, 3, 10, 9 }
};
HPRef_Struct reftet_1f_0e_2v =
{
  HP_TET,
  reftet_1f_0e_2v_splitedges, 
  0, 0,
  reftet_1f_0e_2v_newelstypes, 
  reftet_1f_0e_2v_newels
};








HPRefStruct<HP_TET> reftet_1f_0e_3v
  {
    HP_TET_1F_0E_3V,
    {
      El(HP_TET,   { V1, E21, E31, E41 }),
      El(HP_PRISM_1FA_0E_0V, { F234, F423, F324, E21, E41, E31 }),
      El(HP_PRISM_1FB_1EA_0V, { E32, F324, E31, E23, F234, E21  }),
      El(HP_PRISM_1FB_1EA_0V, { E43, F423, E41, E34, F324, E31  }),
      El(HP_PRISM_1FB_1EA_0V, { E24, F234, E21, E42, F423, E41  }),
      El(HP_TET_1F_0E_0V, { E21, E24, E23, F234 }),
      El(HP_TET_1F_0E_1VA, { E21, V2, E23, E24 }),            
      El(HP_TET_1F_0E_0V, { E31, E32, E34, F324 }),
      El(HP_TET_1F_0E_1VA, { E31, V3, E34, E32 }),            
      El(HP_TET_1F_0E_0V, { E41, E43, E42, F423 }),
      El(HP_TET_1F_0E_1VA, { E41, V4, E42, E43 }),            
    }
  };











// HP_TET_1F_1EA_0V  ... sing edge is 1..2
int reftet_1f_1ea_0v_splitedges[][3] =
{
  { 1, 3, 5 },
  { 1, 4, 6 },
  { 2, 1, 7 },
  { 2, 3, 8 },
  { 2, 4, 9 },
  { 3, 1, 10 },
  { 4, 1, 11 },
  { 0, 0, 0 }
};

int reftet_1f_1ea_0v_splitfaces[][4] =
  {
    { 2, 1, 3, 12 },
    { 2, 1, 4, 13 },
    { 0, 0, 0, 0 }
  };


HPREF_ELEMENT_TYPE reftet_1f_1ea_0v_newelstypes[] =
{
  HP_HEX_1F_0E_0V,
  //  HP_PRISM,
  HP_PYRAMID_1FB_0E_1VA,
  HP_TET_1E_1VA,
  HP_PRISM_SINGEDGE,
  HP_PRISM,
  HP_NONE,
};
int reftet_1f_1ea_0v_newels[][8] =
{
  { 3, 8, 9, 4, 10, 12, 13, 11 },
  // { 2, 9, 8, 7, 13, 12 },
  { 8, 9, 13, 12, 2 },
  { 2, 7, 13, 12 },
  { 7, 13, 12, 1, 6, 5 },
  { 6, 11, 13, 5, 10, 12 }
};
HPRef_Struct reftet_1f_1ea_0v =
{
  HP_TET,
  reftet_1f_1ea_0v_splitedges, 
  reftet_1f_1ea_0v_splitfaces, 
  0, 
  reftet_1f_1ea_0v_newelstypes, 
  reftet_1f_1ea_0v_newels
};








// HP_TET_1F_1EB_0V     singular edge in face, edge is 2-3
int reftet_1f_1eb_0v_splitedges[][3] =
{
  { 2, 1, 5 },
  { 2, 4, 6 },
  { 3, 1, 7 },
  { 3, 4, 8 },
  { 4, 1, 9 },
  { 0, 0, 0 }
};


HPREF_ELEMENT_TYPE reftet_1f_1eb_0v_newelstypes[] =
{
  HP_PRISM_1FB_1EA_0V,
  HP_PRISM_1FA_0E_0V,
  HP_TET,
  HP_NONE,
};
int reftet_1f_1eb_0v_newels[][8] =
{
  // { 2, 5, 6, 3, 7, 8 },
  { 3, 8, 7, 2, 6, 5 },
  { 6, 4, 8, 5, 9, 7 },
  { 5, 9, 7, 1}
};
HPRef_Struct reftet_1f_1eb_0v =
{
  HP_TET,
  reftet_1f_1eb_0v_splitedges, 
  0, 0, 
  reftet_1f_1eb_0v_newelstypes, 
  reftet_1f_1eb_0v_newels
};







// HP_TET_1F_1E_1VA      // 1 sing edge in face e23, sing vert 2
int reftet_1f_1e_1va_splitedges[][3] =
{
  { 2, 1, 5 },
  { 2, 3, 10 }, 
  { 2, 4, 6 },
  { 3, 1, 7 },
  { 3, 4, 8 },
  { 4, 1, 9 },
  { 0, 0, 0 }
};

HPREF_ELEMENT_TYPE reftet_1f_1e_1va_newelstypes[] =
{
  HP_PRISM_1FB_1EA_0V,
  HP_PRISM_1FA_0E_0V,
  HP_TET,
  HP_TET_1F_1E_1VA,
  HP_NONE,
};
int reftet_1f_1e_1va_newels[][8] =
{
  { 3, 8, 7, 10, 6, 5 },
  { 6, 4, 8, 5, 9, 7 },
  { 5, 9, 7, 1},
  { 5, 2, 10, 6 }
};
HPRef_Struct reftet_1f_1e_1va =
{
  HP_TET,
  reftet_1f_1e_1va_splitedges, 
  0, 0, 
  reftet_1f_1e_1va_newelstypes, 
  reftet_1f_1e_1va_newels
};







// HP_TET_1F_1E_1VB      // 1 sing edge in face e24, sing vert 2
int reftet_1f_1e_1vb_splitedges[][3] =
{
  { 2, 1, 5 },
  { 2, 3, 6 },
  { 2, 4, 7 },
  { 3, 1, 8 },
  { 4, 1, 9 },
  { 4, 3, 10 },
  { 0, 0, 0 }
};

HPREF_ELEMENT_TYPE reftet_1f_1e_1vb_newelstypes[] =
{
  HP_PRISM_1FB_1EA_0V,
  HP_PRISM_1FA_0E_0V,  
  HP_TET,
  HP_TET_1F_1E_1VB,  
  HP_NONE,
};
int reftet_1f_1e_1vb_newels[][8] =
{
  { 4, 9, 10, 7, 5, 6 },
  { 3, 6, 10, 8, 5, 9 },
  { 5, 9, 8, 1},
  { 5, 2, 6, 7 }  
};
HPRef_Struct reftet_1f_1e_1vb =
{
  HP_TET,
  reftet_1f_1e_1vb_splitedges, 
  0, 0, 
  reftet_1f_1e_1vb_newelstypes, 
  reftet_1f_1e_1vb_newels
};












// HP_TET_1F_1E_2VA     //  1 sing edge not in face (e12), sing v2,v3    
int reftet_1f_1e_2va_splitedges[][3] =
{
  { 1, 3, 5 },
  { 1, 4, 6 },
  { 2, 1, 7 },
  { 2, 3, 8 },
  { 2, 4, 9 },
  { 3, 1, 10 },
  { 3, 2, 11 },
  { 3, 4, 12 },
  { 4, 1, 13 },
  { 0, 0, 0 }
};

int reftet_1f_1e_2va_splitfaces[][4] =
  {
    { 2, 1, 3, 14 },
    { 2, 1, 4, 15 },
    { 0, 0, 0, 0 }
  };


HPREF_ELEMENT_TYPE reftet_1f_1e_2va_newelstypes[] =
{
  HP_PRISM,
  HP_PRISM_SINGEDGE,
  HP_PRISM_1FB_0E_0V,
  HP_TET_1F_0E_1VA,
  HP_TET_1F_0E_1VA,
  HP_TET_1E_1VA,  
  HP_TET_1E_1VA,
  HP_HEX7_1FB,
  HP_NONE,  
};
int reftet_1f_1e_2va_newels[][8] =
{
  { 5, 14, 10, 6, 15, 13 },
  { 1, 5, 6, 7, 14, 15 },
  { 8, 14, 9, 11, 10, 12 },
  { 10, 3, 12, 11 },
  { 14, 2, 8, 9 },
  { 2, 7, 15, 14 },
  { 2, 9, 14, 15 },
  // { 13, 10, 14, 15, 4, 12, 9 }
  { 10, 13, 15, 14, 12, 4, 9 }
};
HPRef_Struct reftet_1f_1e_2va =
{
  HP_TET,
  reftet_1f_1e_2va_splitedges,
  reftet_1f_1e_2va_splitfaces,   
  0, 
  reftet_1f_1e_2va_newelstypes, 
  reftet_1f_1e_2va_newels
};





// HP_TET_1F_1E_2VB     //  1 sing edge not in face (e12), sing v2,v4    
int reftet_1f_1e_2vb_splitedges[][3] =
{
  { 1, 3, 5 },
  { 1, 4, 6 },
  { 2, 1, 7 },
  { 2, 3, 8 },
  { 2, 4, 9 },
  { 3, 1, 10 },
  { 4, 1, 11 },
  { 4, 2, 12 },
  { 4, 3, 13 },
  { 0, 0, 0 }
};

int reftet_1f_1e_2vb_splitfaces[][4] =
  {
    { 2, 1, 3, 14 },
    { 2, 1, 4, 15 },
    { 0, 0, 0, 0 }
  };


HPREF_ELEMENT_TYPE reftet_1f_1e_2vb_newelstypes[] =
{
  HP_PRISM,
  HP_PRISM_SINGEDGE,
  HP_PRISM_1FB_0E_0V,
  HP_TET_1F_0E_1VA,
  HP_TET_1F_0E_1VA,
  HP_TET_1E_1VA,  
  HP_TET_1E_1VA,
  HP_HEX7_1FB,
  HP_NONE,

};
int reftet_1f_1e_2vb_newels[][8] =
{
  { 5, 14, 10, 6, 15, 11 },
  { 1, 5, 6, 7, 14, 15 },
  { 8, 15, 9, 13, 11, 12 },
  { 11, 4, 12, 13 },
  { 15, 2, 8, 9 },
  { 2, 7, 15, 14 },
  { 2, 8, 14, 15 },
  { 10, 11, 15, 14, 3, 13, 8 }
};
HPRef_Struct reftet_1f_1e_2vb =
{
  HP_TET,
  reftet_1f_1e_2vb_splitedges,
  reftet_1f_1e_2vb_splitfaces,   
  0, 
  reftet_1f_1e_2vb_newelstypes, 
  reftet_1f_1e_2vb_newels
};






//  HP_TET_1F_1EA_3V
HPRefStruct<HP_TET> reftet_1f_1ea_3v
  {
    HP_TET_1F_1EA_3V,
    {
      El(HP_PRISM,  { E14, E41, F214, E13, E31, F213 }),
      El(HP_PRISM_SINGEDGE, { V1, E13, E14, E21, F213, F214 }),
      El(HP_HEX7_1FB, { E31, E41, F214, F213, F324, F423, F234 }),
      El(HP_PRISM_1FB_0E_0V, { F234, E23, F213, F324, E32, E31 }),
      El(HP_PRISM_1FB_0E_0V, { F423, E42, E41, F234, E24, F214 }),
      El(HP_PRISM_1FB_0E_0V, { F324, E34, E31, F423, E43, E41 }),
      El(HP_TET_1F_0E_0V, { E31, E32, E34, F324 }),
      El(HP_TET_1F_0E_1VA, { E31, V3, E34, E32 }),

      El(HP_TET_1F_0E_0V, { E41, E43, E42, F423 }),
      El(HP_TET_1F_0E_1VA, { E41, V4, E42, E43 }),
      El(HP_PYRAMID_1FB_0E_0V, {  E24, E23, F213, F214, F234 }),  // needs check
      El(HP_PYRAMID_1FB_0E_1VA, { E23, E24, F214, F213, V2 }),
      El(HP_TET_1E_1VA, { V2, E21, F214, F213 }),
    }
  };



//  HP_TET_1F_1E_3V     e4 (E23), V2, V3, V4
HPRefStruct<HP_TET> reftet_1f_1e_3v
  {
    HP_TET_1F_1E_3V,
    {
      El(HP_TET,  { V1, E21, E31, E41 }),
      El(HP_HEX7_1FA, { E34, E24, E42, E43, E31, E21, E41 }),
      El(HP_PRISM_1FB_1EA_0V, { E32, E34, E31, E23, E24, E21 }),
      El(HP_TET_1F_0E_1VA, { E41, V4, E42, E43 }),
      El(HP_TET_1F_1E_1VB, { E21, V2, E23, E24 }),
      El(HP_TET_1F_1E_1VA, { E31, V3, E34, E32 }),
    }
  };



HPRefStruct<HP_TET> reftet_1f_2eoo_3v
  {
    HP_TET_1F_2Eoo_3V,
    {
      El(HP_TET,  { E14, E41, F214, F314 }),
      El(HP_PRISM_1FA_0E_0V, { V4, E34, E24, E41, F314, F214 }),
      El(HP_PRISM, { F123, F213, F312, E14, F214, F314 }),
      El(HP_PRISM_SINGEDGE, { E12, F123, E14, E21, F213, F214 }), 
      El(HP_PRISM_SINGEDGE, { E13, E14, F123, E31, F314, F312 }),
      El(HP_HEX_1F_0E_0V, { E32, E23, E24, E34, F312, F213, F214, F314 }),
      El(HP_TET_1E_1VA, { V1, E12, F123, E14 }),
      El(HP_TET_1E_1VA, { V1, E13, E14, F123 }),
      El(HP_PYRAMID_1FB_0E_1VA, { E34, E32, F312, F314, V3 }),
      El(HP_PYRAMID_1FB_0E_1VA, { E23, E24, F214, F213, V2 }),
      El(HP_TET_1E_1VA, { V2, E21, F214, F213 }),
      El(HP_TET_1E_1VA, { V3, E31, F312, F314 }),

    }
  };

// HP_TET_1F_2E_0VA    singular edge in face 234 is 34, and edge not in face is 14
int reftet_1f_2e_0va_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 2, 1, 7 },
  { 3, 1, 8 },
  { 3, 2, 9 },
  { 4, 1, 10 },
  { 4, 2, 11 },
  { 4, 3, 12 },
  { 0, 0, 0 }
};

int reftet_1f_2e_0va_splitfaces[][4] =
{
  { 4, 1, 2, 13 }, 
  { 4, 1, 3, 14 },
  { 0, 0, 0, 0 }
};

HPREF_ELEMENT_TYPE reftet_1f_2e_0va_newelstypes[] =
{
  HP_PRISM,
  HP_PRISM_SINGEDGE,
  HP_HEX7_1FB,
  HP_PRISM_1FB_1EA_0V,
  HP_TET_1F_1E_1VA,
  HP_TET_1E_1VA,
  HP_TET_1E_1VA,
  HP_NONE,
};
int reftet_1f_2e_0va_newels[][8] =
{
  { 6, 8, 14, 5, 7, 13 },
  { 1, 5, 6, 10, 13, 14 },
  { 7, 8, 14, 13, 2, 9, 11 },
  { 3, 8, 9, 12, 14, 11 },
  { 14, 4, 11, 12 },
  { 4, 11, 13, 14 },
  { 4, 10, 14, 13 }
};
HPRef_Struct reftet_1f_2e_0va =
{
  HP_TET,
  reftet_1f_2e_0va_splitedges, 
  reftet_1f_2e_0va_splitfaces, 
  0, 
  reftet_1f_2e_0va_newelstypes, 
  reftet_1f_2e_0va_newels
};





// HP_TET_1F_2E_0VB    singular edge in face 234 is 34, and edge not in face is 13
int reftet_1f_2e_0vb_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 4, 6 },
  { 2, 1, 7 },
  { 3, 1, 8 },
  { 3, 2, 9 },
  { 3, 4, 10 },
  { 4, 1, 11 },
  { 4, 2, 12 },
  { 0, 0, 0 }
};

int reftet_1f_2e_0vb_splitfaces[][4] =
{
  { 3, 1, 2, 13 }, 
  { 3, 1, 4, 14 },
  { 0, 0, 0, 0 }
};

HPREF_ELEMENT_TYPE reftet_1f_2e_0vb_newelstypes[] =
{
  HP_PRISM,
  HP_PRISM_SINGEDGE,
  HP_HEX7_1FB,
  HP_PRISM_1FB_1EA_0V,
  HP_TET_1F_1E_1VA,
  HP_TET_1E_1VA,
  HP_TET_1E_1VA,
  HP_NONE,
};
int reftet_1f_2e_0vb_newels[][8] =
{
  { 6, 14, 11, 5, 13, 7 },
  { 1, 6, 5, 8, 14, 13 },
  { 11, 7, 13, 14, 12, 2, 9 },
  { 4, 12, 11, 10, 9, 14 },
  { 14, 3, 10, 9 },
  { 3, 8, 13, 14  },
  { 3, 9, 14, 13 }
};
HPRef_Struct reftet_1f_2e_0vb =
{
  HP_TET,
  reftet_1f_2e_0vb_splitedges, 
  reftet_1f_2e_0vb_splitfaces, 
  0, 
  reftet_1f_2e_0vb_newelstypes, 
  reftet_1f_2e_0vb_newels
};





//  HP_TET_1F_2E_1V     e4,e5 (E23,E24), V2 
HPRefStruct<HP_TET> reftet_1f_2e_1v
  {
    HP_TET_1F_2E_1V,
    {
      El(HP_TET,  { V1, E21, E31, E41 }),
      El(HP_PRISM_1FA_0E_0V, { F234, E43, E34, E21, E41, E31 }),
      El(HP_PRISM_1FB_1EA_0V, { V3, E34, E31, E23, F234, E21 }),
      El(HP_PRISM_1FB_1EA_0V, { E24, F234, E21, V4, E43, E41 }),
      El(HP_TET_1F_1E_1VA, { E21, V2, E23, F234 }),
      El(HP_TET_1F_1E_1VB, { E21, V2, F234, E24 }),
    }
  };



//  HP_TET_1F_2E_3V     e4,e5 (E23,E24), V2, V3, V4
HPRefStruct<HP_TET> reftet_1f_2e_3v
  {
    HP_TET_1F_2E_3V,
    {
      El(HP_TET,  { V1, E21, E31, E41 }),
      El(HP_PRISM_1FA_0E_0V, { F234, E43, E34, E21, E41, E31 }),
      El(HP_PRISM_1FB_1EA_0V, { E32, E34, E31, E23, F234, E21 }),
      El(HP_PRISM_1FB_1EA_0V, { E24, F234, E21, E42, E43, E41 }),
      El(HP_TET_1F_1E_1VA, { E21, V2, E23, F234 }),
      El(HP_TET_1F_1E_1VB, { E21, V2, F234, E24 }),
      El(HP_TET_1F_1E_1VA, { E31, V3, E34, E32 }),
      El(HP_TET_1F_1E_1VB, { E41, V4, E42, E43 }),
    }
  };





/* ************************ 2 singular faces ******************** */


// HP_TET_2F_0E_0V
int reftet_2f_0e_0v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 2, 1, 6 },
  { 3, 1, 7 },
  { 3, 2, 8 },
  { 4, 1, 9 },
  { 4, 2, 10 },
  { 0, 0, 0 }
};

int reftet_2f_0e_0v_splitfaces[][4] =
  {
    { 3, 1, 2, 11 },
    { 4, 1, 2, 12 },
    { 0, 0, 0, 0 }
  };


HPREF_ELEMENT_TYPE reftet_2f_0e_0v_newelstypes[] =
{
  HP_PRISM_1FA_0E_0V,
  HP_PRISM_1FA_0E_0V,
  HP_PRISM_1FB_1EA_0V,
  HP_PRISM_1FB_1EA_0V,
  HP_TET,
  HP_NONE,
};
int reftet_2f_0e_0v_newels[][8] =
{
  { 2, 10, 8, 6, 12, 11 },
  { 1, 7, 9, 5, 11, 12 },
  //   { 3, 11, 8, 4, 12, 10 },
  { 4, 10, 12, 3, 8, 11 }, 
  { 3, 7, 11, 4, 9, 12 },
  { 5, 6, 11, 12 }
};
HPRef_Struct reftet_2f_0e_0v =
{
  HP_TET,
  reftet_2f_0e_0v_splitedges, 
  reftet_2f_0e_0v_splitfaces, 
  0, 
  reftet_2f_0e_0v_newelstypes, 
  reftet_2f_0e_0v_newels
};





// HP_TET_2F_0E_1V
int reftet_2f_0e_1v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 2, 1, 6 },
  { 3, 1, 7 },
  { 3, 2, 8 },
  { 4, 1, 9 },
  { 4, 2, 10 },
  { 4, 3, 13 },
  { 0, 0, 0 }
};

int reftet_2f_0e_1v_splitfaces[][4] =
  {
    { 3, 1, 2, 11 },
    { 4, 1, 2, 12 },
    { 0, 0, 0, 0 }
  };


HPREF_ELEMENT_TYPE reftet_2f_0e_1v_newelstypes[] =
{
  HP_PRISM_1FA_0E_0V,
  HP_PRISM_1FA_0E_0V,
  HP_PRISM_1FB_1EA_0V,
  HP_PRISM_1FB_1EA_0V,
  HP_TET,
  HP_TET_1F_1E_1VB,  
  HP_TET_1F_1E_1VA,  
  HP_NONE
};
int reftet_2f_0e_1v_newels[][8] =
{
  { 2, 10, 8, 6, 12, 11 },
  { 1, 7, 9, 5, 11, 12 },
  //   { 3, 11, 8, 4, 12, 10 },
  { 13, 10, 12, 3, 8, 11 }, 
  { 3, 7, 11, 13, 9, 12 },
  { 5, 6, 11, 12 },
  { 12, 4, 10, 13 },
  { 12, 4, 13, 9 } 
};
HPRef_Struct reftet_2f_0e_1v =
{
  HP_TET,
  reftet_2f_0e_1v_splitedges, 
  reftet_2f_0e_1v_splitfaces, 
  0, 
  reftet_2f_0e_1v_newelstypes, 
  reftet_2f_0e_1v_newels
};








//  HP_TET_2F_1E_0VA,  // 2 singular faces, sing edge e4
int reftet_2f_1e_0va_splitedges[][3] =
{
  { 1, 2, 5 },
  { 2, 1, 6 },
  { 2, 4, 7 },
  { 3, 1, 8 },
  { 3, 2, 9 },
  { 3, 4, 10 },
  { 4, 1, 11 },
  { 4, 2, 12 },
  { 0, 0, 0 }
};

int reftet_2f_1e_0va_splitfaces[][4] =
  {
    { 3, 2, 1, 13 },
    { 3, 2, 4, 14 },
    { 4, 1, 2, 15 },
    { 0, 0, 0, 0 }
  };

HPREF_ELEMENT_TYPE reftet_2f_1e_0va_newelstypes[] =
{
  HP_TET,
  HP_PRISM_1FA_0E_0V,
  HP_PRISM_1FA_0E_0V,
  HP_PRISM_SINGEDGE, 
  HP_PRISM_SINGEDGE, 
  HP_PRISM_SINGEDGE,
  HP_TET_1F_1E_1VA, 
  HP_TET_1F_1E_1VB, 
  HP_TET_1F_1E_1VB, 
  HP_NONE,
};
int reftet_2f_1e_0va_newels[][8] =
{
  { 5, 6, 13, 15 },
  { 1, 8, 11, 5, 13, 15 },
  { 7, 12, 14, 6, 15, 13 },
  { 2, 6, 7, 9, 13, 14 },
  { 4, 15, 11,  10, 13, 8 },
  { 4, 12, 15,  10, 14, 13, },
  { 13, 3, 10, 14 },
  { 13, 3, 14, 9 },
  { 13, 3, 8, 10 },
};
HPRef_Struct reftet_2f_1e_0va =
{
  HP_TET,
  reftet_2f_1e_0va_splitedges, 
  reftet_2f_1e_0va_splitfaces, 
  0, 
  reftet_2f_1e_0va_newelstypes, 
  reftet_2f_1e_0va_newels
};


//  HP_TET_2F_1E_0VB = 602,  // 2 singular faces f234,f134, sing edge e5=e24
int reftet_2f_1e_0vb_splitedges[][3] =
{
  { 1, 2, 5 },
  { 2, 1, 6 },
  { 2, 3, 7 },
  { 3, 1, 8 },
  { 3, 2, 9 },
  { 4, 1, 10 },
  { 4, 2, 11 },
  { 4, 3, 12 },
  { 0, 0, 0 }
};

int reftet_2f_1e_0vb_splitfaces[][4] =
  {
    { 4, 2, 3, 13 },
    { 4, 1, 2, 14 },
    { 3, 1, 2, 15 },
    { 0, 0, 0, 0 }
  };

HPREF_ELEMENT_TYPE reftet_2f_1e_0vb_newelstypes[] =
{
  HP_TET,
  HP_PRISM_1FA_0E_0V,
  HP_PRISM_1FA_0E_0V,
  HP_PRISM_SINGEDGE, 
  HP_PRISM_SINGEDGE, 
  HP_PRISM_SINGEDGE,
  HP_TET_1F_1E_1VA, 
  HP_TET_1F_1E_1VB, 
  HP_TET_1F_1E_1VB, 
  HP_NONE,
};
int reftet_2f_1e_0vb_newels[][8] =
{
  { 5, 6, 15, 14 },
  { 1, 8, 10, 5, 15, 14 },
  { 7, 13, 9, 6, 14, 15 },
  { 2, 7, 6, 11, 13, 14 },
  { 3, 8, 15, 12, 10, 14 },
  { 3, 15, 9, 12, 14, 13 },
  { 14, 4, 11, 13 },
  { 14, 4, 13, 12 },
  { 14, 4, 12, 10 },
};

HPRef_Struct reftet_2f_1e_0vb =
{
  HP_TET,
  reftet_2f_1e_0vb_splitedges, 
  reftet_2f_1e_0vb_splitfaces, 
  0, 
  reftet_2f_1e_0vb_newelstypes, 
  reftet_2f_1e_0vb_newels
};


// HP_TET_3F_0E_0V = 610,  // 3 singular faces, no additional points or edges
int reftet_3f_0e_0v_splitedges[][3] =
{
  { 1, 2, 5 },
  { 1, 3, 6 },
  { 2, 1, 7 },
  { 2, 3, 8 },
  { 3, 1, 9 },
  { 3, 2, 10 },
  { 4, 1, 11 },
  { 4, 2, 12 },
  { 4, 3, 13 },
  { 0, 0, 0 }
};

int reftet_3f_0e_0v_splitfaces[][4] =
  {
    { 1, 2, 3, 14 },
    { 2, 3, 1, 15 },
    { 3, 1, 2, 16 },
    { 4, 1, 2, 17 },
    { 4, 1, 3, 18 },
    { 4, 2, 3, 19 },
    { 0, 0, 0, 0 }
  };

int reftet_3f_0e_0v_splitelements[][5] =
  {
    { 4, 1, 2, 3, 20 },
    { 0 }
  };
HPREF_ELEMENT_TYPE reftet_3f_0e_0v_newelstypes[] =
{
  HP_TET,
  HP_PRISM_1FA_0E_0V,
  HP_PRISM_1FA_0E_0V,
  HP_PRISM_1FA_0E_0V,
  
  HP_PRISM_1FB_1EA_0V,
  HP_PRISM_1FB_1EA_0V,
  HP_PRISM_1FB_1EA_0V,
  HP_PRISM_1FB_1EA_0V,
  HP_PRISM_1FB_1EA_0V,
  HP_PRISM_1FB_1EA_0V,

  HP_TET_1F_1E_1VA,  // 1E_1VA
  HP_TET_1F_0E_0V,  
  HP_TET_1F_0E_0V,  // 1E_1VA
  HP_TET_1F_0E_0V,  // 1E_1VA
  HP_TET_1F_0E_0V,  // 1E_1VA
  HP_TET_1F_0E_0V,  // 1E_1VA
  HP_NONE,
};
int reftet_3f_0e_0v_newels[][8] =
{
  { 14, 15, 16, 20 },
  { 5, 17, 7, 14, 20, 15 },
  { 6, 9, 18, 14, 16, 20 },
  { 10, 8, 19, 16, 15, 20 },
  
  { 1, 5, 14, 11, 17, 20 },
  { 11, 18, 20, 1, 6, 14 },
  { 2, 8, 15, 12, 19, 20 },
  { 12, 17, 20, 2, 7, 15 },
  { 3, 9, 16, 13, 18, 20 },
  { 13, 19, 20, 3, 10, 16 },

  { 20, 4, 11, 17 },
  { 20, 4, 17, 12 },
  { 20, 4, 12, 19 },
  { 20, 4, 19, 13 },
  { 20, 4, 13, 18 },
  { 20, 4, 18, 11 }
};
HPRef_Struct reftet_3f_0e_0v =
{
  HP_TET,
  reftet_3f_0e_0v_splitedges, 
  reftet_3f_0e_0v_splitfaces, 
  reftet_3f_0e_0v_splitelements, 
  reftet_3f_0e_0v_newelstypes, 
  reftet_3f_0e_0v_newels
};







/*

*/
