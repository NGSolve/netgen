#ifndef WRITEUSER
#define WRITEUSER

/**************************************************************************/
/* File:    writeuser.hh                                                  */
/* Authors: many                                                          */
/* Date:    10. Dec. 97                                                   */
/**************************************************************************/

namespace netgen {

DLL_HEADER extern
void WriteFile (int typ,
                const Mesh & mesh,
                const NetgenGeometry & geom,
                const char * filename,
                const char * geomfile = NULL,
                double h = 0);



DLL_HEADER extern
void ReadFile (Mesh & mesh,
               const string & filename);






extern
void WriteNeutralFormat (const Mesh & mesh,
                         const NetgenGeometry & geom,
                         const string & filename);

extern
void WriteSurfaceFormat (const Mesh & mesh,
                         const string & filename);

extern
void WriteSTLFormat (const Mesh & mesh,
                     const string & filename);


// Philippose - 16 August 2010
// Added the STL Extended format in which
// each face of the geometry is treated as
// a separate "solid" entity in the STL file
extern
void WriteSTLExtFormat (const Mesh & mesh,
                        const string & filename);


extern
void WriteVRMLFormat (const Mesh & mesh,
                      bool faces,
                      const string & filename);

extern
void WriteFEPPFormat (const Mesh & mesh,
                      const NetgenGeometry & geom,
                      const string & filename);

extern
void WriteGmshFormat (const Mesh & mesh,
                      const NetgenGeometry & geom,
                      const string & filename);


// Philippose - 29/01/2009
// Added GMSH v2.xx Mesh Export support
void WriteGmsh2Format (const Mesh & mesh,
                       const NetgenGeometry & geom,
                       const string & filename);


// Philippose - 25/10/2009
// Added OpenFOAM 1.5+ Mesh Export support
extern 
void WriteOpenFOAM15xFormat (const Mesh & mesh, 
                             const string & casename,
                             const bool compressed);


extern
void WriteUserChemnitz (const Mesh & mesh,
                        const string & filename);

extern
void WriteJCMFormat (const Mesh & mesh,
                     const NetgenGeometry & geom,
                     const string & filename);


extern
void WriteDiffPackFormat (const Mesh & mesh,
                          const NetgenGeometry & geom,
                          const string & filename);

extern
void WriteTochnogFormat (const Mesh & mesh,
                         const string & filename);

extern
void WriteTecPlotFormat (const Mesh & mesh,
                         const NetgenGeometry & geom,
                         const string & filename);

extern
void WriteAbaqusFormat (const Mesh & mesh,
                        const string & filename,
                        bool forcesurface,
                        bool forcevolume);

extern
void WriteFluentFormat (const Mesh & mesh,
                        const string & filename);

extern
void WritePermasFormat (const Mesh & mesh,
                        const string & filename);

extern
void WriteFEAPFormat (const Mesh & mesh,
                      const string & filename);

extern
void WriteElmerFormat (const Mesh & mesh,
                       const string & filename);


extern
void WriteEdgeElementFormat (const Mesh & mesh,
                             const NetgenGeometry & geom,
                             const string & filename);



#ifdef OLIVER
extern
void WriteTETFormat (const Mesh & mesh,
                     const string & filename);

#endif

extern void ReadTETFormat (Mesh & mesh,
                           const string & filename);


extern void ReadFNFFormat (Mesh & mesh,
                           const string & filename);



extern void DLL_HEADER ReadCGNSMesh (Mesh & mesh,
                           const string & filename);

// Read Mesh and solutions from CGNS file
extern tuple<shared_ptr<Mesh>, vector<string>, vector<Array<double>>, vector<int>>
DLL_HEADER ReadCGNSFile(string filename, int base);



void WriteDolfinFormat (const Mesh & mesh,
                        const string & filename);


extern void DLL_HEADER RegisterUserFormats (NgArray<const char*> & names,
                                 NgArray<const char*> & extensions);


extern bool DLL_HEADER WriteUserFormat (const string & format,
                                        const Mesh & mesh,
                                        // const NetgenGeometry & geom,
                                        const string & filename);

}

#endif

