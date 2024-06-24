NGX_INLINE DLL_HEADER Ng_Point Ngx_Mesh :: GetPoint (int nr) const
{
  return Ng_Point (&mesh->Point(PointIndex(nr+PointIndex::BASE))(0));
}


template <>
NGX_INLINE DLL_HEADER int Ngx_Mesh :: GetElementIndex<0> (size_t nr) const
{
  return (*mesh).pointelements[nr].index;
}

template <>
NGX_INLINE DLL_HEADER int Ngx_Mesh :: GetElementIndex<1> (size_t nr) const
{
  /*
  if(mesh->GetDimension()==3)
    return (*mesh)[SegmentIndex(nr)].edgenr;
  else
    return (*mesh)[SegmentIndex(nr)].si;
  */
  if(mesh->GetDimension()==3)
    return mesh->LineSegments()[nr].edgenr;
  else
    return mesh->LineSegments()[nr].si;    
}
  
template <>
NGX_INLINE DLL_HEADER int Ngx_Mesh :: GetElementIndex<2> (size_t nr) const
{
  // int ind = (*mesh)[SurfaceElementIndex(nr)].GetIndex(); 
  // return mesh->GetFaceDescriptor(ind).BCProperty();
  const Element2d & el = (*mesh)[SurfaceElementIndex(nr)];
  return mesh->GetFaceDescriptor(el).BCProperty();
}

template <>
NGX_INLINE DLL_HEADER int Ngx_Mesh :: GetElementIndex<3> (size_t nr) const
{
  return (*mesh)[ElementIndex(nr)].GetIndex();
}


template <>
NGX_INLINE DLL_HEADER Ng_Element Ngx_Mesh :: GetElement<0> (size_t nr) const
{
  const Element0d & el = mesh->pointelements[nr];
  
  Ng_Element ret;
  ret.type = NG_PNT;
  ret.index = el.index;
  ret.mat = el.name;
  
  ret.points.num = 1;
  ret.points.ptr = (int*)&el.pnum;
  
  ret.vertices.num = 1;
  ret.vertices.ptr = (int*)&el.pnum;

  /*
  ret.edges.num = 0;
  ret.edges.ptr = NULL;
  */
  ret.edges.Assign ( FlatArray<T_EDGE2> (0, nullptr) );
  /*
  ret.faces.num = 0;
  ret.faces.ptr = NULL;
  */
  ret.faces.Assign ( { 0, nullptr } );  
  
  ret.facets.num = 1;
  ret.facets.base = POINTINDEX_BASE;
  ret.facets.ptr = (int*)&el.pnum;

  if (mesh->GetDimension() == 1)
    ret.mat = *(mesh->GetBCNamePtr(el.index-1));
  else if (mesh->GetDimension() == 2)
    ret.mat = *(mesh->GetCD2NamePtr(el.index-1));
  else
    ret.mat = *(mesh->GetCD3NamePtr(el.index-1));
  
  return ret;
}



template <> 
NGX_INLINE DLL_HEADER Ng_Element Ngx_Mesh :: GetElement<1> (size_t nr) const
{
  // const Segment & el = mesh->LineSegment (SegmentIndex(nr));
  const Segment & el = mesh->LineSegments()[nr];

  Ng_Element ret;
  ret.type = NG_ELEMENT_TYPE(el.GetType());
  if(mesh->GetDimension()==3)
    ret.index = el.edgenr;
  else
    ret.index = el.si;
  if (mesh->GetDimension() == 2)
    ret.mat = *(mesh->GetBCNamePtr(el.si-1));
  else
    {
      if (mesh->GetDimension() == 3)
        ret.mat = *(mesh->GetCD2NamePtr(el.edgenr-1));
      else
        ret.mat = *(mesh->GetMaterialPtr(el.si));
    }

  ret.points.num = el.GetNP();
  ret.points.ptr = (int*)&(el[0]);

  ret.vertices.num = 2;
  ret.vertices.ptr = (int*)&(el[0]);

  /*
  ret.edges.num = 1;
  ret.edges.ptr = mesh->GetTopology().GetSegmentElementEdgesPtr (nr);
  */
  ret.edges.Assign ( FlatArray<T_EDGE2> (1, const_cast<T_EDGE2*>( mesh->GetTopology().GetSegmentElementEdgesPtr (nr))));

  /*
  ret.faces.num = 0;
  ret.faces.ptr = NULL;
  */
  ret.faces.Assign ( { 0, nullptr });
  
  if (mesh->GetDimension() == 3)
    {
      ret.facets.num = 0;
      ret.facets.base = 0;
      ret.facets.ptr = nullptr;
    }
  else if (mesh->GetDimension() == 2)
    {
      ret.facets.num = 1;
      ret.facets.base = 0;
      ret.facets.ptr = ret.edges.Data();
    }
  else
    {
      ret.facets.num = 2;
      ret.facets.base = 1;
      ret.facets.ptr = (int*)&(el[0]);
    }

  // ret.is_curved = mesh->GetCurvedElements().IsSegmentCurved(nr);
  ret.is_curved = el.IsCurved();

  return ret;
}

template <> 
NGX_INLINE DLL_HEADER Ng_Element Ngx_Mesh :: GetElement<2> (size_t nr) const
{
  const Element2d & el = mesh->SurfaceElements()[nr];
  
  Ng_Element ret;
  ret.type = NG_ELEMENT_TYPE(el.GetType());
  const FaceDescriptor & fd = mesh->GetFaceDescriptor(el); // .GetIndex());
  ret.index = fd.BCProperty();
  if (mesh->GetDimension() == 3)
    ret.mat = fd.GetBCName();
  else
    ret.mat = *(mesh -> GetMaterialPtr(ret.index));
  ret.points.num = el.GetNP();
  ret.points.ptr  = (int*)&el[0];

  ret.vertices.num = el.GetNV();
  ret.vertices.ptr = (int*)&(el[0]);

  /*
  ret.edges.num = MeshTopology::GetNEdges (el.GetType());
  ret.edges.ptr = mesh->GetTopology().GetSurfaceElementEdgesPtr (nr);
  */
  ret.edges.Assign (mesh->GetTopology().GetEdges (SurfaceElementIndex(nr)));
  /*
  ret.faces.num = MeshTopology::GetNFaces (el.GetType());
  ret.faces.ptr = mesh->GetTopology().GetSurfaceElementFacesPtr (nr);
  */
  ret.faces.Assign ( { 1, const_cast<int*>(mesh->GetTopology().GetSurfaceElementFacesPtr (nr)) });
  if (mesh->GetDimension() == 3)
    {
      ret.facets.num = ret.faces.Size();
      ret.facets.base = 0;
      ret.facets.ptr = ret.faces.Data();
    }
  else
    {
      ret.facets.num = ret.edges.Size();
      ret.facets.base = 0;      
      ret.facets.ptr = ret.edges.Data();
    }
  ret.is_curved = el.IsCurved();
  ret.newest_vertex = el.NewestVertex();
  return ret;
}

template <> 
NGX_INLINE DLL_HEADER Ng_Element Ngx_Mesh :: GetElement<3> (size_t nr) const
{
  const Element & el = mesh->VolumeElements()[nr];
  
  Ng_Element ret;
  ret.type = NG_ELEMENT_TYPE(el.GetType());
  ret.index = el.GetIndex();
  ret.mat = *(mesh -> GetMaterialPtr(ret.index));
  ret.points.num = el.GetNP();
  ret.points.ptr = (int*)&el[0];

  ret.vertices.num = el.GetNV();
  ret.vertices.ptr = (int*)&(el[0]);

  /*
  ret.edges.num = MeshTopology::GetNEdges (el.GetType());
  ret.edges.ptr = mesh->GetTopology().GetElementEdgesPtr (nr);
  */
  ret.edges.Assign (mesh->GetTopology().GetEdges (ElementIndex(nr)));  

  /*
  ret.faces.num = MeshTopology::GetNFaces (el.GetType());
  ret.faces.ptr = mesh->GetTopology().GetElementFacesPtr (nr);
  */
  ret.faces.Assign (mesh->GetTopology().GetFaces (ElementIndex(nr)));
  
  ret.facets.num = ret.faces.Size();
  ret.facets.base = 0;
  ret.facets.ptr = ret.faces.Data();

  ret.is_curved = el.IsCurved();
  ret.newest_vertex = el.NewestVertex();
  return ret;
}



template <> NGX_INLINE DLL_HEADER
string_view Ngx_Mesh :: GetMaterialCD<0> (int region_nr) const
{
  return mesh->GetMaterial(region_nr+1);
}

template <> NGX_INLINE DLL_HEADER
string_view Ngx_Mesh :: GetMaterialCD<1> (int region_nr) const
{
  return mesh->GetBCName(region_nr);
}

template <> NGX_INLINE DLL_HEADER
string_view Ngx_Mesh :: GetMaterialCD<2> (int region_nr) const
{
  return mesh->GetCD2Name(region_nr);
}

template <> NGX_INLINE DLL_HEADER
string_view Ngx_Mesh :: GetMaterialCD<3> (int region_nr) const
{
  return mesh->GetCD3Name(region_nr);
}





template <> NGX_INLINE DLL_HEADER int Ngx_Mesh :: GetNNodes<1> ()
{
  return mesh->GetTopology().GetNEdges();
}

template <> NGX_INLINE DLL_HEADER int Ngx_Mesh :: GetNNodes<2> ()
{
  return mesh->GetTopology().GetNFaces();
}

template <> NGX_INLINE DLL_HEADER const Ng_Node<0> Ngx_Mesh :: GetNode<0> (int vnr) const
{
  Ng_Node<0> node;
  vnr++;
  switch (mesh->GetDimension())
    {
    case 3:
      {
        auto ia = mesh->GetTopology().GetVertexElements(vnr);
        node.elements.ne = ia.Size();
        node.elements.ptr = (int*)ia.Data();
        
        auto bia = mesh->GetTopology().GetVertexSurfaceElements(vnr);
        node.bnd_elements.ne = bia.Size();
        node.bnd_elements.ptr = (int*)bia.Data();
        break;
      }
    case 2:
      {
        auto ia = mesh->GetTopology().GetVertexSurfaceElements(vnr);
        node.elements.ne = ia.Size();
        node.elements.ptr = (int*)ia.Data();
        
        auto bia = mesh->GetTopology().GetVertexSegments(vnr);
        node.bnd_elements.ne = bia.Size();
        node.bnd_elements.ptr = (int*)bia.Data();
        break;
      }
    case 1:
      {
        auto ia = mesh->GetTopology().GetVertexSegments(vnr);
        node.elements.ne = ia.Size();
        node.elements.ptr = (int*)ia.Data();
        
        auto bia = mesh->GetTopology().GetVertexPointElements(vnr);
        node.bnd_elements.ne = bia.Size();
        node.bnd_elements.ptr = (int*)bia.Data();
        break;
      }
    default:
      ;
    }
  return node;
}
  
template <> NGX_INLINE DLL_HEADER const Ng_Node<1> Ngx_Mesh :: GetNode<1> (int nr) const
{
  Ng_Node<1> node;
  node.vertices.ptr = (const int*)mesh->GetTopology().GetEdgeVerticesPtr(nr);
  return node;
}

template <> NGX_INLINE DLL_HEADER const Ng_Node<2> Ngx_Mesh :: GetNode<2> (int nr) const
{
  Ng_Node<2> node;
  node.vertices.ptr = (const int*)mesh->GetTopology().GetFaceVerticesPtr(nr);
  node.vertices.nv = (node.vertices.ptr[3] == 0) ? 3 : 4;
  node.surface_el = mesh->GetTopology().GetFace2SurfaceElement (nr+1)-1;
  return node;
}


NGX_INLINE DLL_HEADER Ng_Buffer<int[2]> Ngx_Mesh :: GetPeriodicVertices(int idnr) const
{
  NgArray<INDEX_2> apairs;
  mesh->GetIdentifications().GetPairs (idnr+1, apairs);
  for(auto& ind : apairs)
    {
      ind.I1()--;
      ind.I2()--;
    }
  typedef int ti2[2];
  return { apairs.Size(), (ti2*)(void*)apairs.Release() };
}


NGX_INLINE void Ngx_Mesh :: GetParentNodes (int ni, int * parents) const
{
  ni++;
  if (ni <= mesh->mlbetweennodes.Size())
    {
      parents[0] = mesh->mlbetweennodes.Get(ni).I1()-1;
      parents[1] = mesh->mlbetweennodes.Get(ni).I2()-1;
      }
  else
    parents[0] = parents[1] = -1;
}

inline bool Ngx_Mesh :: HasParentEdges() const
{
  return mesh->GetTopology().HasParentEdges();
}

inline tuple<int, std::array<int,3>> Ngx_Mesh :: GetParentEdges (int enr) const
{
  return mesh->GetTopology().GetParentEdges(enr);
}

inline tuple<int, std::array<int,4>> Ngx_Mesh :: GetParentFaces (int fnr) const
{
  return mesh->GetTopology().GetParentFaces(fnr);
}


inline auto Ngx_Mesh :: GetTimeStamp() const { return mesh->GetTimeStamp(); }
