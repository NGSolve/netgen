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
  if(mesh->GetDimension()==3)
    return (*mesh)[SegmentIndex(nr)].edgenr;
  else
    return (*mesh)[SegmentIndex(nr)].si;
}
  
template <>
NGX_INLINE DLL_HEADER int Ngx_Mesh :: GetElementIndex<2> (size_t nr) const
{
  int ind = (*mesh)[SurfaceElementIndex(nr)].GetIndex(); 
  return mesh->GetFaceDescriptor(ind).BCProperty();
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
  
  ret.points.num = 1;
  ret.points.ptr = (int*)&el.pnum;
  
  ret.vertices.num = 1;
  ret.vertices.ptr = (int*)&el.pnum;
  
  ret.edges.num = 0;
  ret.edges.ptr = NULL;
  
  ret.faces.num = 0;
  ret.faces.ptr = NULL;

  ret.facets.num = 1;
  ret.facets.base = 1;
  ret.facets.ptr = (int*)&el.pnum;
  
  return ret;
}



template <> 
NGX_INLINE DLL_HEADER Ng_Element Ngx_Mesh :: GetElement<1> (size_t nr) const
{
  const Segment & el = mesh->LineSegment (SegmentIndex(nr));

  Ng_Element ret;
  ret.type = NG_ELEMENT_TYPE(el.GetType());
  if(mesh->GetDimension()==3)
    ret.index = el.edgenr;
  else
    ret.index = el.si;
  if (mesh->GetDimension() == 2)
    ret.mat = mesh->GetBCNamePtr(el.si-1);
  else
    {
      if (mesh->GetDimension() == 3)
        ret.mat = mesh->GetCD2NamePtr(el.edgenr-1);
      else
        ret.mat = mesh->GetMaterialPtr(el.si);
    }

  ret.points.num = el.GetNP();
  ret.points.ptr = (int*)&(el[0]);

  ret.vertices.num = 2;
  ret.vertices.ptr = (int*)&(el[0]);

  ret.edges.num = 1;
  ret.edges.ptr = (T_EDGE2*)mesh->GetTopology().GetSegmentElementEdgesPtr (nr);

  ret.faces.num = 0;
  ret.faces.ptr = NULL;

  if (mesh->GetDimension() == 2)
    {
      ret.facets.num = 1;
      ret.facets.base = 0;
      ret.facets.ptr = (int*)ret.edges.ptr;
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
  // const Element2d & el = mesh->SurfaceElement (SurfaceElementIndex (nr));
  const Element2d & el = mesh->SurfaceElements()[nr];
  
  Ng_Element ret;
  ret.type = NG_ELEMENT_TYPE(el.GetType());
  const FaceDescriptor & fd = mesh->GetFaceDescriptor(el.GetIndex());
  ret.index = fd.BCProperty();
  if (mesh->GetDimension() == 3)
    ret.mat = &fd.GetBCName();
  else
    ret.mat = mesh -> GetMaterialPtr(ret.index);
  ret.points.num = el.GetNP();
  ret.points.ptr  = (int*)&el[0];

  ret.vertices.num = el.GetNV();
  ret.vertices.ptr = (int*)&(el[0]);

  ret.edges.num = MeshTopology::GetNEdges (el.GetType());
  ret.edges.ptr = (T_EDGE2*)mesh->GetTopology().GetSurfaceElementEdgesPtr (nr);

  ret.faces.num = MeshTopology::GetNFaces (el.GetType());
  ret.faces.ptr = (T_FACE2*)mesh->GetTopology().GetSurfaceElementFacesPtr (nr);

  if (mesh->GetDimension() == 3)
    {
      ret.facets.num = ret.faces.num;
      ret.facets.base = 0;
      ret.facets.ptr = (int*)ret.faces.ptr;
    }
  else
    {
      ret.facets.num = ret.edges.num;
      ret.facets.base = 0;      
      ret.facets.ptr = (int*)ret.edges.ptr;
    }
  ret.is_curved = el.IsCurved();
  return ret;
}

template <> 
NGX_INLINE DLL_HEADER Ng_Element Ngx_Mesh :: GetElement<3> (size_t nr) const
{
  // const Element & el = mesh->VolumeElement (ElementIndex (nr));
  const Element & el = mesh->VolumeElements()[nr];
  
  Ng_Element ret;
  ret.type = NG_ELEMENT_TYPE(el.GetType());
  ret.index = el.GetIndex();
  ret.mat = mesh -> GetMaterialPtr(ret.index);
  ret.points.num = el.GetNP();
  ret.points.ptr = (int*)&el[0];

  ret.vertices.num = el.GetNV();
  ret.vertices.ptr = (int*)&(el[0]);

  ret.edges.num = MeshTopology::GetNEdges (el.GetType());
  ret.edges.ptr = (T_EDGE2*)mesh->GetTopology().GetElementEdgesPtr (nr);

  ret.faces.num = MeshTopology::GetNFaces (el.GetType());
  ret.faces.ptr = (T_FACE2*)mesh->GetTopology().GetElementFacesPtr (nr);

  ret.facets.num = ret.faces.num;
  ret.facets.base = 0;
  ret.facets.ptr = (int*)ret.faces.ptr;

  ret.is_curved = el.IsCurved();
  return ret;
}



template <> NGX_INLINE DLL_HEADER
const string &  Ngx_Mesh :: GetMaterialCD<0> (int region_nr) const
{
  return mesh->GetMaterial(region_nr+1);
}

template <> NGX_INLINE DLL_HEADER
const string &  Ngx_Mesh :: GetMaterialCD<1> (int region_nr) const
{
  return mesh->GetBCName(region_nr);
}

template <> NGX_INLINE DLL_HEADER
const string &  Ngx_Mesh :: GetMaterialCD<2> (int region_nr) const
{
  return mesh->GetCD2Name(region_nr);
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
        FlatArray<ElementIndex> ia = mesh->GetTopology().GetVertexElements(vnr);
        node.elements.ne = ia.Size();
        node.elements.ptr = (int*)&ia[0];
        
        FlatArray<SurfaceElementIndex> bia = mesh->GetTopology().GetVertexSurfaceElements(vnr);
        node.bnd_elements.ne = bia.Size();
        node.bnd_elements.ptr = (int*)&bia[0];
        break;
      }
    case 2:
      {
        FlatArray<SurfaceElementIndex> ia = mesh->GetTopology().GetVertexSurfaceElements(vnr);
        node.elements.ne = ia.Size();
        node.elements.ptr = (int*)&ia[0];
        
        FlatArray<SegmentIndex> bia = mesh->GetTopology().GetVertexSegments(vnr);
        node.bnd_elements.ne = bia.Size();
        node.bnd_elements.ptr = (int*)&bia[0];
        break;
      }
    case 1:
      {
        FlatArray<SegmentIndex> ia = mesh->GetTopology().GetVertexSegments(vnr);
        node.elements.ne = ia.Size();
        node.elements.ptr = (int*)&ia[0];
        
        FlatArray<int> bia = mesh->GetTopology().GetVertexPointElements(vnr);
        node.bnd_elements.ne = bia.Size();
        node.bnd_elements.ptr = (int*)&bia[0];
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
  node.vertices.ptr = mesh->GetTopology().GetEdgeVerticesPtr(nr);
  return node;
}

template <> NGX_INLINE DLL_HEADER const Ng_Node<2> Ngx_Mesh :: GetNode<2> (int nr) const
{
  Ng_Node<2> node;
  node.vertices.ptr = mesh->GetTopology().GetFaceVerticesPtr(nr);
  node.vertices.nv = (node.vertices.ptr[3] == 0) ? 3 : 4;
  node.surface_el = mesh->GetTopology().GetFace2SurfaceElement (nr+1)-1;
  return node;
}

