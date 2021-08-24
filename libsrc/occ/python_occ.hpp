class DirectionalInterval
{
public:
  gp_Vec dir;
  double minval = -1e99;
  double maxval = 1e99;
  bool openmin = false, openmax = false;
  
  DirectionalInterval (gp_Vec adir) : dir(adir) { ; }
  DirectionalInterval (const DirectionalInterval & i2)
    : dir(i2.dir), minval(i2.minval), maxval(i2.maxval) { ; }
  
  DirectionalInterval operator< (double val) const
  {
    cout << "create interval with < " << ", val = " << val << endl;
    cout << "me = " << this->minval << " - " << this->maxval << endl;    
    DirectionalInterval i2 = *this;
    i2.maxval = val;

    cout << "resulting i = " << i2.minval << " - " << i2.maxval << endl;
    return i2;
  }

  DirectionalInterval operator> (double val) const
  {
    cout << "create interval with > " << ", val = " << val << endl;
    cout << "me = " << this->minval << " - " << this->maxval << endl;
    
    DirectionalInterval i2 = *this;
    i2.minval = val;
    cout << "resulting i = " << i2.minval << " - " << i2.maxval << endl;
    
    return i2;
  }

  
  bool Contains (gp_Pnt p, double eps = 1e-8)
  {
    // cout << "Contains point " << p.X() << "," << p.Y() << "," << p.Z() << " ? " << endl;
    double val = dir.X()*p.X() + dir.Y()*p.Y() + dir.Z() * p.Z();
    // cout << "minval = " << minval << ", val = " << val << " maxval = " << maxval << endl;
    if (openmin) {
      if (val < minval+eps) return false;
    } else {
      if (val < minval-eps) return false;
    }
    if (openmax) {
      if (val > maxval-eps) return false;
    } else {
      if (val > maxval+eps) return false;
    }
    return true;
  }
};


inline gp_Pnt Center (TopoDS_Shape shape)
{
  GProp_GProps props;
  switch (shape.ShapeType())
    {
    case TopAbs_FACE:
      BRepGProp::SurfaceProperties (shape, props); break;
    default:
      BRepGProp::LinearProperties(shape, props);
    }
  return props.CentreOfMass();
}

