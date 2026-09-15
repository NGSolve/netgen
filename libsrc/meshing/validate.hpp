#ifndef VALIDATE_HPP
#define VALIDATE_HPP

namespace netgen
{
  
  void GetPureBadness(Mesh & mesh, Array<double, PointIndex> & pure_badness,
                      const TBitArray<PointIndex> & isnewpoint);
  double Validate(const Mesh & mesh, Array<ElementIndex> & bad_elements,
                  const Array<double, PointIndex> & pure_badness, 
                  double max_worsening, const bool uselocalworsening,
                  Array<double, ElementIndex> * quality_loss = NULL);
  void RepairBisection(Mesh & mesh, Array<ElementIndex> & bad_elements, 
                       const TBitArray<PointIndex> & isnewpoint, const Refinement & refinement,
                       const Array<double, PointIndex> & pure_badness, 
                       double max_worsening, const bool uselocalworsening,
                       const Array< idmap_type* > & idmaps);

}

#endif // VALIDATE_HPP
