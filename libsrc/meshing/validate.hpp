#ifndef VALIDATE_HPP
#define VALIDATE_HPP

namespace netgen
{
  
  void GetPureBadness(Mesh & mesh, NgArray<double> & pure_badness,
		      const TBitArray<PointIndex> & isnewpoint);
  double Validate(const Mesh & mesh, NgArray<ElementIndex> & bad_elements,
		  const NgArray<double> & pure_badness, 
		  double max_worsening, const bool uselocalworsening,
		  NgArray<double> * quality_loss = NULL);
  void RepairBisection(Mesh & mesh, NgArray<ElementIndex> & bad_elements, 
		       const TBitArray<PointIndex> & isnewpoint, const Refinement & refinement,
		       const NgArray<double> & pure_badness, 
		       double max_worsening, const bool uselocalworsening,
		       const NgArray< idmap_type* > & idmaps);

}

#endif // VALIDATE_HPP
