#ifndef VALIDATE_HPP
#define VALIDATE_HPP

namespace netgen
{
  
  void GetPureBadness(Mesh & mesh, NgArray<double> & pure_badness,
		      const NgBitArray & isnewpoint);
  double Validate(const Mesh & mesh, NgArray<ElementIndex> & bad_elements,
		  const NgArray<double> & pure_badness, 
		  double max_worsening, const bool uselocalworsening,
		  NgArray<double> * quality_loss = NULL);
  void RepairBisection(Mesh & mesh, NgArray<ElementIndex> & bad_elements, 
		       const NgBitArray & isnewpoint, const Refinement & refinement,
		       const NgArray<double> & pure_badness, 
		       double max_worsening, const bool uselocalworsening,
		       const NgArray< NgArray<int,PointIndex::BASE>* > & idmaps);

}

#endif // VALIDATE_HPP
