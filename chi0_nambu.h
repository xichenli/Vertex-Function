#ifndef CHI0_NAMBU_H
#define CHI0_NAMBU_H
#include"vertex.h"

///Class for the  \f$ \chi_0\f$ (bare susceptibility) vertex in nambu formalism
class chi0_magnetic:public vertex_magnetic{
public:
  ///Constructor
  chi0_magnetic(const alps::params &p, const k_space_structure &ks):vertex_magnetic(p,ks){}
 };

///Class for the bare cluster \f$ \chi_0\f$ (bare susceptibility) vertex in nambu formalism
class chi0_cluster_magnetic:public chi0_magnetic{
public:
  ///Constructor
  chi0_cluster_magnetic(const alps::params &p, const k_space_structure &ks):chi0_magnetic(p,ks){}
  
  ///get the chi0 out of the Green's function
  void compute_chi0(const nambu_fermionic_green_function &g2);

};

///Class for the bare coarse-grained lattice \f$ \overline\chi_0\f$ (bare susceptibility) vertex in particle hole-particle particle notation
class chi0_cglattice_magnetic:public chi0_magnetic{
public:
  ///Constructor
  chi0_cglattice_magnetic(const alps::params &p, const k_space_structure &ks);
  ~chi0_cglattice_magnetic();
  ///get the chi0 out of the Green's function
  void compute_chi0(const nambu_fermionic_self_energy &Sigma);
private:
  void precompute_integration_weights();
  std::vector<std::vector<double> >epsilon_kckl_;
  std::vector<std::vector<double> >epsilon_kckl_plus_pipi_;
  std::vector<double> weight_kl_;
  std::vector<std::vector<double> >symmetry_kckl_;//For example, the symmetry factor for dwave is symmetry_kckl_=cos(kx)-cos(ky)
  double mu_;
};

#endif
