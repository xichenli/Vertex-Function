#ifndef CHI_NAMBU_H
#define CHI_NAMBU_H
#include "vertex.h"
#include "chi0_nambu.h"

///Class for the chi (susceptibility) vertex in particle hole-particle particle notation.
///Note that there are two 'chi'-s: the CLUSTER susceptibility and the coarse-grained lattice susceptibility. We have classes for both of them.
class chi_magnetic:public vertex_magnetic{
public:
  ///Constructor
  chi_magnetic(const alps::params &p, const k_space_structure &ks):vertex_magnetic(p, ks){}
//  void symmetrize_chi();

};

///This is a class for the CLUSTER susceptibility.
class chi_cluster_magnetic:public chi_magnetic{
public:
  ///Constructor
  chi_cluster_magnetic(const alps::params &p, const k_space_structure &ks):chi_magnetic(p, ks){}
  ///get the susceptibility out of the green's funciton and the Monte Carlo vertex
  void compute_chi(const mc_vertex &g4_phuu, const mc_vertex &g4_phud, alps::hdf5::archive &ar);
  void compute_chi(const mc_vertex &g4_m, alps::hdf5::archive &ar);
};

///This is a class for the coarse-grained LATTICE susceptibility. In the Jarrell literature it is usually called \f$ \overline{\chi} \f$
class chi_cglattice_magnetic:public chi_magnetic{
public:
  ///Constructor
  chi_cglattice_magnetic(const alps::params &p, const k_space_structure &ks):chi_magnetic(p, ks){}

  /// \todo compute the coarse grained LATTICE susceptibility from the cluster susceptibility
  //this should be done by solving chibar^-1 = chi_c^-1 - chi_c^0^-1 + chi0bar^-1
  void compute_chi_cglattice_directly(const chi0_cglattice_magnetic &chi0, const chi_cluster_magnetic &chic, const chi0_cluster_magnetic &chi0c);
};

#endif
