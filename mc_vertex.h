#ifndef MCVERTEX
#define MCVERTEX
#include"vertexaux.h"
#include"alps/hdf5/archive.hpp"
#include"alps/hdf5/vector.hpp"
/// \brief A class for describing vertices read from Monte Carlo.
///
/// This class knows how to get vertex data from Monte Carlo. It is mainly used by the vertex_phpp class to extract the \f$\chi\f$.
class mc_vertex{
public:
  mc_vertex(const alps::params &p, const k_space_structure &ks);
  static void set_parameters(alps::params &p);
  void read(alps::hdf5::archive &ar,std::string type,std::string e_or_v);
  int freqindex3(int nu, int omega1, int omega2) const{
    return nu*n_omega4_*n_omega4_*4+(omega1+n_omega4_)*n_omega4_*2+(omega2+n_omega4_);
  }

  vertex_t vertex00;
  vertex_t vertex01;

private:
  int n_omega4_;
  int n_omega4_bose_;
  int n_sites_;
  double beta_;
  const k_space_structure *ks_;
};

#endif
