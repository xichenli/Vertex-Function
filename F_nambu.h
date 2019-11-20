#ifndef F_NAMBU_H
#define F_NAMBU_H
#include"vertex.h"
#include "chi0_nambu.h"
#include "chi_nambu.h"
class F_magnetic:public vertex_magnetic{
public:
  F_magnetic(const alps::params &p, const k_space_structure &ks):vertex_magnetic(p, ks){}
  //void invert_F_nambu(const Gamma_dmst &Gamma, const vertex_dmst &chi_times_chi0inv_dmst);
  void invert_F_magnetic(const chi_magnetic &chi, const chi0_magnetic &chi0);
 // void check_crossing_symmetry();
  void equation_of_motion(const nambu_fermionic_green_function &g2);
};
#endif
