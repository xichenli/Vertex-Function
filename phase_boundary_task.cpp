#include"task.h"
#include "chi_nambu.h"
#include "chi0_nambu.h"
#include "F_nambu.h"
#include"alps/numeric/vector_functions.hpp"
#include <sstream>
#include <string>

using namespace alps::numeric;
using namespace std;

void nambu_lattice_susceptibility_task::work(){
  mc_vertex g4(p_,ks_);
  chi0_cluster_magnetic  chi0m_vertex(p_, ks_);
  chi_cluster_magnetic   chim_vertex(p_, ks_);

  nambu_fermionic_green_function sfg2=nambu_read_single_freq_gf();
  g4.read(ar_,"m","value");
  nambu_fermionic_self_energy sigma=nambu_import_self_energy();

  chi0m_vertex.compute_chi0(sfg2);
  chim_vertex.compute_chi(g4, ar_);
  F_magnetic Fm_vertex(p_, ks_);
  Fm_vertex.invert_F_magnetic(chim_vertex,chi0m_vertex);
  Fm_vertex.equation_of_motion(sfg2);
}


void sigma_task::work(){std::cout<<"to be implemented"<<std::endl;}
void vertex_construction_task::work(){std::cout<<"to be implemented"<<std::endl;}
void BS_eigenvals_task::work(){std::cout<<"to be implemented"<<std::endl;}
