#ifndef TASK_HPP
#define TASK_HPP

#include"vertexaux.h"
#include"mc_vertex.h"
#include"chi_nambu.h"
#include"chi0_nambu.h"


#include"alps/hdf5/archive.hpp"

enum task_type{
  sigma_task_type,
  vertex_construction_task_type,
  BS_eigenvalues_task_type,
  nambu_lattice_susceptibility_task_type,
};

///This program has a couple of purposes: create lambdas, create chis for alessandro, compute phase boundaries, etc. This selector picks out the tasks that should be performed via polymorphism.
class task{
public:
  static task* task_factory(int argc, char **argv);
  task(const alps::params &p);
  virtual ~task(){}
  virtual void work()=0;
private:
  static void factory_parse_command_line(int argc, char ** argv, std::string &out_file_basename, std::string &sim_file_name, std::string &param_file_name, std::string &lattice_file_name, std::string &sigma_file_name, int &n_omega4, int &n_omega4_bose, int &lattice_size, std::string &lattice_name, task_type &task_at_hand);
  static alps::params factory_read_parameters(const std::string &sim_file_name, const std::string &param_file_name, const std::string &lattice_file_name, const std::string&out_file_basename, const std::string &sigma_file_name, int &n_omega4, int &n_omega4_bose,  int &output_lattice_size, std::string &output_lattice_name );
  static void define_parameters(alps::params &p);
  void initialize_common_parameters();

protected:

  ///reading a single frequency nambu green's funciton from file
  nambu_fermionic_green_function nambu_read_single_freq_gf();
  nambu_fermionic_green_function nambu_read_two_freq_gf();
  nambu_fermionic_self_energy nambu_import_self_energy();

  ///Command line argument
  int argc_;
  char *argv_;
  ///ALPS parameter file that contains all the important parameters
  alps::params p_;
  ///k-space information: adding and subtracting momenta, k-points, etc.
  k_space_structure ks_;
  ///base name for output files generated throughout the code
  std::string out_file_basename_;
  ///hdf5 input file
  alps::hdf5::archive ar_;

  ///interaction strength
  double U_;
  ///inverse temperature
  double beta_;
  ///number of cluster sites
  int n_sites_;
  ///number of matsubara frequencies
  int n_freq_;
  ///range of frequencies for the single particle quantities and the measurements.
  int w4_range_;
  ///chemical potential
  double mu_;
  ///translation vectors to the vienna format and back to our format
  std::vector<int> index_to_vienna_, vienna_to_index_;

};

class sigma_task: public task{
public:
  sigma_task(const alps::params &p):task(p){}
  virtual ~sigma_task(){};
  void work();
};

class nambu_lattice_susceptibility_task: public task{
public:
  nambu_lattice_susceptibility_task(const alps::params &p):task(p){}
  virtual ~nambu_lattice_susceptibility_task(){};
  void work();
};

class vertex_construction_task: public task{
public:
  vertex_construction_task(const alps::params &p):task(p){}
  virtual ~vertex_construction_task(){};
  void work();
};
class BS_eigenvals_task: public task{
public:
  BS_eigenvals_task(const alps::params &p):task(p){}
  virtual ~BS_eigenvals_task(){};
  void work();
};

#endif //TASK_HPP
