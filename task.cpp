#include<boost/program_options.hpp>
#include"task.h"
#include"gtest/gtest.h"

namespace po = boost::program_options;

///Use polymorphism to create different tasks. We switch between evaluating chis (e.g. for Alessandro) to computing second order sigmas and running GW, to computing phase boundaries etc.
task* task::task_factory(int argc, char **argv){
  std::string out_file_basename, sim_file_name, param_file_name, lattice_file_name,sigma_file_name, output_lattice_name="";
  int n_omega4=-1, n_omega4_bose=-1;
  int output_lattice_size=-1;
  task_type task_at_hand;
  factory_parse_command_line(argc, argv, out_file_basename, sim_file_name, param_file_name, lattice_file_name,  sigma_file_name, n_omega4, n_omega4_bose, output_lattice_size, output_lattice_name,  task_at_hand);
  const alps::params p=factory_read_parameters(sim_file_name, param_file_name, lattice_file_name, out_file_basename, sigma_file_name, n_omega4,  n_omega4_bose,output_lattice_size, output_lattice_name);
  if(task_at_hand==sigma_task_type){
    return new sigma_task(p);
  }else if(task_at_hand==nambu_lattice_susceptibility_task_type){
    return new nambu_lattice_susceptibility_task(p);
  }else if(task_at_hand==vertex_construction_task_type){
    return new vertex_construction_task(p);
  }else if(task_at_hand==BS_eigenvalues_task_type){
    return new BS_eigenvals_task(p);
  }else{
    throw std::runtime_error("task not understood.");
  }
}

void task::factory_parse_command_line(int argc, char ** argv, std::string &out_file_basename, std::string &sim_file_name, std::string &param_file_name, std::string &lattice_file_name, std::string &sigma_file_name, int &n_omega4, int &n_omega4_bose, int &output_lattice_size, std::string &output_lattice_name, task_type &task_at_hand) {
  namespace po = boost::program_options;
  /*! \brief Parses command line options
   *
   *  This function deals with the external parameter input, in particular finds the location of the Monte Carlo <CODE>sim.h5</CODE> file, finds the ALPS lattice file, and sets the prefix for the output files.
   */

  ///Parameter input
  std::string task_type_string;
  po::options_description desc( "Allowed options\n Please have a look at the doxygen documentation.");
  desc.add_options()("lattice",po::value < std::string> (&lattice_file_name)->default_value("cluster.xml"),"lattice file cluster.xml")
  ("basename",po::value < std::string> (&out_file_basename)->default_value("vert"),"base name of text output file")
  ("sim",po::value < std::string > (&sim_file_name)->default_value("sim.h5"), "hdf5 simulation file")
  ("param",po::value < std::string > (&param_file_name), "parameter file")
  ("help", "produce the help message")
  ("sigma",po::value < std::string > (&sigma_file_name), "self energy file (if required)")
  ("test", "run unit tests")
  ("n_omega", po::value<int>(&n_omega4), "modified number of fermionic frequencies for vertex")
  ("n_Omega", po::value<int>(&n_omega4_bose), "modified number of bosonic frequencies for vertex")
  ("lattice_size", po::value<int>(&output_lattice_size), "nmr task: Size of output for lattice quantities")
  ("lattice_name", po::value< std::string>(&output_lattice_name), "nmr task: Structure of lattice quantity for k_space_structure.")  
  ("task", po::value < std::string > (&task_type_string),"type of task to be performed: chi, sigma, doubleocc, lattice_susc, nmr_lattice_susc, neutron_lattice_susc, BS_eigenvals, or vertex_construction");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  if (vm.count("help")) {
    std::cout << desc;
    exit(1);
  }
  if (!vm.count("lattice"))
    throw std::runtime_error("you need to specify the cluster.xml file as lattice.");
  if (vm.count("test")){
    ::testing::InitGoogleTest(&argc, argv);
    exit(RUN_ALL_TESTS());
  }if(!vm.count("task")){ throw std::runtime_error("you need to specify a task."); }
  if (task_type_string=="sigma") task_at_hand=sigma_task_type;
  else if (task_type_string=="nambu_lattice_susc") task_at_hand=nambu_lattice_susceptibility_task_type;
  else if (task_type_string=="BS_eigenvals") task_at_hand=BS_eigenvalues_task_type;
  else if (task_type_string=="vertex_construction") task_at_hand=vertex_construction_task_type;
  else throw std::invalid_argument("not understood which task to do");
}
void task::define_parameters(alps::params &p){
  k_space_structure::define_parameters(p);
  p.define<double>("U", "interaction strength");
  p.define<int>("NOMEGA4_BOSE_ORIG", "number of bosonic frequencies in sim file");
  p.define<int>("NOMEGA4_ORIG", "number of fermionic frequencies in sim file");
}
alps::params task::factory_read_parameters(const std::string &sim_file_name, const std::string &param_file_name, const std::string &lattice_file_name, const std::string&out_file_basename, const std::string &sigma_file_name, int &n_omega4, int &n_omega4_bose,  int &output_lattice_size, std::string &output_lattice_name ){
  std::cout << "reading simulation output file" << std::endl;
  { std::ifstream param_file(param_file_name.c_str()); if(!param_file.good()) throw std::invalid_argument("please supply parameter file");}
  alps::params p(param_file_name);
  define_parameters(p);
  p["dca.LATTICE_LIBRARY"] = lattice_file_name;
  p["SIM_FILE_NAME"]=sim_file_name;
  p["OUT_FILE_NAME"]=out_file_basename;
  p["OUT_LATTICE_NAME"]=output_lattice_name;
  p["OUT_LATTICE_SIZE"]=output_lattice_size;
  if(sigma_file_name !="")
    p["SIGMA"]=sigma_file_name;
  if(n_omega4 !=-1){
    if(p.exists("ctaux.NOMEGA4")){ int n_omega4_orig=p["ctaux.NOMEGA4"]; p["NOMEGA4"]=n_omega4_orig;}
    int n_omega4_orig=p["NOMEGA4"];
    if(n_omega4_orig < n_omega4) throw std::runtime_error("can't have more fermionic frequencies than measured!");
    p["NOMEGA4_ORIG"]=n_omega4_orig;
    p["NOMEGA4"]=n_omega4;
  }else{
    if(p.exists("ctaux.NOMEGA4")){ int n_omega4_orig=p["ctaux.NOMEGA4"]; p["NOMEGA4"]=n_omega4_orig;}
    p["NOMEGA4_ORIG"]=p["NOMEGA4"];
  }
  if(n_omega4_bose !=-1){
    if(p.exists("ctaux.NOMEGA4_BOSE")){ int n_omega4_orig=p["ctaux.NOMEGA4_BOSE"]; p["NOMEGA4_BOSE"]=n_omega4_orig;}
    int n_omega4_bose_orig=p["NOMEGA4_BOSE"];
    p["NOMEGA4_BOSE_ORIG"]=n_omega4_bose_orig;
    if(n_omega4_bose_orig < n_omega4_bose) throw std::runtime_error("can't have more bosonic frequencies than measured!");
    p["NOMEGA4_BOSE"]=n_omega4_bose;
  }else{
    if(p.exists("ctaux.NOMEGA4_BOSE")){ int n_omega4_orig=p["ctaux.NOMEGA4_BOSE"]; p["NOMEGA4_BOSE"]=n_omega4_orig;}
    p["NOMEGA4_BOSE_ORIG"]=p["NOMEGA4_BOSE"];
  }
  std::cout << "done reading parameters" << std::endl;
  std::cout<<"p is: "<<p<<std::endl;
  std::cout<<"p nomega4 is: "<<p["NOMEGA4"]<<std::endl;
  return p;
}
///general constructor for task, this will initialize everything that is common to all the different tasks we need to perform
task::task(const alps::params &p):p_(p),ks_(p_),ar_(p["SIM_FILE_NAME"], alps::hdf5::archive::READ)
{
  std::cout<<"task construction."<<std::endl;
  initialize_common_parameters();

  out_file_basename_=p_["OUT_FILE_NAME"].as<std::string>();
  std::cout<<"done task construction."<<std::endl;
}
void task::initialize_common_parameters(){
 U_ = p_["U"];
 beta_ = p_["BETA"];
 n_sites_ = p_["dca.SITES"];
 n_freq_ = p_["NMATSUBARA"];
 w4_range_ = (int)p_["ctaux.NOMEGA4"]+(int)(p_["ctaux.NOMEGA4_BOSE"]);
 std::cout<<"w4_range="<<(int)p_["ctaux.NOMEGA4"]<<"+"<<(int)p_["ctaux.NOMEGA4_BOSE"]<<w4_range_<<std::endl;
 mu_ = p_["MU"];
}

//*********************  NAMBU ****************************//
nambu_fermionic_green_function task::nambu_read_single_freq_gf(){
  std::cout << "reading single-particle g2." << std::endl;
  std::cout<<"nfreq: "<<n_freq_<<" nsites: "<<n_sites_<<" beta: "<<beta_<<std::endl;
  nambu_fermionic_green_function sfg2(n_freq_, n_sites_, 2, beta_, beta_conv_single_freq);
  sfg2.read_single_freq_ff("G_omega_0");
  sfg2.write_single_freq(out_file_basename_ + "_sfgf");
  sfg2.write_single_freq_hr(out_file_basename_ + "_sfgfhr");
  //sfg2.write_single_freq_hr(out_file_basename_ + "_sym_sfgfhr");
  return sfg2;
}
nambu_fermionic_green_function task::nambu_read_two_freq_gf(){
  std::cout << "reading two-frequency g2." << std::endl;
  nambu_fermionic_green_function g2(w4_range_, n_sites_, 2, beta_, beta_conv_two_freq);
  std::cout<<"nfreq: "<<w4_range_<<" nsites: "<<n_sites_<<" beta: "<<beta_<<std::endl;
  g2.read_two_freq(ar_);
  // g2.write_single_freq_vienna(vienna_to_index, out_file_basename + "_gf");
  g2.write_single_freq( out_file_basename_ + "_gf");
  g2.write_single_freq_hr(out_file_basename_ + "_gfhr");
  //g2.symmetrize(ks_);
  //g2.write_single_freq_hr(out_file_basename_ + "_sym_gfhr");
  return g2;
}
nambu_fermionic_self_energy task::nambu_import_self_energy(){
  nambu_fermionic_self_energy sigma(w4_range_, n_sites_, 2, beta_);
  sigma.read_single_freq_ff(p_["SIGMA"]);
  return sigma;
}
