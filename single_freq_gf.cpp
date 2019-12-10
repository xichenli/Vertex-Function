#include"vertexaux.h"
#include"single_freq_gf.h"
#include<alps/hdf5.hpp>

#ifdef ALPSCore_HAS_MPI
#include<mpi.h>
#endif

void nambu_fermionic_green_function::read_density(alps::hdf5::archive&ar){
  std::cout<<"reading single particle density from density estimator"<<std::endl;
  density_[0]=0;//00
  density_[1]=0;//01
  density_[2]=0;//10
  density_[3]=0;//11
  std::vector<double> d_00_buff(ns_);
  std::vector<double> d_01_buff(ns_);
  std::vector<double> d_10_buff(ns_);
  std::vector<double> d_11_buff(ns_);
  double sign;
  ar>>alps::make_pvp("/simulation/results/Sign/mean/value",sign);
  ar>>alps::make_pvp("/simulation/results/densities_0_0_times_sign/mean/value", d_00_buff);
  ar>>alps::make_pvp("/simulation/results/densities_0_1_times_sign/mean/value", d_01_buff);
  ar>>alps::make_pvp("/simulation/results/densities_1_0_times_sign/mean/value", d_10_buff);
  ar>>alps::make_pvp("/simulation/results/densities_1_1_times_sign/mean/value", d_11_buff);
  for (int i=0; i< ns_; i++)
  {
    density_[0]+=d_00_buff[i]/ns_/sign; //this stores the average real space density
    density_[1]+=d_01_buff[i]/ns_/sign;
    density_[2]+=d_10_buff[i]/ns_/sign; //this stores the average real space density
    density_[3]+=d_11_buff[i]/ns_/sign;

  }
  set_hubbard_hifreq_tails(ar);
}

void nambu_fermionic_green_function::set_hubbard_hifreq_tails(alps::hdf5::archive&ar){
  std::cout<<"setting Hubbard high frequency tails for gf"<<std::endl;
  if(nf_!=2) throw std::logic_error("Hubbard with more than two spins?");
  
  double U;
  double doubleocc=0.;
  std::vector<double> doubleocc_buff(ns_);
  double sign;
  ar>>alps::make_pvp("/simulation/results/Sign/mean/value",sign); 
  ar>>alps::make_pvp("/parameters/U", U);
  ar>>alps::make_pvp("/simulation/results/doubleoc_times_sign/mean/value", doubleocc_buff);
  for (int i=0; i< ns_; i++)
  {
    doubleocc+=doubleocc_buff[i]/ns_/sign;
  }
  set_hubbard_hifreq_tails(U);
}
void nambu_fermionic_green_function::set_hubbard_hifreq_tails(double U){
  for (int i=0; i< ns_; i++){
    for (int f=0; f< nf_; f++){
      c1_[i*nf_*nf_+f*nf_+f]=1; //Only 00 and 11 have c1=1. For 01 and 10, c1_=0
      c2_[i*nf_*nf_+f*nf_+f]=(f==0?U*(density_[0]-0.5):(-U*(density_[3]-0.5)));//density_ are stored as density_[0],density_[1],density_[2],density_[3]
      //std::cout<<"setting c1_ to: "<<c1_[i*nf_+f]<<" and c2 to: "<<c2_[i*nf_+f]<<std::endl;
      //c3_[i*nf_+f]=U*U*(doubleocc-density_[spinup]*density_[spindn]);
    }
  }
}
 
void nambu_fermionic_green_function::read_single_freq(alps::hdf5::archive &ar)
{
  /**
   * IO function for single particle Green's function. receives the CT-AUX sim hdf file and reads data from <CODE>/simulation/results/G_omega_up_re"  <<k<<"/mean/value</CODE> and related locations. In addition reads density data out of <CODE>/simulation/results/density_up/mean/value</CODE>
   *
   */
  
  std::cout<<"reading sign from sim.h5 file"<<std::endl;
  double sign_result;
  ar>>alps::make_pvp("/simulation/results/Sign/mean/value", sign_result);
  std::cout<<"Using value of Sign= "<< sign_result<<std::endl;
  //In the old ALPS, single frequency GF is a signed observable, here we will need to divide the value in sim.h5 by sign
  std::cout<<"reading single particle GF from single-frequency estimator"<<std::endl;
  if(conv_!=beta_conv_single_freq) throw std::logic_error("this function reads single frequency GF");
  for(int k=0;k<ns_;++k){
    std::stringstream G_00_re_name, G_00_im_name;
    std::stringstream G_01_re_name, G_01_im_name;
    std::stringstream G_10_re_name, G_10_im_name;
    std::stringstream G_11_re_name, G_11_im_name;
    
    G_00_re_name<<"/simulation/results/G_omega_re"<<k<<"_0_0_times_sign/mean/value";
    G_00_im_name<<"/simulation/results/G_omega_im"<<k<<"_0_0_times_sign/mean/value";
    G_01_re_name<<"/simulation/results/G_omega_re"<<k<<"_0_1_times_sign/mean/value";
    G_01_im_name<<"/simulation/results/G_omega_im"<<k<<"_0_1_times_sign/mean/value";
    G_10_re_name<<"/simulation/results/G_omega_re"<<k<<"_1_0_times_sign/mean/value";
    G_10_im_name<<"/simulation/results/G_omega_im"<<k<<"_1_0_times_sign/mean/value";
    G_11_re_name<<"/simulation/results/G_omega_re"<<k<<"_1_1_times_sign/mean/value";
    G_11_im_name<<"/simulation/results/G_omega_im"<<k<<"_1_1_times_sign/mean/value";

    std::vector<double> G_00_re_mean; ar>>alps::make_pvp(G_00_re_name.str(), G_00_re_mean);
    std::vector<double> G_00_im_mean; ar>>alps::make_pvp(G_00_im_name.str(), G_00_im_mean);
    std::vector<double> G_01_re_mean; ar>>alps::make_pvp(G_01_re_name.str(), G_01_re_mean);
    std::vector<double> G_01_im_mean; ar>>alps::make_pvp(G_01_im_name.str(), G_01_im_mean);
    std::vector<double> G_10_re_mean; ar>>alps::make_pvp(G_10_re_name.str(), G_10_re_mean);
    std::vector<double> G_10_im_mean; ar>>alps::make_pvp(G_10_im_name.str(), G_10_im_mean);
    std::vector<double> G_11_re_mean; ar>>alps::make_pvp(G_11_re_name.str(), G_11_re_mean);
    std::vector<double> G_11_im_mean; ar>>alps::make_pvp(G_11_im_name.str(), G_11_im_mean);
    for(int i=0;i<nt_;++i){
      //std::cout<<i<<" "<<nt_<<" "<<k<<" "<<G_up_re_mean[i]<<" "<<G_up_im_mean[i]<<std::endl;
      operator()( i  ,k,0,0)=std::complex<double>(G_00_re_mean[i]  , G_00_im_mean[i]  )/sign_result;
      operator()( i  ,k,0,1)=std::complex<double>(G_01_re_mean[i]  , G_01_im_mean[i]  )/sign_result;
      operator()( i  ,k,1,0)=std::complex<double>(G_10_re_mean[i]  , G_10_im_mean[i]  )/sign_result;
      operator()( i  ,k,1,1)=std::complex<double>(G_11_re_mean[i]  , G_11_im_mean[i]  )/sign_result;
    }
  }
  read_density(ar);
}
//NAMBU: read single freq ferminonic matsubara function from file
void nambu_fermionic_matsubara_function::read_single_freq_ff(const std::string &in_file_name_gf){
  std::ifstream infile_gf(in_file_name_gf.c_str());
  std::cout<<"read single freq from file"<<in_file_name_gf<<std::endl;
  if(!infile_gf.is_open()) throw std::runtime_error("tried to read file that could not be opened: "+in_file_name_gf);
  double ignored;
  std::string ignored_line;
  //There are four lines of grid information
  char c;infile_gf>>c;
  if(c=='#'){
    for(int i=0;i<4;i++)
      std::getline(infile_gf,ignored_line);
  }
  for(int omega=0;omega<nt_;++omega){
    infile_gf>>ignored;
    for(int K=0;K<ns_;++K){
      double g00_in_real, g00_in_imag, g01_in_real, g01_in_imag, g10_in_real, g10_in_imag, g11_in_real, g11_in_imag;
      infile_gf>>g00_in_real>>g00_in_imag>>g01_in_real>>g01_in_imag>>g10_in_real>>g10_in_imag>>g11_in_real>>g11_in_imag;
      operator()(omega,K,0,0)=std::complex<double>(g00_in_real, g00_in_imag);
      operator()(omega,K,0,1)=std::complex<double>(g01_in_real, g01_in_imag);
      operator()(omega,K,1,0)=std::complex<double>(g10_in_real, g10_in_imag);
      operator()(omega,K,1,1)=std::complex<double>(g11_in_real, g11_in_imag);
//      std::cout<<" read in: "<<operator()(omega,K,0,0).real()<<" "<<operator()(omega,K,0,0).imag()<<" "<<operator()(omega,K,1,1).real()<<" "<<operator()(omega,K,1,1).imag()<<std::endl;
//      std::cout<<" read in: "<<f00_in_real<<" "<<f00_in_imag<<" "<<f01_in_real<<" "<<f01_in_imag<<" "<<f10_in_real<<" "<<f10_in_imag<<" "<<f11_in_real<<" "<<f11_in_imag<<std::endl;
    }
  }
}

void nambu_fermionic_green_function::read_two_freq(alps::hdf5::archive&ar){
  /**
   * IO function for single particle Green's function. receives the CT-AUX sim hdf file and reads data from <CODE>/simulation/results/G_omega_omega_up_re"  <<k<<"/mean/value</CODE> and related locations. This function does not read the density and it assumes that the Green's function is defined as a two-frequency GF.
   *
   */
  std::cout<<"reading sign from sim.h5 file"<<std::endl;
  double sign_result;
  ar>>alps::make_pvp("/simulation/results/Sign/mean/value", sign_result);
  std::cout<<"Using value of Sign= "<< sign_result<<std::endl;
  std::cout<<"reading single particle GF from two-frequency estimator"<<std::endl;
  
  if(conv_!=beta_conv_two_freq) throw std::logic_error("this function reads two frequency GF");
  for(int k=0;k<ns_;++k){
    std::stringstream G_00_re_name, G_00_im_name;
    std::stringstream G_01_re_name, G_01_im_name;
    std::stringstream G_10_re_name, G_10_im_name;
    std::stringstream G_11_re_name, G_11_im_name;
    G_00_re_name<<"/simulation/results/G_omega_omega_00_re"  <<k<<"_"<<k<<"_times_sign/mean/value";
    G_00_im_name<<"/simulation/results/G_omega_omega_00_im"  <<k<<"_"<<k<<"_times_sign/mean/value";
    G_01_re_name<<"/simulation/results/G_omega_omega_01_re"  <<k<<"_"<<k<<"_times_sign/mean/value";
    G_01_im_name<<"/simulation/results/G_omega_omega_01_im"  <<k<<"_"<<k<<"_times_sign/mean/value";
    G_10_re_name<<"/simulation/results/G_omega_omega_10_re"  <<k<<"_"<<k<<"_times_sign/mean/value";
    G_10_im_name<<"/simulation/results/G_omega_omega_10_im"  <<k<<"_"<<k<<"_times_sign/mean/value";
    G_11_re_name<<"/simulation/results/G_omega_omega_11_re"  <<k<<"_"<<k<<"_times_sign/mean/value";
    G_11_im_name<<"/simulation/results/G_omega_omega_11_im"  <<k<<"_"<<k<<"_times_sign/mean/value";
    
    std::vector<double> G_00_re_mean; ar>>alps::make_pvp(G_00_re_name.str(), G_00_re_mean);
    std::vector<double> G_00_im_mean; ar>>alps::make_pvp(G_00_im_name.str(), G_00_im_mean);
    std::vector<double> G_01_re_mean; ar>>alps::make_pvp(G_01_re_name.str(), G_01_re_mean);
    std::vector<double> G_01_im_mean; ar>>alps::make_pvp(G_01_im_name.str(), G_01_im_mean);
    std::vector<double> G_10_re_mean; ar>>alps::make_pvp(G_10_re_name.str(), G_10_re_mean);
    std::vector<double> G_10_im_mean; ar>>alps::make_pvp(G_10_im_name.str(), G_10_im_mean);
    std::vector<double> G_11_re_mean; ar>>alps::make_pvp(G_11_re_name.str(), G_11_re_mean);
    std::vector<double> G_11_im_mean; ar>>alps::make_pvp(G_11_im_name.str(), G_11_im_mean);
    
    for(int i=0;i<nt_;++i){
      int index=2*nt_*(i+nt_)+i+nt_;
      
      operator()(i,k,0,0)=std::complex<double>(G_00_re_mean[index], G_00_im_mean[index])/(beta_*ns_*sign_result);
      operator()(i,k,0,1)=std::complex<double>(G_01_re_mean[index], G_01_im_mean[index])/(beta_*ns_*sign_result);
      operator()(i,k,1,0)=std::complex<double>(G_10_re_mean[index], G_10_im_mean[index])/(beta_*ns_*sign_result);
      operator()(i,k,1,1)=std::complex<double>(G_11_re_mean[index], G_11_im_mean[index])/(beta_*ns_*sign_result);
    }
  }
  read_density(ar);
}


void nambu_fermionic_green_function::read_two_freq_err(alps::hdf5::archive&ar){
  /**
   * IO function for single particle Green's function. receives the CT-AUX sim hdf file and reads data from <CODE>/simulation/results/G_omega_omega_up_re"  <<k<<"/mean/value</CODE> and related locations. This function does not read the density and it assumes that the Green's function is defined as a two-frequency GF.
   *
   */
  std::cout<<"reading sign from sim.h5 file"<<std::endl;
  double sign_result;
  ar>>alps::make_pvp("/simulation/results/Sign/mean/value", sign_result);
  std::cout<<"Using value of Sign= "<< sign_result<<std::endl;
  std::cout<<"reading single particle GF from two-frequency estimator"<<std::endl;
  
  if(conv_!=beta_conv_two_freq) throw std::logic_error("this function reads two frequency GF");
  for(int k=0;k<ns_;++k){
    std::stringstream G_00_re_name, G_00_im_name;
    std::stringstream G_01_re_name, G_01_im_name;
    std::stringstream G_10_re_name, G_10_im_name;
    std::stringstream G_11_re_name, G_11_im_name;
    G_00_re_name<<"/simulation/results/G_omega_omega_00_re"  <<k<<"_"<<k<<"_times_sign/mean/err";
    G_00_im_name<<"/simulation/results/G_omega_omega_00_im"  <<k<<"_"<<k<<"_times_sign/mean/err";
    G_01_re_name<<"/simulation/results/G_omega_omega_01_re"  <<k<<"_"<<k<<"_times_sign/mean/err";
    G_01_im_name<<"/simulation/results/G_omega_omega_01_im"  <<k<<"_"<<k<<"_times_sign/mean/err";
    G_10_re_name<<"/simulation/results/G_omega_omega_10_re"  <<k<<"_"<<k<<"_times_sign/mean/err";
    G_10_im_name<<"/simulation/results/G_omega_omega_10_im"  <<k<<"_"<<k<<"_times_sign/mean/err";
    G_11_re_name<<"/simulation/results/G_omega_omega_11_re"  <<k<<"_"<<k<<"_times_sign/mean/err";
    G_11_im_name<<"/simulation/results/G_omega_omega_11_im"  <<k<<"_"<<k<<"_times_sign/mean/err";
    
    std::vector<double> G_00_re_mean; ar>>alps::make_pvp(G_00_re_name.str(), G_00_re_mean);
    std::vector<double> G_00_im_mean; ar>>alps::make_pvp(G_00_im_name.str(), G_00_im_mean);
    std::vector<double> G_01_re_mean; ar>>alps::make_pvp(G_01_re_name.str(), G_01_re_mean);
    std::vector<double> G_01_im_mean; ar>>alps::make_pvp(G_01_im_name.str(), G_01_im_mean);
    std::vector<double> G_10_re_mean; ar>>alps::make_pvp(G_10_re_name.str(), G_10_re_mean);
    std::vector<double> G_10_im_mean; ar>>alps::make_pvp(G_10_im_name.str(), G_10_im_mean);
    std::vector<double> G_11_re_mean; ar>>alps::make_pvp(G_11_re_name.str(), G_11_re_mean);
    std::vector<double> G_11_im_mean; ar>>alps::make_pvp(G_11_im_name.str(), G_11_im_mean);
    
    for(int i=0;i<nt_;++i){
      int index=2*nt_*(i+nt_)+i+nt_;
      
      operator()(i,k,0,0)=std::complex<double>(G_00_re_mean[index], G_00_im_mean[index])/(beta_*ns_*sign_result);
      operator()(i,k,0,1)=std::complex<double>(G_01_re_mean[index], G_01_im_mean[index])/(beta_*ns_*sign_result);
      operator()(i,k,1,0)=std::complex<double>(G_10_re_mean[index], G_10_im_mean[index])/(beta_*ns_*sign_result);
      operator()(i,k,1,1)=std::complex<double>(G_11_re_mean[index], G_11_im_mean[index])/(beta_*ns_*sign_result);
    }
  }
  read_density(ar);
}

///write single frequency green's function or self energy to file.
void nambu_fermionic_matsubara_function::write_single_freq(const std::string &out_file_name_gf) const{
  std::ofstream outfile_gf(out_file_name_gf.c_str());
  outfile_gf<<std::setprecision(15);
  for(int omega=-nt_;omega<nt_;++omega){
    for(int K=0;K<ns_;++K){
      outfile_gf<<(2*(omega)+1)*M_PI/beta_ << " "<<operator()(omega,K,0,0).real()<<" "<<operator()(omega,K,0,0).imag()<< " "<<operator()(omega,K,0,1).real()<<" "<<operator()(omega,K,0,1).imag()<< " "<<operator()(omega,K,1,0).real()<<" "<<operator()(omega,K,1,0).imag()<< " "<<operator()(omega,K,1,1).real()<<" "<<operator()(omega,K,1,1).imag()<<std::endl;
    }
  }
}

///human readable verion of write_single_freq
void nambu_fermionic_matsubara_function::write_single_freq_hr(const std::string &out_file_name_gf) const{
  std::ofstream outfile_gf(out_file_name_gf.c_str());
  outfile_gf<<std::setprecision(15);
  std::cout<<"writing "<<out_file_name_gf<<" # of omega="<<nt_<<" # of site="<<ns_;
  for(int omega=0;omega<nt_;++omega){
    outfile_gf<<(2*(omega)+1)*M_PI/beta_<<" ";
    for(int K=0;K<ns_;++K){
      outfile_gf<<operator()(omega,K,0,0).real()<<" "<<operator()(omega,K,0,0).imag()<< " "<<operator()(omega,K,0,1).real()<<" "<<operator()(omega,K,0,1).imag()<< " "<<operator()(omega,K,1,0).real()<<" "<<operator()(omega,K,1,0).imag()<< " "<<operator()(omega,K,1,1).real()<<" "<<operator()(omega,K,1,1).imag();
    }
    outfile_gf<<std::endl;
  }
}


 std::complex<double> nambu_fermionic_matsubara_function::hifreq( int n,  int site,  int f1,int f2) const{
 double wn=(2*n+1)*M_PI/beta_;
 std::complex<double> iwn(0., wn);
 
 return c1(site,f1,f2)/iwn+c2(site,f1,f2)/(iwn*iwn)+c3(site,f1,f2)/(iwn*iwn*iwn);
 }

