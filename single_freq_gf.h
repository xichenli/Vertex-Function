#ifndef SINGLE_FREQ_GF_H
#define SINGLE_FREQ_GF_H
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <complex>
#include <vector>
#include <alps/hdf5/archive.hpp>
#include <alps/params.hpp>
#include "k_space_structure.h"


/// \brief Enum type to designate if Green's function is stored as single frequency or two frequencies.
///
/// The same single particle Green's functions can be stored in two different formats... which one is it?
/// Vertex functions should check this enum to make sure the right type of GF is being used.
enum beta_convention{
  ///Frequency convention of single frequency GF, which is 1/energy
  beta_conv_single_freq,
  ///Frequency convention of two frequency GF, which is 1/energy^2
  beta_conv_two_freq,
  ///Frequency convention not determined (yet)
  beta_conv_not_set
};

///enum for spin up and spin down
enum spinflavors {spinup=0, spindn=1} ;

template <typename T> class nambu_single_freq_gf{
public:
  //construction and destruction, assignement and copy constructor
  ///constructor: how many time slices, how many sites, how many flavors
  nambu_single_freq_gf( int ntime,  int nsite,  int nflavor, double beta):nt_(ntime), ns_(nsite), nf_(nflavor), ntns_(ntime*nsite), beta_(beta){
    val_.resize(nt_*ns_*nf_*nf_);//Now instead of an up and a down, there are 00, 01, 10 ,11
    density_.resize(nf_*nf_, 0.);
    c0_.resize(ns_*nf_*nf_,0.);
    c1_.resize(ns_*nf_*nf_,0.);
    c2_.resize(ns_*nf_*nf_,0.);
    c3_.resize(ns_*nf_*nf_,0.);
  }
  void clear(){ memset(&(val_[0]), 0, ns_*nt_*nf_*nf_*sizeof(T)); }
  //access of vectors and elements
  
  ///access element with given time, site, and flavor
  inline T &operator()( int t,  int site,  int f1,int f2){
    if(t<0 || t >= nt_) throw std::logic_error("Here 3 GF access out of bounds!");
    return val_[t+nt_*site+ntns_*f1+ntns_*nf_*f2];
  }
  ///access element with given time, site, and flavor (const reference)
  inline const T operator()( int t,  int site,  int f1,int f2)const{
    if(t<0 || t >= nt_) throw std::logic_error("const GF access out of bounds!");
    return val_[t+nt_*site+ntns_*f1+ntns_*nf_*f2];
  }
  
  //get the highfreq coefficients
  ///high frequency coefficient of the \f$ \frac{1}{i \omega_n} \f$ term
  double &c1(int s1, int f1,int f2){ return c1_[s1*nf_*nf_+nf_*f1+f2]; }
  ///high frequency coefficient of the \f$ \frac{1}{i \omega_n} \f$ term
  const double &c1(int s1,  int f1,int f2) const{ return c1_[s1*nf_*nf_+nf_*f1+f2]; }
  ///high frequency coefficient of the \f$ \frac{1}{(i \omega_n)^2} \f$ term
  double &c2(int s1, int f1,int f2){ return c2_[s1*nf_*nf_+nf_*f1+f2]; }
  ///high frequency coefficient of the \f$ \frac{1}{(i \omega_n)^2} \f$ term
  const double &c2(int s1, int f1,int f2) const{ return c2_[s1*nf_*nf_+nf_*f1+f2]; }
  ///high frequency coefficient of the \f$ \frac{1}{(i \omega_n)^3} \f$ term
  double &c3(int s1, int f1,int f2){ return c3_[s1*nf_*nf_+nf_*f1+f2]; }
  ///high frequency coefficient of the \f$ \frac{1}{(i \omega_n)^3} \f$ term
  const double &c3(int s1, int f1,int f2) const{ return c3_[s1*nf_*nf_+nf_*f1+f2]; }
  ///constant high frequency coefficient (only needed for self energies etc)
  const double &c0(int s1, int f1,int f2) const{ return c0_[s1*nf_*nf_+nf_*f1+f2]; }
  std::vector<double> &c1(){ return c1_; }
  ///high frequency coefficient of the \f$ \frac{1}{(i \omega_n)^2} \f$ term
  std::vector<double> &c2(){ return c2_; }
  ///high frequency coefficient of the \f$ \frac{1}{(i \omega_n)^3} \f$ term
  std::vector<double> &c3(){ return c3_; }
  const std::vector<double> &c1()const { return c1_; }
  const std::vector<double> &c2()const { return c2_; }
  const std::vector<double> &c3()const { return c3_; }
  //size information
  ///how many flavors do we have? (flavors are usually spins, GF of different flavors are zero)
  inline const  int &nflavor()const{return nf_;}
  ///return # of sites
  inline const  int &nsite()const{return ns_;}
  ///return # of imaginary time values
  inline const  int &ntime()const{return nt_;}
  ///return # of matsubara frequencies. Exactly equivalent to ntime().
  ///In case of a Matsubara GF 'ntime' sounds odd -> define 'nfreq' instead.
  inline const  int &nfreq()const{return nt_;} //nfreq is an alias to ntime - more intuitive use for Matsubara GF
  ///get out the inverse temperature
  const double &beta() const{return beta_;}
  ///compute maximum difference to another GF
  double max_abs_diff(const nambu_single_freq_gf &other)const{double d=0.; for(int i=0;i<val_.size();++i){double m=std::abs(val_[i]-other.val_[i]); d=d>m?d:m;} return d;}
  
protected:  //const values
  const  int nt_; ///imag time points
  const  int ns_; ///number of sites
  const  int nf_; ///number of flavors
  const  int ntns_; ///nt*ns
  // the actual values and errors.
  std::vector<T> val_;
  std::vector<double> c0_, c1_, c2_, c3_;
  std::vector<double> density_;
  const double beta_;
  
};
/// \brief base class for green's functions and self-energies, combining the common concepts for the two.
///
/// This class stores a fermionic Green's function.
class nambu_fermionic_matsubara_function:public nambu_single_freq_gf<std::complex<double> >{
  
public:
  ///Constructor
  nambu_fermionic_matsubara_function(int ntime,  int nsite,  int nflavor, double beta):nambu_single_freq_gf<std::complex<double> >(ntime,nsite,nflavor,beta){}
  ///element access function taking care of antisymmetry for negative elements: This function is fermionic, the imag part flips under negation of the argument, the real part stays the same: G(-Q)=(Re G(Q), -Im G(Q))
  inline std::complex<double> operator()( int t,  int site,  int f1,int f2)const{
    if(t<-nt_ || t >= nt_) throw std::logic_error("Here 2 matsubara const GF access out of bounds!");
    if(t<0){
      return std::complex<double>(val_[-t-1+nt_*site+ntns_*f1+ntns_*nf_*f2].real(), -val_[-t-1+nt_*site+ntns_*f1+ntns_*nf_*f2].imag());
    }else{
      return val_[t+nt_*site+ntns_*f1+ntns_*nf_*f2];
    }
  }
  inline std::complex<double> &operator()( int t,  int site,  int f1,int f2){
    //  if(t<0 || t >= nt_) throw std::logic_error("Here 1 matsubara GF access out of bounds!");
    return nambu_single_freq_gf<std::complex<double> >::operator()(t,site,f1,f2);
  }
  
  ///file output function
  void write_single_freq(const std::string &out_file_name_gf) const;
  ///file output function, better suited for plotting
  void write_single_freq_hr(const std::string &out_file_name_gf) const;
  ///file input function
  void read_single_freq_ff(const std::string &in_file_name_gf);
  /*
  ///symmetrize the matsubara function
  void symmetrize(const k_space_structure &ks);
  ///spin symmetrization only
  void symmetrize_spin();
  ///k-space symmetrization only
  void symmetrize_k_space(const k_space_structure &ks);
  
   */
  ///evaluate the high frequency behavior
  std::complex<double> hifreq( int t,  int site,  int f1,int f2) const;
};

/// \brief fermionic Green's function, as a function of a fermionic frequency.
///
/// This class stores a fermionic Green's function.
class nambu_fermionic_green_function:public nambu_fermionic_matsubara_function{
public:
  //construction and destruction, assignement and copy constructor
  ///constructor: how many time slices, how many sites, how many flavors
  nambu_fermionic_green_function(int ntime,  int nsite,  int nflavor, double beta, beta_convention conv):
  nambu_fermionic_matsubara_function(ntime,nsite,nflavor,beta),
  conv_(conv){}
  
  ///file input function
  void read_single_freq(alps::hdf5::archive &ar);
  ///file input function
  void read_two_freq(alps::hdf5::archive &ar);
  void read_two_freq_err(alps::hdf5::archive &ar);
  
  ///read the density from the ALPS density estimator for the total density
  void read_density(alps::hdf5::archive &ar);
  
  ///set Hubbard high frequency tails
  void set_hubbard_hifreq_tails(alps::hdf5::archive&ar);
  ///set Hubbard high frequency tails from known U
  void set_hubbard_hifreq_tails(double U);
  ///const access to density up
  const std::vector<double> &density() const{ return density_;}
  ///const access to total density
  double total_density() const{ return std::accumulate(density_.begin(), density_.end(), 0.);}
  ///compute density
  ///is this a single or a double frequency Green's function, i.e. dimension 1/Energy or dimension 1/Energy**2?
  beta_convention conv() const{return conv_;}
  /*
  ///this is a test/debug function to set the GF to a local (i.e. k-independent) Bethe lattice gf.
  void initialize_to_local_bethe(double mu);
  */
private:
  const beta_convention conv_;
};

class nambu_fermionic_self_energy:public nambu_fermionic_matsubara_function{
public:
  
  nambu_fermionic_self_energy(int ntime, int nsite, int nflavor, double beta):nambu_fermionic_matsubara_function(ntime,nsite,nflavor,beta){}
};

#endif
