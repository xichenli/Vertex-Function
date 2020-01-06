#ifndef VERTEX_H
#define VERTEX_H
#include"vertexaux.h"
#include"single_freq_gf.h"
#include"k_space_structure.h"

/// \brief The main vertex class.
///
/// This class stores vertices. The two main formats, ph/pp and dm/st are derived from this class.
class vertex{
public:
  vertex(const alps::params &p,const k_space_structure &ks);
  virtual ~vertex(){}
  ///compute maximum absolute difference between v and *this.
  double max_diff(const vertex &v) const;
  double max_rel_diff(const vertex &v) const;
  ///number of bosonic frequencies
  int n_omega4_bose() const{return n_omega4_bose_;}
  ///number of fermionic frequencies
  int n_omega4()const{return n_omega4_;}
  ///number of sites
  int n_sites()const{return n_sites_;}

  ///matrix EV routine. There ahs to be a better place to put those...
  static void inplace_eigenvalues_eigenvectors(matrix_t &A, matrix_t &U, vector_t &EV);

  //matrix output/debug routine. Only used here so far.
  static std::ostream& matlab_print(std::ostream &os, const matrix_t &A);

  static void triple_invert_inv_minv_pinv(matrix_t &res, const matrix_t &A, const matrix_t &B, const matrix_t &C);
  static void invert_oneplussmall_chic(matrix_t &res, const matrix_t &A, const matrix_t &B, const matrix_t &C);
  static void solve_a_eq_b_X_c(const matrix_t &A, const matrix_t &B, matrix_t &X, const matrix_t &C) ;
  //equation solver routine for the equation \f$ A = B X C\f$ using inversion instead of solver
  static void solve_a_eq_b_X_c_directly(const matrix_t &A, const matrix_t &B, matrix_t &X, const matrix_t &C) ;
  virtual void broadcast_data();
  matrix_t combine_matrix(const matrix_t &x00, const matrix_t &x01);
  void separate_matrix(const matrix_t &matrix_whole, matrix_t &x00, matrix_t &x01);


  std::complex<double> &v00(int K, int Kprime, int Q, int omega, int omegaprime, int nu){
    return a_[bindex(Q,nu)](findex(K,omega),findex(Kprime,omegaprime));
  }
  ///access function for individual element, anomalous
  std::complex<double> &v01(int K, int Kprime, int Q, int omega, int omegaprime, int nu){
    return b_[bindex(Q,nu)](findex(K,omega),findex(Kprime,omegaprime));
  }
  ///const access function for individual element, normal
  const std::complex<double> &v00(int K, int Kprime, int Q, int omega, int omegaprime, int nu)const{
    return a_[bindex(Q,nu)](findex(K,omega),findex(Kprime,omegaprime));
  }
  ///const access function for individual element, anomalous
  const std::complex<double> &v01(int K, int Kprime, int Q, int omega, int omegaprime, int nu)const{
    return b_[bindex(Q,nu)](findex(K,omega),findex(Kprime,omegaprime));
  }
  ///const access function for matrix for same frequency and momentum transfer, normal
  const matrix_t &v00(int Q, int nu)const{
    return a_[bindex(Q,nu)];
  }
  ///const access function for matrix for same frequency and momentum transfer, anomalous
  const matrix_t &v01(int Q, int nu)const{
    return b_[bindex(Q,nu)];
  }
  /// access function for matrix for same frequency and momentum transfer, normal
  matrix_t &v00(int Q, int nu){
    return a_[bindex(Q,nu)];
  }
  /// access function for matrix for same frequency and momentum transfer, anomalous
  matrix_t &v01(int Q, int nu){
    return b_[bindex(Q,nu)];
  }

protected:
  int bindex(int Q, int nu)const{ return Q*(2*n_omega4_bose_+1)+(nu+n_omega4_bose_);} //bosonic multi-index
  int findex(int K, int omega)const{ return K*(2*n_omega4_)+(omega+n_omega4_);} //fermionic multi-index

  vertex_matrix_t a_;
  vertex_matrix_t b_;

  int n_omega4_;
  int n_omega4_bose_;
  int n_sites_;
  double beta_;
  const k_space_structure *ks_;
  double U_;
};
#endif
#ifndef VERTEX_NAMBU_H
#define VERTEX_NAMBU_H
#include"vertex.h"
/// \brief ph_uu/ph_ud/pp_uu/pp_ud notation vertex class.
///
/// This class stores vertices in the dm/st format. Most of the vertex calculations are performed in this format.
class vertex_magnetic:public vertex{
public:
  vertex_magnetic(const alps::params &p, const k_space_structure &ks):vertex(p, ks){}

  ///access function for individual element, normal
  std::vector<std::vector<std::complex<double> > >compute_static_magnetic_susceptibility(int n_omega) const;
};
#endif
