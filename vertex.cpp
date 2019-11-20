#include"vertex.h"

#include"blasheader.h"

#ifdef ALPSCore_HAS_MPI
#include<mpi.h>
#endif

vertex::vertex(const alps::params &p, const k_space_structure &ks):ks_(&ks){
  n_omega4_=p["NOMEGA4"];
  n_omega4_bose_=p["NOMEGA4_BOSE"];
  n_sites_=p["dca.SITES"];
  beta_=p["BETA"];
  U_ = p["U"];

  for(int i=0;i<n_sites_*(2*n_omega4_bose_+1);++i){
    a_.push_back(Eigen::MatrixXcd::Zero(n_omega4_*n_sites_*2,n_omega4_*n_sites_*2));
    b_.push_back(Eigen::MatrixXcd::Zero(n_omega4_*n_sites_*2,n_omega4_*n_sites_*2));
  }
}

double vertex::max_diff(const vertex &v) const{
  double max_diff=0;
  for(int i=0;i<a_.size();++i){
    for(int j=0;j<a_[i].rows();++j){
      for(int k=0;k<a_[i].rows();++k){
        max_diff=max_diff>(std::abs(a_[i](j,k)-v.a_[i](j,k)))?max_diff:(std::abs(a_[i](j,k)-v.a_[i](j,k)));
        max_diff=max_diff>(std::abs(b_[i](j,k)-v.b_[i](j,k)))?max_diff:(std::abs(b_[i](j,k)-v.b_[i](j,k)));
      }
    }
  }
  return max_diff;
}
/// Compute the maximum relative difference between vertex v and *this. This is useful to see if two vertices are equal (up to roundoff errors).
double vertex::max_rel_diff(const vertex &v) const{
  double max_diff=0;
  for(int i=0;i<a_.size();++i){
    for(int j=0;j<a_[i].rows();++j){
      for(int k=0;k<a_[i].rows();++k){
        double rel_a=std::abs((a_[i](j,k)-v.a_[i](j,k))/a_[i](j,k));
        double rel_b=std::abs((b_[i](j,k)-v.b_[i](j,k))/b_[i](j,k));
        max_diff=max_diff>rel_a?max_diff:rel_a;
        max_diff=max_diff>rel_b?max_diff:rel_b;
      }
    }
  }
  return max_diff;
}

///This function combines the 4 submatrix of a nambu vertex into a whole matrix
matrix_t vertex::combine_matrix(const matrix_t &x00, const matrix_t &x01){
  int size=int(std::sqrt(x00.size()));
  matrix_t full_vertex(size*2,size*2);
  for(int i=0;i<size;i++){
    for(int j=0;j<size;j++){
      full_vertex(i,j) = x00(i,j);
      full_vertex(i,j+size) = x01(i,j);
      full_vertex(i+size,j) = std::conj(x01(i,j));
      full_vertex(i+size,j+size) = std::conj(x00(i,j));
    }
  }
  return full_vertex;
}
void vertex::separate_matrix(const matrix_t &whole_matrix,matrix_t &x00, matrix_t &x01){
  int size=int(std::sqrt(whole_matrix.size()));
  for(int i=0;i<size/2;i++){
    for(int j=0;j<size/2;j++){
      x00(i,j) = whole_matrix(i,j);
      x01(i,j) = whole_matrix(i,j+size/2);
    }
  }
}
///This function will give you back the eigenvectors and eigenvalues of A, such that
/// inv(U) * A * U = EV
/// This is validated as a unit test
void vertex::inplace_eigenvalues_eigenvectors(matrix_t &A, matrix_t &U, vector_t &EV){
  int size=EV.size();
  int wsize=4*size*size;
  std::vector<std::complex<double> > lwork(wsize);
  std::vector<double > rwork(wsize);
  char Vec='V';
  char NoVec='N';
  int info=0;
  matrix_t ignored(size, size);

  //std::cout<<"on entry zgeev"<<A.rows()<<" "<<A.num_cols()<<" "<<U.rows()<<" "<<U.cols()<<" "<<ignored.rows()<<" "<<ignored.cols()<<std::endl;
  // make sure to pass in matrices for left and right eigenvectors
  lapack::zgeev_(&NoVec,&Vec,&size, &(A(0,0)), &size, &(EV[0]), &(ignored(0,0)), &size, &(U(0,0)), &size, &(lwork[0]), &wsize,&(rwork[0]),&info);
}

std::ostream& vertex::matlab_print(std::ostream &os, const matrix_t &A){
  for(int i=0;i<A.rows();++i){
    for(int j=0;j<A.cols();++j){
      os<<A(i,j).real()<<" "<<A(i,j).imag()<<" ";
    }
    os<<std::endl;
  }
  os<<std::endl;
  return os;
}
void vertex::solve_a_eq_b_X_c(const matrix_t &A, const matrix_t &B, matrix_t &X, const matrix_t &C) {
  ///solves the equation A = B X C for X
  matrix_t Abuf(A), Bbuf(B), Cbuf(C);
  std::vector<int> ipiv(A.rows());

  gesv(Bbuf, ipiv,Abuf);
  matrix_t CtransXtrans=Abuf.transpose();
  matrix_t Ctrans=C.transpose();

  gesv(Ctrans,ipiv,CtransXtrans);
  X=CtransXtrans.transpose();
}
void vertex::solve_a_eq_b_X_c_directly(const matrix_t &A, const matrix_t &B, matrix_t &X, const matrix_t &C) {
  ///solves the equation A = B X C for X
  std::vector<int> ipiv(A.rows());
  matrix_t Cinv(C),Binv(B);
  matrix_t Ccpy(C),Bcpy(B);

  invert_matrix_inplace(Ccpy,Cinv,ipiv);
  invert_matrix_inplace(Bcpy,Binv,ipiv);

  std::complex<double> one=1., zero=0.;
  matrix_t BinvA(Cinv.cols(),Cinv.rows());
  BinvA=Binv*A;
  X=BinvA*Cinv;
}
/// Compute \f{align} \chi^{-1} = \chi_c^{-1} -(\chi_c^0)^-1 + {\overline{\chi_0^{-1}}} \f}
void vertex::triple_invert_inv_minv_pinv(matrix_t &res, const matrix_t &chic, const matrix_t &chi0c, const matrix_t &chi0bar){
  matrix_t chic_buf=chic;
  matrix_t chi0c_buf=chi0c;
  matrix_t chi0bar_buf=chi0bar;
  std::vector<int> ipiv(chic.rows());

  matrix_t chic_inv(chic.rows(), chic.rows());
  matrix_t chi0c_inv(chic.rows(), chic.rows());
  matrix_t chi0bar_inv(chic.rows(), chic.rows());

  invert_matrix_inplace(chic_buf, chic_inv,ipiv);
  invert_matrix_inplace(chi0c_buf, chi0c_inv,ipiv);
  invert_matrix_inplace(chi0bar_buf, chi0bar_inv,ipiv);
  matrix_t res_inv=chic_inv - chi0c_inv + chi0bar_inv;

  invert_matrix_inplace(res_inv, res,ipiv);

}
///Compute chi = [chic(chi0^{-1}-chi0c^{-1})+1]^{-1}*chic
void vertex::invert_oneplussmall_chic(matrix_t &res, const matrix_t &chic, const matrix_t &chi0c, const matrix_t &chi0bar){
  matrix_t chi0c_buf=chi0c;
  matrix_t chi0bar_buf=chi0bar;
  std::vector<int> ipiv(chic.rows());
  matrix_t chi0c_inv(chic.rows(), chic.rows());
  matrix_t chi0_inv(chic.rows(), chic.rows());
  matrix_t oneplussmall_inv(chic.rows(), chic.rows());
  matrix_t Identity;
  Identity.setIdentity(chic.rows(), chic.rows());
  
  invert_matrix_inplace(chi0c_buf, chi0c_inv,ipiv);
  invert_matrix_inplace(chi0bar_buf, chi0_inv,ipiv);
  
  matrix_t oneplussmall=chic*(chi0_inv-chi0c_inv)+Identity;
  invert_matrix_inplace(oneplussmall,oneplussmall_inv,ipiv);
  res = oneplussmall_inv*chic;
}


//-------------------- Magnetic Vertex -------------------------------------
std::vector<std::vector<std::complex<double> > >vertex_magnetic::compute_static_magnetic_susceptibility(int n_omega) const{
  std::vector<std::vector<std::complex<double> > >chiloc(n_sites_,std::vector<std::complex<double> >(n_omega4_bose_,0.));
  for(int nu=n_omega4_bose_;nu<n_omega4_bose_;++nu){
    for(int Q=0;Q<n_sites_;++Q){
      for(int omega=-n_omega;omega<n_omega;++omega){
        for(int K=0;K<n_sites_;++K){
          for(int omegaprime=-n_omega;omegaprime<n_omega;++omegaprime){
            for(int Kprime=0;Kprime<n_sites_;++Kprime){
              chiloc[Q][n_omega4_bose_]+=2.*(v00(K,Kprime, Q, omega,omegaprime, nu)+v01(K,Kprime, Q, omega,omegaprime, nu));
            }
          }
        }
      }
      chiloc[Q][nu+n_omega4_bose_]/=(beta_*beta_*n_sites_*n_sites_);
    }
  }
  return chiloc;
}
void vertex::broadcast_data(){
#ifdef ALPSCore_HAS_MPI
  for(int i=0;i<a_.size();++i){
    MPI_Bcast(&a_[i](0,0), a_[i].size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(&b_[i](0,0), b_[i].size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  }
#else
  return;
#endif
}
