#include"vertexaux.h"
#include"blasheader.h"

void invert_matrix_inplace(matrix_t &matrix_to_be_inverted, matrix_t &inverted_matrix, std::vector<int> &ipiv){
  inverted_matrix=matrix_t::Identity(matrix_to_be_inverted.rows(),matrix_to_be_inverted.rows());
  ipiv.resize(matrix_to_be_inverted.rows());
  gesv(matrix_to_be_inverted ,ipiv,inverted_matrix);
}

void gemm(std::complex<double> &a, const matrix_t &A, const matrix_t &B, std::complex<double> &b, matrix_t &C){
  char notrans='N';
  int size=A.rows();
  //std::cout<<"on entry zgemm"<<A.num_rows()<<" "<<A.num_cols()<<" "<<B.num_rows()<<" "<<B.num_cols()<<" "<<C.num_rows()<<" "<<C.num_cols()<<std::endl;
  blas::zgemm_(&notrans, &notrans, &size, &size, &size, &a, &(A(0,0)), &size, &(B(0,0)), &size, &b,&(C(0,0)), &size);
  //std::cout<<"on exit zgemm"<<std::endl;
}

void gesv(matrix_t &A, std::vector<int> &ipiv, matrix_t &B){
    int size=A.rows();
    int info;
    lapack::zgesv_(&size, &size, &(A(0,0)), &size, &(ipiv[0]), &(B(0,0)), &size, &info);
    if(info !=0) throw std::runtime_error("problem in gesv, info !=0");
}
