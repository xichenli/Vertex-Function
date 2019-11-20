#ifndef VERTEXAUX_H
#define VERTEXAUX_H
#include<iostream>
#include<complex>
#include<fstream>
#include<vector>
#include<cmath>
#include<complex>
#include <Eigen/Dense>

//#include "matrix.h"
typedef Eigen::Matrix< std::complex<double> , Eigen::Dynamic , Eigen::Dynamic > matrix_t;
typedef Eigen::Matrix< int , Eigen::Dynamic , Eigen::Dynamic > int_matrix_t;
typedef std::vector<std::complex<double> > vector_t;

//ALPS headers
#include"alps/params.hpp"


//Green's functions
//#include"k_space_structure.h"
//#include"single_freq_gf.h"
//#include"gfpower.h"

typedef std::vector<std::vector<std::complex<double> > > onefreq_gf_t;
typedef std::vector<std::vector<std::vector<std::complex<double> > > > twofreq_gf_t;
typedef std::vector<std::vector<std::vector<std::vector<std::complex<double> > > > > vertex_t;
typedef std::vector<matrix_t> vertex_matrix_t;


//forward declarations
class vertex_dmst;
class Gamma_dmst;
class chi0_dmst;
class chi0_phpp;
class chi0_dmst;
class chi_phpp;
class chi_dmst;
class chi_cluster_dmst;
class chi0_cluster_dmst;
class F_dmst;
class vertexTest;
class k_space_structure;
class mc_vertex;
class chi_cglattice_dmst;
class chi0_cglattice_dmst;
class k_space_structure;


typedef std::vector<int> index_map_t;

void invert_matrix_inplace(matrix_t &matrix_to_be_inverted, matrix_t &inverted_matrix, std::vector<int> &ipiv);
void gesv(matrix_t &A, std::vector<int> &ipiv, matrix_t &B);
void gemm(std::complex<double> &a, const matrix_t &A, const matrix_t &B, std::complex<double> &b, matrix_t &C);
#endif //vertexaux
