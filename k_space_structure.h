#ifndef K_SPACE_STRUCTURE_H
#define K_SPACE_STRUCTURE_H
#include<set>
#include"vertexaux.h"

class NormalStateClusterTransformer;
/// \brief A class for describing the momentum structure.
///
/// This class knows about how to add and subtract momenta. The necessary info is extracted from the ALPS DCA framework.
/// The class can also handle the symmetry classes of the vertex and green's functions – this, too is extracted from the ALPS DMFT framework.
class k_space_structure{
public:
  //types
  typedef unsigned int site_t;
  typedef std::set<site_t> site_class_type;
  typedef std::pair<site_t, site_t> pair_type;
  typedef std::set<pair_type> pair_class_type;
  typedef std::set<pair_class_type> pair_class_set_type;

  k_space_structure(alps::params &p):p_(p){
    n_sites_=p["dca.SITES"];
    find_vertex_symmetries();
  }
  static void define_parameters(alps::params &p);

  /// return P + Q
  int momenta_sum (int P,int Q)const{ return momenta_sum_ [P][Q];}
  /// return P - Q
  int momenta_diff(int P,int Q)const{ return momenta_diff_[P][Q];}
  ///return -K
  int momenta_neg(int P) const{ return momenta_diff_[zero_][P];}
  int symmetrized_K_index     (int K, int Kprime, int Q) const {return vertex_symmetry_map_[K*n_sites_*n_sites_+Kprime*n_sites_+Q]/(n_sites_*n_sites_);}
  int symmetrized_Kprime_index(int K, int Kprime, int Q) const {return (vertex_symmetry_map_[K*n_sites_*n_sites_+Kprime*n_sites_+Q]/n_sites_)%n_sites_;}
  int symmetrized_Q_index     (int K, int Kprime, int Q) const {return vertex_symmetry_map_[K*n_sites_*n_sites_+Kprime*n_sites_+Q]%n_sites_;}
  int vertex_multiplicity     (int K, int Kprime, int Q) const {return vertex_multiplicity_map_[K*n_sites_*n_sites_+Kprime*n_sites_+Q];}
  
  ///the momentum (0,0) is special. This gives the index that points to this momentum.
  int zero_momentum() const{return zero_;}
  int pipi_momentum() const{return pipi_;}
  
  ///access function for pair_class_set_k_space used for symmetrization of, e.g., Green's functions.
  const pair_class_set_type &pair_class_set_k_space() const{return pair_class_set_k_space_;}
  
  ///return number of sites
  int n_sites() const{return n_sites_;}
  double momentum(int K,int d)const{return momenta_[K][d];}
  void precompute_dispersion_and_weights(std::vector<double> &weights, std::vector<std::vector<double> > &epsilon_kckl, std::vector<std::vector<double> > &epsilon_kckl_plus_pipi,std::vector<std::vector<double> > &symmetry_kckl) const;
  int find_cluster_momentum(std::vector<double>& k) const;
private:
  void find_vertex_symmetries();
  void find_momenta_sum_diff();
  void find_zero_momentum();
  void find_momenta_table(const NormalStateClusterTransformer& clusterhandler);
  void find_vertex_symmetry_map(const NormalStateClusterTransformer &clusterhandler);
  std::vector<std::vector<int> > momenta_sum_;
  std::vector<std::vector<int> > momenta_diff_;
  std::vector<int> vertex_symmetry_map_;
  std::vector<int> vertex_multiplicity_map_;
  std::vector<std::vector<double> > momenta_;
  
  int_matrix_t symmetry_table_k_space_;
  //site_class_set_type site_class_set_k_space_;
  pair_class_set_type pair_class_set_k_space_;

  int n_sites_;
  int dimensionality_;
  int zero_;
  int pipi_;

  const alps::params p_;
};

#endif //vertexaux
