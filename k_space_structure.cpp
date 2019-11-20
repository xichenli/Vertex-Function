//DMFT headers
#include"cluster.h"
#include"k_space_structure.h"
#include"general_matrix.h"


void k_space_structure::define_parameters(alps::params &p){
  p.define<double>("t", "Hopping element for semi circular & square lattice density of states ");
  p.define<double>("tprime", 0., "Next-nearest Hopping element for square lattice density of states ");
  p.define<double>("MU", "Chemical Potential");
  p.define<int>("N", "Number of imaginary time discretization points");
  p.define<int>("NMATSUBARA", "Number of matsubara frequency discretization points");
  p.define<int>("NOMEGA4", "Number of two-particle matsubara frequency points for fermionic frequencies");
  p.define<int>("NOMEGA4_BOSE", "Number of two-particle matsubara frequency points for bosonic frequencies");
  p.define<int>("ctaux.NOMEGA4", "Number of two-particle matsubara frequency points for fermionic frequencies, in later ctaux codes");
  p.define<int>("ctaux.NOMEGA4_BOSE", "Number of two-particle matsubara frequency points for bosonic frequencies, in later ctaux codes");
  p.define<int>("FLAVORS", 2, "Number of spins or (diagonal) orbitals");
  p.define<double>("BETA", "Inverse temperature");
  p.define<double>("H", "Staggered (antiferromagnetic) field");
  define_dca_parameters(p);
}

///find momenta sums and differences (modulo 2 PI)
void k_space_structure::find_momenta_sum_diff(){
  momenta_diff_.resize(n_sites_);
  momenta_sum_.resize(n_sites_);
  for(int i=0;i<n_sites_;++i){
    momenta_diff_[i].resize(n_sites_, -1);
    momenta_sum_[i].resize(n_sites_, -1);
  }
  double tol=1.e-14;
  for(int K=0;K<n_sites_;++K){
    for(int Kprime=0;Kprime<n_sites_;++Kprime){
      //std::cout<<"running forK "<<K<<" Kprime: "<<Kprime<<std::endl;
      double sum_x=momenta_[K][0]+momenta_[Kprime][0]; if(sum_x>M_PI+tol) sum_x-=2*M_PI; if(sum_x<-M_PI+tol) sum_x+=2*M_PI;
      double diff_x=momenta_[K][0]-momenta_[Kprime][0]; if(diff_x<-M_PI+tol) diff_x+=2*M_PI; if(diff_x>M_PI+tol) diff_x-=2*M_PI;
      double sum_y, diff_y;
      if(dimensionality_>1){
        sum_y=momenta_[K][1]+momenta_[Kprime][1]; if(sum_y>M_PI+tol) sum_y-=2*M_PI; if(sum_y<-M_PI+tol) sum_y+=2*M_PI;
        diff_y=momenta_[K][1]-momenta_[Kprime][1]; if(diff_y<-M_PI+tol) diff_y+=2*M_PI; if(diff_y>M_PI+tol) diff_y-=2*M_PI;
      }
      if(dimensionality_>2) throw std::logic_error("implement summing up and taking differences of 3d momenta.");
      for(int Q=0;Q<n_sites_;++Q){
        if(std::abs(sum_x-momenta_[Q][0])<tol && (dimensionality_==1?true:std::abs(sum_y-momenta_[Q][1])<tol))
          momenta_sum_[K][Kprime]=Q;
        if(std::abs(diff_x-momenta_[Q][0])<tol && (dimensionality_==1?true:std::abs(diff_y-momenta_[Q][1])<tol))
          momenta_diff_[K][Kprime]=Q;
      }
      if(momenta_sum_[K][Kprime]==-1){
        std::cout<<"K: "<<momenta_[K][0]<<" "<<(dimensionality_==1?0.:momenta_[K][1])<<" kprime: "<<momenta_[Kprime][0]<<" "<<(dimensionality_==1?0.:momenta_[Kprime][1])<<" sum: "<<sum_x<<" "<<(dimensionality_==1?0.:sum_y)<<std::endl;
        throw std::logic_error("problem with periodicity of cluster!");
      }
      if(momenta_diff_[K][Kprime]==-1){
        std::cout<<"K: "<<momenta_[K][0]<<" "<<(dimensionality_==1?0.:momenta_[K][1])<<" kprime: "<<momenta_[Kprime][0]<<" "<<(dimensionality_==1?0.:momenta_[Kprime][1])<<" diff: "<<diff_x<<" "<<diff_y<<std::endl;
        throw std::logic_error("problem with periodicity of cluster!");
      }
    }
  }

/*  std::cout<<"k space structure:"<<std::endl;
  for(int K=0;K<n_sites_;++K){
    std::cout<<K<<" "<<momenta_[K][0]<<" "<<momenta_[K][1]<<std::endl;
  }*/
}
//find the momenta patch a k vector belong to
int k_space_structure::find_cluster_momentum(std::vector<double>& k) const{
  double tol = 1.e-15;
  double distance=1000,tmp_distance=1000;
  int the_K=0;
  for(int K=0;K<n_sites_;++K){
    tmp_distance=(k[0]-momenta_[K][0])*(k[0]-momenta_[K][0])+(k[1]-momenta_[K][1])*(k[1]-momenta_[K][1]);
    if(tmp_distance<distance){
      the_K=K;
      distance=tmp_distance;
    }
    if(std::abs(momenta_[K][0]-M_PI)<tol||std::abs(momenta_[K][0]-M_PI)<tol){
      tmp_distance=(k[0]+M_PI)*(k[0]+M_PI)+(k[1]-momenta_[K][1])*(k[1]-momenta_[K][1]);
      if(tmp_distance<distance){
      the_K=K;
      distance=tmp_distance;
      }
      tmp_distance=(k[0]+M_PI)*(k[0]+M_PI)+(k[1]+M_PI)*(k[1]+M_PI);
      if(tmp_distance<distance){
      the_K=K;
      distance=tmp_distance;
      }
      tmp_distance=(k[0]-momenta_[K][0])*(k[0]-momenta_[K][0])+(k[1]+M_PI)*(k[1]+M_PI);
      if(tmp_distance<distance){
      the_K=K;
      distance=tmp_distance;
    }
  }
} 
return the_K;
    } 
void k_space_structure::find_zero_momentum() {
  //std::cout<<"done finding momenta sum diff."<<std::endl;
  double tol = 1.e-15;
  zero_ = -1;
  pipi_ = -1;
  for (int i = 0; i < n_sites_; ++i) {
    if (std::abs(momenta_[i][0]) < tol) {
      if ((dimensionality_ < 2) || std::abs(momenta_[i][1]) < tol) {
        if ((dimensionality_ < 3) || std::abs(momenta_[i][2]) < tol) {
          zero_ = i;
        }
      }
    }
  }
  for (int i = 0; i < n_sites_; ++i) {
    if (std::abs(momenta_[i][0]-M_PI) < tol) {
      if ((dimensionality_ < 2) || std::abs(momenta_[i][1]-M_PI) < tol) {
        if ((dimensionality_ < 3) || std::abs(momenta_[i][2]-M_PI) < tol) {
          pipi_ = i;
        }
      }
    }
  }

  if (zero_ == -1)
    throw std::logic_error("could not find momentum zero.");

  if (n_sites_!=1&&pipi_ == -1)//For paramagnetic calculation, 1 site does not have pipi vector
    throw std::logic_error("could not find AFM (pi,pi,...) vector.");
}

void k_space_structure::find_momenta_table(const NormalStateClusterTransformer &clusterhandler) {
  //std::cout<<"dimension is: "<<dimensionality_<<std::endl;
  momenta_.resize(n_sites_, std::vector<double>(dimensionality_, 0.));
  for (int i = 0; i < n_sites_; ++i) {
    momenta_[i][0] = clusterhandler.cluster_momentum(i, 0);
    if (dimensionality_ > 1)
      momenta_[i][1] = clusterhandler.cluster_momentum(i, 1);

    if (dimensionality_ > 2)
      throw std::logic_error("check that all of this indeed works for 3d.");
    //std::cout<<"momenta: "<<i<<" "; for(int j=0;j<momenta_[i].size();++j){ std::cout<< momenta_[i][j]<<" ";} std::cout<<std::endl;
  }
}

void k_space_structure::find_vertex_symmetry_map(const NormalStateClusterTransformer &clusterhandler) {
  //get the symmetries for pairs
  NormalStateClusterTransformer::pair_class_set_type pcst=clusterhandler.pair_class_set_k_space();
  for(NormalStateClusterTransformer::pair_class_set_type::const_iterator it=pcst.begin(); it!=pcst.end();++it){
    k_space_structure::pair_class_type pct;
    for(NormalStateClusterTransformer::pair_class_type::const_iterator it2=(*it).begin();  it2!=(*it).end();++it2){
       pct.insert(*it2);
    }
    pair_class_set_k_space_.insert(pct);
  }
  vertex_symmetry_map_.resize(n_sites_*n_sites_*n_sites_,n_sites_*n_sites_*n_sites_);
  vertex_multiplicity_map_.resize(n_sites_*n_sites_*n_sites_,0);
  //std::cout<<"finding equivalent k-vectors"<<std::endl;
  int symmetry_counter = 0;
  for (int K = 0; K < n_sites_; ++K) {
    for (int Kprime = 0; Kprime < n_sites_; ++Kprime) {
      for (int Q = 0; Q < n_sites_; ++Q) {
        for (std::size_t s = 0; s < symmetry_table_k_space_.rows(); ++s) {
          vertex_symmetry_map_[K * n_sites_ * n_sites_ + Kprime * n_sites_ + Q] =
              std::min(
                  vertex_symmetry_map_[K * n_sites_ * n_sites_
                      + Kprime * n_sites_ + Q],
                  symmetry_table_k_space_(s, K) * n_sites_ * n_sites_
                      + symmetry_table_k_space_(s, Kprime) * n_sites_
                      + symmetry_table_k_space_(s, Q));
        }
        //vertex_symmetry_map_[K*n_sites_*n_sites_+Kprime*n_sites_+Q]=K*n_sites_*n_sites_+Kprime*n_sites_+Q; //NO SYMMETRIES
        vertex_multiplicity_map_[vertex_symmetry_map_[K * n_sites_ * n_sites_
            + Kprime * n_sites_ + Q]]++;
        if (vertex_symmetry_map_[K * n_sites_ * n_sites_ + Kprime * n_sites_ + Q]
            == K * n_sites_ * n_sites_ + Kprime * n_sites_ + Q)
          symmetry_counter++;

        // Debug, no vertex symmetry
//        vertex_symmetry_map_[K * n_sites_ * n_sites_ + Kprime * n_sites_ + Q] = K * n_sites_ * n_sites_ + Kprime * n_sites_ + Q;
 //       vertex_multiplicity_map_[K * n_sites_ * n_sites_ + Kprime * n_sites_ + Q] = 1;
  //      symmetry_counter = n_sites_*n_sites_*n_sites_;
      }
    }
  }
}

///find the internal symmetries of the vertex. also initialize the index of the zero momentum (0,0,....,0)
void k_space_structure::find_vertex_symmetries(){
  //std::cout<<"computing vertex function symmetry table."<<std::endl;
  NormalStateClusterTransformer clusterhandler(p_);
  //std::cout<<"symmetry table k space."<<std::endl;
  //std::cout<<clusterhandler->symmetry_table_k_space()<<std::endl;
  //get the symmetry table
  const blas::general_matrix<int> symmetry_table_k_space=clusterhandler.symmetry_table_k_space();
  symmetry_table_k_space_.resize(symmetry_table_k_space.size1(), symmetry_table_k_space.size2());
  for(int i=0;i<symmetry_table_k_space.size1();++i){
    for(int j=0;j<symmetry_table_k_space.size2();++j){
      symmetry_table_k_space_(i,j)=symmetry_table_k_space(i,j);
    }
  }

  ///creating momenta table
  dimensionality_=clusterhandler.dimension();
  //std::cout<<"dimension is: "<<dimensionality_<<std::endl;
  find_momenta_table(clusterhandler);
  //std::cout<<"finding momenta sum diff."<<std::endl;
  find_momenta_sum_diff();
  //std::cout<<"done finding momenta sum diff."<<std::endl;
  find_zero_momentum();
  //std::cout<<"finding equivalent k-vectors"<<std::endl;
  find_vertex_symmetry_map(clusterhandler);
  //std::cout<<"done finding equivalent k-vectors"<<std::endl;
}

void k_space_structure::precompute_dispersion_and_weights(std::vector<double> &weight_kl, std::vector<std::vector<double> > &epsilon_kckl, std::vector<std::vector<double> > &epsilon_kckl_plus_pipi,std::vector<std::vector<double> > &symmetry_kckl) const{
  NormalStateClusterTransformer ct(p_);
  epsilon_kckl.resize(n_sites_, std::vector<double>(ct.n_lattice_momenta()));
  epsilon_kckl_plus_pipi.resize(n_sites_, std::vector<double>(ct.n_lattice_momenta()));
  symmetry_kckl.resize(n_sites_, std::vector<double>(ct.n_lattice_momenta()));
  weight_kl.resize(ct.n_lattice_momenta());

  for (site_t kc=0; kc<n_sites_; ++kc) {
    for (site_t kl=0; kl<ct.n_lattice_momenta(); ++kl) {
      std::vector<double> k(ct.dimension());
      std::vector<double> k_plus_pipi(ct.dimension());
      for (unsigned int d=0; d<ct.dimension(); ++d){
        k[d] = ct.lattice_momentum(kl,d) + ct.cluster_momentum(kc,d);
        k_plus_pipi[d] = ct.lattice_momentum(kl,d) + ct.cluster_momentum(kc,d)+M_PI;
      }
      epsilon_kckl[kc][kl]=ct.epsilon(k);
      epsilon_kckl_plus_pipi[kc][kl]=ct.epsilon(k_plus_pipi);
      symmetry_kckl[kc][kl]=cos(k[0])-cos(k[1]);
    }
  }

  for (site_t kl=0; kl<ct.n_lattice_momenta(); ++kl) {
    weight_kl[kl]= ct.lattice_momentum(kl,ct.dimension());
    //this is dangerous and may fail if we don't have 'L' explicitly specified in the parameter file.
    double L=p_["dca.L"];
    weight_kl[kl]/=L*L; //normalization to one.
  }

}


