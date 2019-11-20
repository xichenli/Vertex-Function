#include "mc_vertex.h"
#include "k_space_structure.h"
#include "blasheader.h"

mc_vertex::mc_vertex(const alps::params &p, const k_space_structure &ks):ks_(&ks){
  //read in basic parameters for vertices
  n_omega4_=p["NOMEGA4_ORIG"];
  n_omega4_bose_=p["NOMEGA4_BOSE_ORIG"];
  n_sites_=p["dca.SITES"];
  beta_=p["BETA"];

  vertex00.resize(n_sites_,twofreq_gf_t(n_sites_, std::vector<std::vector<std::complex<double> > >(n_sites_)));
  vertex01.resize(n_sites_,twofreq_gf_t(n_sites_, std::vector<std::vector<std::complex<double> > >(n_sites_)));
}
void mc_vertex::set_parameters(alps::params &p){
  p.define<int>("inverter.NOMEGA4_ORIG", "number of fourpoint frequencies in input file");
  p.define<int>("inverter.NOMEGA4", "number of fourpoint frequencies to be used in code");
  p.define<int>("inverter.NOMEGA4_BOSE_ORIG", "number of bosonic frequencies in input file");
  p.define<int>("inverter.NOMEGA4_BOSE", "number of bosonic frequencies to be used in code");
  p.define<int>("dca.SITES", "number of cluster sites");
  p.define<double>("BETA", "inverse temperature");
}


void mc_vertex::read(alps::hdf5::archive &ar,std::string type,std::string e_or_v){

  /*  std::cout<<"reading sign from sim.h5 file"<<std::endl;
  double sign_result;
  ar>>alps::make_pvp("/simulation/results/Sign/mean/value", sign_result);
  std::cout<<"Using value of Sign= "<< sign_result<<std::endl;
   */
  int tot_freq_size=n_omega4_bose_*(2*n_omega4_)*(2*n_omega4_);
  std::cout<<"tot_freq_size"<<tot_freq_size<<std::endl;
  std::vector<double> buffer(tot_freq_size);
  for(int Q=0;Q<n_sites_;++Q){
    for(int K=0;K<n_sites_;++K){
      for(int Kprime=0;Kprime<n_sites_;++Kprime){
        int symm_Kprime_index=ks_->symmetrized_Kprime_index(K,Kprime,Q);
        int symm_K_index     =ks_->symmetrized_K_index(K,Kprime,Q);
        int symm_Q_index     =ks_->symmetrized_Q_index(K,Kprime,Q);
        vertex00[symm_Q_index][symm_K_index][symm_Kprime_index].resize(tot_freq_size);
        vertex01[symm_Q_index][symm_K_index][symm_Kprime_index].resize(tot_freq_size);
      }
    }
  }
  for(int Q=0;Q<n_sites_;++Q){
    for(int K=0;K<n_sites_;++K){
      for(int Kprime=0;Kprime<n_sites_;++Kprime){
        if(ks_->vertex_multiplicity(K,Kprime,Q)!=0){
          std::stringstream vertex00_re_name, vertex01_re_name;
          std::stringstream vertex00_im_name, vertex01_im_name;
          vertex00_re_name<<"/simulation/results/G4_Q_K_Kprime_nu_omega_omegaprime_"<<type<<"_re_"<<Q<<"_"<<K<<"_"<<Kprime<<"_times_sign/mean/"<<e_or_v;
          vertex00_im_name<<"/simulation/results/G4_Q_K_Kprime_nu_omega_omegaprime_"<<type<<"_im_"<<Q<<"_"<<K<<"_"<<Kprime<<"_times_sign/mean/"<<e_or_v;
          vertex01_re_name<<"/simulation/results/F4_Q_K_Kprime_nu_omega_omegaprime_"<<type<<"_re_"<<Q<<"_"<<K<<"_"<<Kprime<<"_times_sign/mean/"<<e_or_v;
          vertex01_im_name<<"/simulation/results/F4_Q_K_Kprime_nu_omega_omegaprime_"<<type<<"_im_"<<Q<<"_"<<K<<"_"<<Kprime<<"_times_sign/mean/"<<e_or_v;
          //std::cout<<vertex00_re_name.str()<<" "<<vertex00_im_name.str()<<" "<<vertex01_re_name.str()<<" "<<vertex01_im_name.str()<<std::endl;
          int one=1, two=2;
          double *real_location, *imag_location;
         //00
          ar>>alps::make_pvp(vertex00_re_name.str(), buffer);
          real_location=(double*)&(vertex00[Q][K][Kprime][0]);
          imag_location=real_location+1;
          blas::dcopy_(&tot_freq_size, &(buffer[0]),&one, real_location,&two);
          ar>>alps::make_pvp(vertex00_im_name.str(), buffer);
          blas::dcopy_(&tot_freq_size, &(buffer[0]),&one, imag_location,&two);
         //01
          ar>>alps::make_pvp(vertex01_re_name.str(), buffer);
          real_location=(double*)&(vertex01[Q][K][Kprime][0]);
          imag_location=real_location+1;
          blas::dcopy_(&tot_freq_size, &(buffer[0]),&one, real_location,&two);
          ar>>alps::make_pvp(vertex01_im_name.str(), buffer);
          blas::dcopy_(&tot_freq_size, &(buffer[0]),&one, imag_location,&two);
        }
      }
    }
  }
}
