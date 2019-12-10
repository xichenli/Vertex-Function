#include"chi_nambu.h"
#include"mc_vertex.h"

/// \!brief Calculation of the susceptibility \f$ \chi \f$.
///  Paper Reference. PRB 86, 125114 (2012).
/// Overview of Derivation: Eq.(5) ->  Eq.(7a) for ph or pp pieces.  Fourier transform green's functions in tau to frequencies.  Focusing on only the single particle terms, their contribution to chi is:
///
/// \f$  \frac{-1}{\beta^2}\sum\limits_{\omega_1, \omega_2} \int\limits_0^\beta d\tau_1 d\tau_2 d\tau_3 G(i\omega_1) G(i\omega_2) e^{-i\tau_1(\nu+\omega_1)} e^{i\tau_2(\nu+\omega +\omega_1)} e^{-i \tau_3(\nu^\prime+\omega+\omega_2)} \f$
///
/// This simplifies to
/// \f$  - \delta_{\nu,\omega_1} \delta_{\nu+\omega, \omega_1} \delta_{\nu^\prime+\omega,\omega_2} G(i\omega_1) G(i\omega_2) \f$
///
/// If you simplify the delta functions, you can get \f$\delta_{\omega,0} G(\nu) G(\nu^\prime)  \f$.  However, we work with a multi index which seems to make the symmetries less straightforward.  Therefore, we implement in the ph channel the full expression, which writes the frequency and momentum dependencies explicitly which is:
///
/// \f$  \sum\limits_{\omega_1, \omega_2} - \delta_{\nu,\omega_1} \delta_{K, K_1} \delta_{\nu+\omega, \omega_1} \delta_{K+Q, K_1} \delta_{\nu^\prime+\omega,\omega_2} \delta_{K^\prime+Q,K_2} G(i\omega_1, K_1) G(i\omega_2, K_2)    \f$
///
/// It is likely that symmetry issues crop up due to some subtlety of adding momenta in the brillouin zone.
//void chi_cluster_nambu::nambu_compute_chi(const mc_vertex &g4, const nambu_fermionic_green_function &g2, alps::hdf5::archive &ar){

//This function only calculate the magnetic channel of chi , which is why there's no -G2*G2.
void chi_cluster_magnetic::compute_chi(const mc_vertex &g4_phuu,const mc_vertex &g4_phud, alps::hdf5::archive &ar){
  double sign_result;
  ar>>alps::make_pvp("/simulation/results/Sign/mean/value", sign_result);
  std::cout<<"Using value of Sign= "<< sign_result<<std::endl;
  /// \todo PUT IN THE BETA CONVENTION FOR THE SINGLE PARTICLE GF
  for(int omega=-n_omega4_;omega<n_omega4_;++omega){
    for(int K=0;K<n_sites_;++K){
      for(int nu=0;nu<n_omega4_bose_;++nu){
        for(int Q=0;Q<n_sites_;++Q){
           for(int omegaprime=-n_omega4_;omegaprime<n_omega4_;++omegaprime){
            for(int Kprime=0;Kprime<n_sites_;++Kprime){
              int symmetrized_K_index     = ks_->symmetrized_K_index     (K,Kprime,Q);
              int symmetrized_Kprime_index= ks_->symmetrized_Kprime_index(K,Kprime,Q);
              int symmetrized_Q_index     = ks_->symmetrized_Q_index     (K,Kprime,Q);
              int freqindex3_buf          = g4_phuu.freqindex3(nu, omega, omegaprime);
              //construct the chi here:
              //the dimension of chi is [beta^3]
              //v01=F4(phuu)-F4(phud)
              v00(K,Kprime,Q,omega,omegaprime,nu)=g4_phuu.vertex00[symmetrized_Q_index][symmetrized_K_index][symmetrized_Kprime_index][freqindex3_buf]/sign_result-g4_phud.vertex00[symmetrized_Q_index][symmetrized_K_index][symmetrized_Kprime_index][freqindex3_buf]/sign_result;
              v01(K,Kprime,Q,omega,omegaprime,nu)=g4_phuu.vertex01[symmetrized_Q_index][symmetrized_K_index][symmetrized_Kprime_index][freqindex3_buf]/sign_result-g4_phud.vertex01[symmetrized_Q_index][symmetrized_K_index][symmetrized_Kprime_index][freqindex3_buf]/sign_result;
            }
          }
        }
      }
    }
  }
}
void chi_cluster_magnetic::compute_chi(const mc_vertex &g4_m, alps::hdf5::archive &ar){
  double sign_result;
  ar>>alps::make_pvp("/simulation/results/Sign/mean/value", sign_result);
  std::cout<<"Using value of Sign= "<< sign_result<<std::endl;
  std::cout<<"Compute Cluster Magnetic Chi"<<std::endl;
  std::ofstream file2;
  file2.open("chi_separate2.dat");
 /// \todo PUT IN THE BETA CONVENTION FOR THE SINGLE PARTICLE GF
  for(int omega=-n_omega4_;omega<n_omega4_;++omega){
    for(int K=0;K<n_sites_;++K){
      for(int nu=0;nu<n_omega4_bose_;++nu){
        for(int Q=0;Q<n_sites_;++Q){
          std::complex<double>chi_nuQ00(0.,0.);
          std::complex<double>chi_nuQ01(0.,0.);
          for(int omegaprime=-n_omega4_;omegaprime<n_omega4_;++omegaprime){
            for(int Kprime=0;Kprime<n_sites_;++Kprime){
              int symmetrized_K_index     = ks_->symmetrized_K_index     (K,Kprime,Q);
              int symmetrized_Kprime_index= ks_->symmetrized_Kprime_index(K,Kprime,Q);
              int symmetrized_Q_index     = ks_->symmetrized_Q_index     (K,Kprime,Q);
              int freqindex3_buf          = g4_m.freqindex3(nu, omega, omegaprime);
              //construct the chi here:
              //the dimension of chi is [beta^3]
              //v01=F4(phuu)-F4(phud)
              v00(K,Kprime,Q,omega,omegaprime,nu)=g4_m.vertex00[symmetrized_Q_index][symmetrized_K_index][symmetrized_Kprime_index][freqindex3_buf]/sign_result;
              v01(K,Kprime,Q,omega,omegaprime,nu)=g4_m.vertex01[symmetrized_Q_index][symmetrized_K_index][symmetrized_Kprime_index][freqindex3_buf]/sign_result;
              chi_nuQ00 += v00(K,Kprime,Q,omega,omegaprime,nu);
              chi_nuQ01 += v01(K,Kprime,Q,omega,omegaprime,nu);
            }
          }
          if(omega>=0)
            file2<<omega<<" "<<K<<" "<<nu<<" "<<Q<<" "<<std::real(chi_nuQ00)<<" "<<std::imag(chi_nuQ00)<<" "<<std::real(chi_nuQ01)<<" "<<std::imag(chi_nuQ01)<<std::endl;
        }
      }
    }
  }
  file2.close();
}
void chi_cglattice_magnetic::compute_chi_cglattice_directly(const chi0_cglattice_magnetic &chi0, const chi_cluster_magnetic &chic, const chi0_cluster_magnetic &chi0c){
  int size = int(std::sqrt(chic.v00(0,0).size()));
  for(int nu=0;nu<n_omega4_bose_;++nu){
    for(int Q=0;Q<n_sites_;++Q){
      std::vector<int> ipiv(size*2);
      matrix_t chi_all = combine_matrix(chic.v00(Q,nu),chic.v01(Q,nu)); 
      matrix_t chic_inv(size*2,size*2);
      invert_matrix_inplace(chi_all,chic_inv,ipiv);
      for(int i=0;i<size;i++){
        std::complex<double>delta_chi0 = std::norm(chi0.v00(Q,nu)(i,i))-std::norm(chi0.v01(Q,nu)(i,i));
        std::complex<double>delta_chi0c = std::norm(chi0c.v00(Q,nu)(i,i))-std::norm(chi0c.v01(Q,nu)(i,i));

        chic_inv(i,i) += std::conj(chi0.v00(Q,nu)(i,i))/delta_chi0-std::conj(chi0c.v00(Q,nu)(i,i))/delta_chi0c;
        chic_inv(i,i+size) += -chi0.v01(Q,nu)(i,i)/delta_chi0+chi0c.v01(Q,nu)(i,i)/delta_chi0c;
        chic_inv(i+size,i) += std::conj(-chi0.v01(Q,nu)(i,i)/delta_chi0+chi0c.v01(Q,nu)(i,i)/delta_chi0c);
        chic_inv(i+size,i+size) += chi0.v00(Q,nu)(i,i)/delta_chi0-chi0c.v00(Q,nu)(i,i)/delta_chi0c;
      }
      invert_matrix_inplace(chic_inv,chi_all,ipiv);
      separate_matrix(chi_all,v00(Q,nu),v01(Q,nu));
    }
  }
}

