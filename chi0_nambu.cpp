#include"chi0_nambu.h"
#include"vertexaux.h"
#include <sstream>
#include <string>

/// \!brief Calculation of the bare cluster susceptibility \f$ \chi_0 \f$.
///
/// For the particle-hole channel see the explicit derivation in Alessandro Eq. 10.
/// The particle particle channel is not documented there, but should be \f$ G(K)G(Q-K)\delta_{KK'} \f$.
/// \image latex chi0_ph.pdf " chi_0  in the particle-hole notation" width=10cm
/// \image html chi0_ph.jpeg " chi_0  in the particle-hole notation" width=10cm
/// \image latex chi0_pp.pdf " chi_0  in the particle-particle notation" width=10cm
/// \image html chi0_pp.jpeg " chi_0  in the particle-particle notation" width=10cm
///
/// the dimension of \f$ \chi \f$ is \f$\frac{1}{\text{energy}^3} \f$.

/// \!brief Calculation of some random initial starting values, mostly for debug purposes so that we don't need to read in actual data.
chi0_cglattice_magnetic::chi0_cglattice_magnetic(const alps::params &p, const k_space_structure &ks):chi0_magnetic(p,ks){
  mu_=p["MU"];
  precompute_integration_weights();
}

/// \!brief computation of k-summation integration weights
/// This function computes all the epsilon (k) and their (Simpson) weights for the k-space integration. This is done by extracting them out of the clustertransformer, which is stored as a private variable.
void chi0_cglattice_magnetic::precompute_integration_weights(){
  //buffer lattice momenta and weights
  ks_->precompute_dispersion_and_weights(weight_kl_, epsilon_kckl_,epsilon_kckl_plus_pipi_,symmetry_kckl_);

}
chi0_cglattice_magnetic::~chi0_cglattice_magnetic(){
}

/// \!brief calculation of the coarse-grained bare lattice susceptibility.
/// This is done according to Thomas maier's review, equation 127, but note that we do this ONLY for Q vectors which are on cluster vectors. In principle we could generalize this to any q-vector.
///
/// The equation is: \f$ \chi^{0,ph})_{\sigma \sigma'}(K, Kprime, Q, omega, omegaprime, nu) = \delta_{\sigma \sigma'} \delta{KK'}\delta_{\omega\omega'} \frac{N_c}{N} \sum_{\tilde k\in K} G(K+\tilde{k},\omega)G(K+Q+\tilde k, \omega+\nu) \f$
/// This function will need to be tested thoroughly.... how?
/// the dimension of this quantity is beta^3, same as the dimension of the other chis.


void chi0_cluster_magnetic::compute_chi0(const nambu_fermionic_green_function &g2){
  //this gets the dimensions right: GF dimension is 1/energy if single frequency convention, 1/energy^2 if double frequency convention.
  double prefactor=g2.conv()==beta_conv_single_freq?beta_*n_sites_:1./beta_/n_sites_;
  std::cout<<"compute cluster chi0. prefactor is: "<<prefactor<<" beta is: "<<beta_<<std::endl;
  //std::cout<<"g 0 0 is: "<<g2(0,0,spinup)<<" "<<g2(0,0,spindn)<<std::endl;
  for(int nu=-n_omega4_bose_;nu<=n_omega4_bose_;++nu){
    for(int Q=0;Q<n_sites_;++Q){
      for(int omega=-n_omega4_;omega<n_omega4_;++omega){
        for(int K=0;K<n_sites_;++K){
          for(int omegaprime=-n_omega4_;omegaprime<n_omega4_;++omegaprime){
            for(int Kprime=0;Kprime<n_sites_;++Kprime){
              if(omega==omegaprime && K == Kprime){
                int KpQ=ks_->momenta_sum(K,Q);
                int QmK=ks_->momenta_diff(Q,K);
                v00(K,Kprime,Q,omega,omegaprime,nu)=-g2(omega,K,0,0)*g2(omega+nu,KpQ,0,0)*prefactor;
                v01(K,Kprime,Q,omega,omegaprime,nu)=-g2(omega,K,1,0)*g2(omega+nu,KpQ,0,1)*prefactor;
              }
            }
          }
        }
      }
    }
  }
  for(int Q=2;Q<=3;Q++){
    for(int nu=-1;nu<=1;nu++){
      std::stringstream name00,name01;
      name00<<"chi0_m00_Q"<<Q<<"nu"<<nu<<".dat";
      name01<<"chi0_m01_Q"<<Q<<"nu"<<nu<<".dat";
      std::ofstream file_m00,file_m01;
      file_m00.open(name00.str().c_str());file_m01.open(name01.str().c_str());
      vertex::matlab_print(file_m00,v00(Q,nu));
      vertex::matlab_print(file_m01,v01(Q,nu));
      file_m00.close();file_m01.close();
    }
  }
}

void chi0_cglattice_magnetic::compute_chi0(const nambu_fermionic_self_energy &sigma){
  double prefactor=beta_*n_sites_;
  for(int nu=-n_omega4_bose_;nu<=n_omega4_bose_;++nu){
    for(int Q=0;Q<n_sites_;++Q){
      for(int omega=-n_omega4_;omega<n_omega4_;++omega){
        for(int K=0;K<n_sites_;++K){
          for(int omegaprime=-n_omega4_;omegaprime<n_omega4_;++omegaprime){
            for(int Kprime=0;Kprime<n_sites_;++Kprime){
              if(omega==omegaprime && K == Kprime){
                
                v00(K,Kprime,Q,omega,omegaprime,nu)=0.;
                v01(K,Kprime,Q,omega,omegaprime,nu)=0.;
                
                //ph and pp
                std::complex<double> iomega_n(0, (2.*omega+1.)*M_PI/beta_);
                std::complex<double> zeta_omega = iomega_n + mu_;
                
                //ph only
                std::complex<double> i_nupomega(0, (2.*(omega+nu)+1.)*M_PI/beta_);
                std::complex<double> zeta_omegapnu = i_nupomega + mu_;
                int KpQ=ks_->momenta_sum(K,Q);
                
                //pp only
                int QmK=ks_->momenta_diff(Q,K);
                std::complex<double> i_numomega(0, (2.*(nu-omega-1)+1.)*M_PI/beta_);
                std::complex<double> zeta_numomega = i_numomega + mu_;
                
                double weight_sum=0.;
                double symmetry_sum=0.;
                for (int kl=0; kl<weight_kl_.size(); ++kl) {
                  std::complex<double> G1_omega_latt_inv[4]=
                  {iomega_n+mu_-epsilon_kckl_[K][kl]-sigma(omega,K,0,0),
                    -sigma(omega,K,0,1),
                    -sigma(omega,K,1,0),
                    iomega_n-mu_+epsilon_kckl_[K][kl]-sigma(omega,K,1,1)};
                  std::complex<double> G1_latt_inv_det=G1_omega_latt_inv[0]*G1_omega_latt_inv[3]-G1_omega_latt_inv[1]*G1_omega_latt_inv[2];
                  
                  std::complex<double> G2_omega_latt_inv[4]=
                  {i_nupomega+mu_-epsilon_kckl_[KpQ][kl]-sigma(omega+nu,KpQ,0,0),
                    -sigma(omega+nu,KpQ,0,1),
                    -sigma(omega+nu,KpQ,1,0),
                    i_nupomega-mu_+epsilon_kckl_[KpQ][kl]-sigma(omega+nu,KpQ,1,1)};
                  std::complex<double> G2_latt_inv_det=G2_omega_latt_inv[0]*G2_omega_latt_inv[3]-G2_omega_latt_inv[1]*G2_omega_latt_inv[2];
                  
                  v00(K,Kprime,Q,omega,omegaprime,nu) += weight_kl_[kl]*G1_omega_latt_inv[3]*G2_omega_latt_inv[3]/G1_latt_inv_det/G2_latt_inv_det;
                  v01(K,Kprime,Q,omega,omegaprime,nu) += weight_kl_[kl]*G1_omega_latt_inv[2]*G2_omega_latt_inv[1]/G1_latt_inv_det/G2_latt_inv_det;
                  weight_sum+=weight_kl_[kl];
                }

                if(std::abs(weight_sum -1) > 1.e-12) throw std::logic_error("problem with the weight sum"); //this should be a unit test
                v00(K,Kprime,Q,omega,omegaprime,nu) *=-prefactor;
                v01(K,Kprime,Q,omega,omegaprime,nu) *=-prefactor;
              }
            }
          }
        }
      }
    }
  }
}
