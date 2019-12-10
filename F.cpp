#include"F_nambu.h"
#include"chi_nambu.h"
#include"chi0_nambu.h"
//#include"Gamma_dmst.h"
#include <fstream>
#include <iostream>
#include<sstream>
template <typename T>std::string to_string(T value)
{
   std::ostringstream os ;
   os << value ;
   return os.str() ;
}
void F_magnetic::invert_F_magnetic(const chi_magnetic &chi, const chi0_magnetic &chi0){
  double beta2=beta_*beta_*n_sites_*n_sites_;
  for(int nu=0;nu<n_omega4_bose_;++nu){
    for(int Q=0;Q<n_sites_;++Q){
      std::vector<int> ipiv(4*n_sites_*n_omega4_);
      matrix_t chi_m=combine_matrix(chi.v00(Q,nu),chi.v01(Q,nu));
      matrix_t chi0_m=combine_matrix(chi0.v00(Q,nu),chi0.v01(Q,nu));
      matrix_t betasq_chi0_m_chi_m=beta2*(chi0_m-chi_m);
      chi_m.resize(0,0);
      gesv(chi0_m,ipiv,betasq_chi0_m_chi_m);
      matrix_t chi0_m_trans=(chi0_m).transpose();
      chi0_m.resize(0,0);
      matrix_t chi0_F_trans=(betasq_chi0_m_chi_m).transpose();
      gesv(chi0_m_trans,ipiv,chi0_F_trans);
      chi0_m_trans.resize(0,0);
      matrix_t F_m=(chi0_F_trans).transpose();
      chi0_F_trans.resize(0,0);
      separate_matrix(F_m,v00(Q,nu),v01(Q,nu));
    }
  }
  std::ofstream file_m00,file_m01;
  file_m00.open("F_m00.dat");file_m01.open("F_m01.dat");
  vertex::matlab_print(file_m00,v00(0,0));
  vertex::matlab_print(file_m01,v01(0,0));
  file_m00.close();file_m01.close();
}
//Calculate selfenergy by equation of motion (does not include constant Hatree term Un/2)
void F_magnetic::equation_of_motion(const nambu_fermionic_green_function &g2){
  std::cout<<"Equation of Motion"<<std::endl;
  double prefactor = U_/(n_sites_*n_sites_*beta_*beta_);
  std::ofstream EOMfile;
  EOMfile.open("selfenergy_EOM2.dat");
  std::ofstream file2;
  file2.open("selfenergy_separate2.dat");

  for(int omega=0;omega<n_omega4_;++omega){
    EOMfile<<(2.0*omega+1)*M_PI/beta_<<" ";
    for(int K=0;K<n_sites_;++K){
      std::complex<double> sigma_00(0.,0.);
      std::complex<double> sigma_01(0.,0.);
      for(int nu=-n_omega4_bose_+1;nu<n_omega4_bose_;++nu){
        for(int Q=0;Q<n_sites_;++Q){
          int KpQ=ks_->momenta_sum(K,Q);
          std::complex<double> sigma_nuQ00(0.,0.);
          std::complex<double> sigma_nuQ01(0.,0.);
          for(int omegaprime=-n_omega4_;omegaprime<n_omega4_;++omegaprime){
            for(int Kprime=0;Kprime<n_sites_;++Kprime){
              int KprimepQ = ks_->momenta_sum(Kprime,Q);
              std::complex<double>tmp(0.,0.);
              if(nu<0){
                tmp = std::conj(v00(K,Kprime,Q,-omega-1,-omegaprime-1,-nu))*(g2(omegaprime,Kprime,0,0)*g2(omegaprime+nu,KprimepQ,0,0)+g2(omegaprime,Kprime,1,0)*g2(omegaprime+nu,KprimepQ,0,1))
                     +std::conj(v01(K,Kprime,Q,-omega-1,-omegaprime-1,-nu))*(g2(omegaprime,Kprime,1,1)*g2(omegaprime+nu,KprimepQ,1,1)+g2(omegaprime,Kprime,0,1)*g2(omegaprime+nu,KprimepQ,1,0));
              }else{
                tmp = v00(K,Kprime,Q,omega,omegaprime,nu)*(g2(omegaprime,Kprime,0,0)*g2(omegaprime+nu,KprimepQ,0,0)+g2(omegaprime,Kprime,1,0)*g2(omegaprime+nu,KprimepQ,0,1))
                     +v01(K,Kprime,Q,omega,omegaprime,nu)*(g2(omegaprime,Kprime,1,1)*g2(omegaprime+nu,KprimepQ,1,1)+g2(omegaprime,Kprime,0,1)*g2(omegaprime+nu,KprimepQ,1,0));
              } 
              sigma_00 += tmp*g2(omega+nu,KpQ,0,0);
              sigma_01 += tmp*g2(omega+nu,KpQ,1,0);
              sigma_nuQ00 += tmp*g2(omega+nu,KpQ,0,0)*prefactor;
              sigma_nuQ01 += tmp*g2(omega+nu,KpQ,1,0)*prefactor;
/*
              std::complex<double> line1,line2,line3,line4;
              if(nu<0){
                line1 = std::conj(v00(K,Kprime,Q,-omega-1,-omegaprime-1,-nu))*g2(omegaprime,Kprime,0,0)*g2(omegaprime+nu,KprimepQ,0,0)*g2(omega+nu,KpQ,1,0);
                line2 = std::conj(v00(K,Kprime,Q,-omega-1,-omegaprime-1,-nu))*g2(omegaprime,Kprime,1,0)*g2(omegaprime+nu,KprimepQ,0,1)*g2(omega+nu,KpQ,1,0);
                line3 = std::conj(v01(K,Kprime,Q,-omega-1,-omegaprime-1,-nu))*g2(omegaprime,Kprime,1,1)*g2(omegaprime+nu,KprimepQ,1,1)*g2(omega+nu,KpQ,1,0);
                line4 = std::conj(v01(K,Kprime,Q,-omega-1,-omegaprime-1,-nu))*g2(omegaprime,Kprime,0,1)*g2(omegaprime+nu,KprimepQ,1,0)*g2(omega+nu,KpQ,1,0);
              }
              else{
                line1 = v00(K,Kprime,Q,omega,omegaprime,nu)*g2(omegaprime,Kprime,0,0)*g2(omegaprime+nu,KprimepQ,0,0)*g2(omega+nu,KpQ,1,0);
                line2 = v00(K,Kprime,Q,omega,omegaprime,nu)*g2(omegaprime,Kprime,1,0)*g2(omegaprime+nu,KprimepQ,0,1)*g2(omega+nu,KpQ,1,0);
                line3 = v01(K,Kprime,Q,omega,omegaprime,nu)*g2(omegaprime,Kprime,1,1)*g2(omegaprime+nu,KprimepQ,1,1)*g2(omega+nu,KpQ,1,0);
                line4 = v01(K,Kprime,Q,omega,omegaprime,nu)*g2(omegaprime,Kprime,0,1)*g2(omegaprime+nu,KprimepQ,1,0)*g2(omega+nu,KpQ,1,0);
              }
              if(K==6&&omega==0)
                file2<<nu<<" "<<Q<<" "<<omegaprime<<" "<<Kprime<<" "<<line1<<line2<<line3<<line4<<std::endl;
*/
/* 
              if(K==6&&omega==0&&nu<0)
                file2<<nu<<" "<<Q<<" "<<omegaprime<<" "<<Kprime<<" "<<std::conj(v00(K,Kprime,Q,-omega-1,-omegaprime-1,-nu))<<" "<<std::conj(v01(K,Kprime,Q,-omega-1,-omegaprime-1,-nu))<<" "<<tmp<<" "<<tmp*g2(omega+nu,KpQ,1,0)<<std::endl;
              if(K==6&&omega==0&&nu>=0)
                file2<<nu<<" "<<Q<<" "<<omegaprime<<" "<<Kprime<<" "<<v00(K,Kprime,Q,omega,omegaprime,nu)<<" "<<v01(K,Kprime,Q,omega,omegaprime,nu)<<" "<<tmp<<" "<<tmp*g2(omega+nu,KpQ,1,0)<<std::endl;
              if(K==6&&omega==-1&&nu<0)
                file3<<nu<<" "<<Q<<" "<<omegaprime<<" "<<Kprime<<" "<<std::conj(v00(K,Kprime,Q,-omega-1,-omegaprime-1,-nu))<<" "<<std::conj(v01(K,Kprime,Q,-omega-1,-omegaprime-1,-nu))<<" "<<tmp<<" "<<tmp*g2(omega+nu,KpQ,1,0)<<std::endl;
              if(K==6&&omega==-1&&nu>=0)
                file3<<nu<<" "<<Q<<" "<<omegaprime<<" "<<Kprime<<" "<<v00(K,Kprime,Q,omega,omegaprime,nu)<<" "<<v01(K,Kprime,Q,omega,omegaprime,nu)<<" "<<tmp<<" "<<tmp*g2(omega+nu,KpQ,1,0)<<std::endl;
*/
            }
          }
          file2<<omega<<" "<<K<<" "<<nu<<" "<<Q<<" "<<std::real(sigma_nuQ00)<<" "<<std::imag(sigma_nuQ00)<<" "<<std::real(sigma_nuQ01)<<" "<<std::imag(sigma_nuQ01)<<std::endl;
        }
      }
      sigma_00 *= prefactor;
      sigma_01 *= prefactor;
      EOMfile<<sigma_00.real()<<" "<<sigma_00.imag()<<" "<<sigma_01.real()<<" "<<sigma_01.imag()<<" ";
    }
    EOMfile<<"\n";
  }
  EOMfile.close();
  file2.close();
}

