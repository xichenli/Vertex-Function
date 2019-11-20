#include"vertexaux.h"
#include"k_space_structureTest.h"

TEST_F(k_space_structureTest, MomentaNeg) {
  for(int i=0;i<k_4site_ptr->n_sites();++i){
    EXPECT_EQ(i, k_4site_ptr->momenta_neg(k_4site_ptr->momenta_neg(i)));
  }
}
TEST_F(k_space_structureTest, MomentaSumDiff) {
  for(int i=0;i<k_4site_ptr->n_sites();++i){
    for(int j=0;j<k_4site_ptr->n_sites();++j){
      EXPECT_EQ(i, k_4site_ptr->momenta_diff(k_4site_ptr->momenta_sum(i,j),j));
      EXPECT_EQ(j, k_4site_ptr->momenta_diff(k_4site_ptr->momenta_sum(i,j),i));
    }
  }
}
TEST_F(k_space_structureTest,FindMomenta){
  std::vector<double> q(3);
  q[0]=M_PI/8;q[1]=M_PI/8;q[2]=0;
  EXPECT_EQ(k_4site_ptr->find_cluster_momentum(q),0);
  q[0]=M_PI/2+0.5;q[1]=M_PI;q[2]=0;
  EXPECT_EQ(k_4site_ptr->find_cluster_momentum(q),1);
  q[0]=-M_PI+0.01;q[1]=-M_PI+0.01;q[2]=0;
  EXPECT_EQ(k_4site_ptr->find_cluster_momentum(q),1);
}  
