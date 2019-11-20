#ifndef K_SPACE_STRUCTURE_TEST_H
#define K_SPACE_STRUCTURE_TEST_H
#include"gtest/gtest.h"
#include"k_space_structure.h"


class k_space_structureTest:public ::testing::Test{
protected:
  virtual void SetUp() {
    alps::params p("UnitTests/sample_4site.param");
    k_space_structure::define_parameters(p);
    k_4site_ptr=new k_space_structure(p);
  }
  virtual void TearDown() {
    delete k_4site_ptr;
  }
  
  k_space_structure *k_4site_ptr;
};
#endif //vertexaux
