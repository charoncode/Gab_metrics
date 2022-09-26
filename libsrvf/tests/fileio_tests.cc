#include <boost/test/unit_test.hpp>

#include <srvf/fileio.h>
#include <srvf/matrix.h>

#include <vector>
#include <fstream>

BOOST_AUTO_TEST_SUITE(fileio_tests)

BOOST_AUTO_TEST_CASE(load_csv_test1)
{
  double exp_data[2][4]=
  {
    {0.0, 1.0, 2.0, 3.0},
    {4.0, 5.0, 6.0, 7.0}
  };
  std::ifstream is("data/load_csv_test1.csv");
  std::vector<srvf::Matrix> mats = srvf::io::load_csv(is);
  BOOST_REQUIRE_EQUAL(mats.size(),1);
  BOOST_REQUIRE_EQUAL(mats[0].rows(),2);
  BOOST_REQUIRE_EQUAL(mats[0].cols(),4);
  
  for (int i=0; i<2; ++i)
  {
    for (int j=0; j<4; ++j)
    {
      BOOST_CHECK_CLOSE(mats[0](i,j),exp_data[i][j],1e-4);
    }
  }
}

BOOST_AUTO_TEST_CASE(load_csv_test2)
{
  double exp_data1[2][2]=
  {
    {0.0, 1.0},
    {2.0, 3.0}
  };
  double exp_data2[2][3]=
  {
    {4.0, 5.0, 6.0},
    {7.0, 8.0, 9.0}
  };
  std::ifstream is("data/load_csv_test2.csv");
  std::vector<srvf::Matrix> mats = srvf::io::load_csv(is);
  BOOST_REQUIRE_EQUAL(mats.size(),2);
  BOOST_REQUIRE_EQUAL(mats[0].rows(),2);
  BOOST_REQUIRE_EQUAL(mats[0].cols(),2);
  BOOST_REQUIRE_EQUAL(mats[1].rows(),2);
  BOOST_REQUIRE_EQUAL(mats[1].cols(),3);
  
  for (int i=0; i<2; ++i)
  {
    for (int j=0; j<2; ++j)
    {
      BOOST_CHECK_CLOSE(mats[0](i,j),exp_data1[i][j],1e-4);
    }
  }

  for (int i=0; i<2; ++i)
  {
    for (int j=0; j<3; ++j)
    {
      BOOST_CHECK_CLOSE(mats[1](i,j),exp_data2[i][j],1e-4);
    }
  }
}

BOOST_AUTO_TEST_CASE(save_csv_test1)
{
  std::vector<srvf::Matrix> original_mats;
  double mat1_data[] = 
  {
    0.1, 0.2, 0.3, 0.4,
    0.5, 0.6, 0.7, 0.8
  };
  original_mats.push_back(srvf::Matrix(2,4,mat1_data));

  std::ofstream os("data/save_csv_test1.csv");
  srvf::io::save_csv(os,original_mats);
  os.close();

  std::ifstream is("data/save_csv_test1.csv");
  std::vector<srvf::Matrix> loaded_mats = srvf::io::load_csv(is);
  BOOST_REQUIRE_EQUAL(loaded_mats.size(),original_mats.size());
  for (size_t i=0; i<loaded_mats.size(); ++i)
  {
    BOOST_CHECK(loaded_mats[i].equals( original_mats[i] ));
  }
}

BOOST_AUTO_TEST_SUITE_END()
