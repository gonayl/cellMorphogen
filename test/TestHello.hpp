#ifndef TESTHELLO_HPP_
#define TESTHELLO_HPP_

#include <cxxtest/TestSuite.h>
#include "FixedBoundaryCondition.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "Hello.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
using namespace std ;

class TestHello : public CxxTest::TestSuite
{
public:
    void TestHelloClass()
    {
      std::cout << "Hello world" << std::endl ;

      // Permet d'importer un fichier test et s'en servir pour taguer les cellules, voir ci-après
      std::cout << "Importing x pos data from txt" << std::endl ;
      ifstream inFilex ;
      double x ;
      std::vector<double> pos_nodes_x;

      inFilex.open("testoutput/test_nodes_x.txt") ;
      if(!inFilex)
      {
        cout << "Unable to open file" << endl ;
      }
      while(inFilex >> x)
      {
        pos_nodes_x.push_back(x) ;
      }
      inFilex.close() ;

      std::cout << "Importing y pos data from txt" << std::endl ;
      ifstream inFiley ;
      double y ;
      std::vector<double> pos_nodes_y;

      inFiley.open("testoutput/test_nodes_y.txt") ;
      if(!inFiley)
      {
        cout << "Unable to open file" << endl ;
      }
      while(inFiley >> y)
      {
        pos_nodes_y.push_back(y) ;
      }
      inFiley.close() ;

    }

};

#endif /*TESTHELLO_HPP_*/
