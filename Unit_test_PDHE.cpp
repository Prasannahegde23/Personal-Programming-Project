#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include "PDHE_element_routine.h"
#include "PDHE_material_routine.h"
#include "PDHE_boundary_effects.h"
 
using namespace std;
namespace fs = std::filesystem;

int main()
{   
    string outputFileName = "Unit_test_log.txt";
    fs::path currentDir = fs::current_path();
    fs::path outputPath = currentDir / outputFileName;
    std::ofstream outFile(outputPath);

    if(outFile.is_open())
    {
      std::vector<double> vec_1(2);
      std::vector<double> vec_2(2);
    
      vec_1[0] = 1.0; vec_1[1] = 2.0;
      vec_2[0] = 3.0; vec_2[1] = 4.0;

      // Unit test for formation of xi vector
      std::vector<double> result1 = xi_vec(vec_1[0], vec_1[1], vec_2[0], vec_2[1]);
      outFile << "//----------------------------------------------------------//" << endl;
      outFile << "        Unit test for formation of xi vector" << endl;
      outFile << "//----------------------------------------------------------//" << endl;
      outFile << "xi vector: " << endl;
      outFile << "x-component " << result1[0] << endl;
      outFile << "y-component " << result1[1] << endl;
      outFile << "--------------------------------------------------------------------------" << endl << endl;

      // Unit test for magnitude of xi vector
      double result2 = mod_xi(result1);
      outFile << "//----------------------------------------------------------//" << endl;
      outFile << "          Unit test for magnitude of xi vector" << endl;
      outFile << "//----------------------------------------------------------//" << endl;
      outFile << "Magnitude of xi vector: " << result2 << endl;
      outFile << "--------------------------------------------------------------------------" << endl << endl;

      // Unit test for magnitude of eta vector
      std::vector<double> result3 = eta_vec(vec_1, vec_2);
      outFile << "//----------------------------------------------------------//" << endl;
      outFile << "          Unit test for formation of eta vector" << endl;
      outFile << "//----------------------------------------------------------//" << endl;
      outFile << "eta vector: " << endl;
      outFile << "x-component " << result3[0] << endl;
      outFile << "y-component " << result3[1] << endl;
      outFile << "--------------------------------------------------------------------------" << endl << endl;

      
    }
  return 0;
}
