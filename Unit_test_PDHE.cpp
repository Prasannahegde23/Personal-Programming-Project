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
      outFile << "x-component = " << result1[0] << endl;
      outFile << "y-component = " << result1[1] << endl;
      outFile << "--------------------------------------------------------------------------" << endl << endl << endl;
      
      // Unit test for magnitude of xi vector
      double result2 = mod_xi(result1);
      outFile << "//----------------------------------------------------------//" << endl;
      outFile << "          Unit test for magnitude of xi vector" << endl;
      outFile << "//----------------------------------------------------------//" << endl;
      outFile << "Magnitude of xi vector: " << result2 << endl;
      outFile << "--------------------------------------------------------------------------" << endl << endl << endl;

      // Unit test for magnitude of eta vector
      std::vector<double> result3 = eta_vec(vec_1, vec_2);
      outFile << "//----------------------------------------------------------//" << endl;
      outFile << "          Unit test for formation of eta vector" << endl;
      outFile << "//----------------------------------------------------------//" << endl;
      outFile << "eta vector: " << endl;
      outFile << "x-component = " << result3[0] << endl;
      outFile << "y-component = " << result3[1] << endl;
      outFile << "--------------------------------------------------------------------------" << endl << endl << endl;

      // Unit test for surface correction factor
      double k_n = 1.0; double k_t = 2.0; double m_horizon = 3.0; double V_i = 4.0; double V_j = 5.0; double m_h = 6.0;
      PDParameter result4 = surface_correction(k_n, k_t, m_horizon, V_i, V_j, m_h);
      outFile << "//----------------------------------------------------------//" << endl;
      outFile << "          Unit test for surface correction factor" << endl;
      outFile << "//----------------------------------------------------------//" << endl;
      outFile << "PD parameter before correction: " << endl;
      outFile << "k_n = " << k_n << endl;
      outFile << "k_t = " << k_t << endl << endl;

      outFile << "PD parameter after correction: " << endl;
      outFile << "k_n = " << result4.K_n << endl;
      outFile << "k_t = " << result4.K_t << endl;
      outFile << "--------------------------------------------------------------------------" << endl << endl << endl;

      // Unit test for volume correction factor
      double m_min_grid_spacing = 4.0;
      outFile << "//----------------------------------------------------------//" << endl;
      outFile << "          Unit test for volume correction factor" << endl;
      outFile << "//----------------------------------------------------------//" << endl;
      outFile << "PD-force vector before correction: " << endl;
      outFile << "x-component = " << vec_1[0] << endl;
      outFile << "y-component = " << vec_1[1] << endl << endl;

      std::vector<double> result5 = volume_correction(vec_1, m_horizon, result2, m_min_grid_spacing);

      outFile << "PD-force vector after correction: " << endl;
      outFile << "x-component = " << result5[0] << endl;
      outFile << "y-component = " << result5[1] << endl;
      outFile << "--------------------------------------------------------------------------" << endl << endl << endl;

      // Unit test for magnitude of xi vector
      vec_1[0] = 1.0; vec_1[1] = 2.0;
      double result6 = vector_magnitude(vec_1);
      outFile << "//----------------------------------------------------------//" << endl;
      outFile << "          Unit test for magnitude of a vector" << endl;
      outFile << "//----------------------------------------------------------//" << endl;
      outFile << "Magnitude of xi vector: " << result6 << endl;
      outFile << "--------------------------------------------------------------------------" << endl << endl << endl;
      
      // Unit test for magnitude of eta vector
      std::vector<double> result7 = addVectors(vec_1, vec_2);
      outFile << "//----------------------------------------------------------//" << endl;
      outFile << "          Unit test for addition of two vectors" << endl;
      outFile << "//----------------------------------------------------------//" << endl;
      outFile << "Vector_1: " << "[" << vec_1[0] << " " << vec_1[1] << "]" << endl;
      outFile << "Vector_2: " << "[" << vec_2[0] << " " << vec_2[1] << "]" << endl << endl;

      outFile << "Addition of Vector_1 and Vector_2: " << endl;
      outFile << "x-component = " << result7[0] << endl;
      outFile << "y-component = " << result7[1] << endl;
      outFile << "--------------------------------------------------------------------------" << endl << endl << endl;

    }
  return 0;
}
