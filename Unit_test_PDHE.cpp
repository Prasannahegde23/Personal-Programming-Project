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
    ofstream outFile(outputPath, std::ios::out | std::ios::app);

    if(outFile.is_open())
    {
      std::vector<double> vec_1(2);
      std::vector<double> vec_2(2);
    
      vec_1[0] = 1.0; vec_1[1] = 2.0;
      vec_2[0] = 1.0; vec_2[1] = 2.0;

      double result = mod_xi(vec_1);
      outFile << "Magnitude of xi vector: " << result << endl;
      
    }
  return 0;
}
