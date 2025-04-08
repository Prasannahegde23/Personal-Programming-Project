#include <iostream>
#include <vector>
#include <cmath>  // or <math.h>
#include <fstream>         // For file I/O
#include <filesystem>      // For current_path()
#include <Teuchos_ParameterList.hpp>
#include "Peridigm_Field.hpp"
#include "material_utilities.h"
#include "Peridigm_PDHE.hpp"
#include "Peridigm_DataManager.hpp"
#include "PDHE_element_routine.h"
 
using namespace std;
namespace fs = std::filesystem;
namespace PeridigmNS 
{ 

  //---------------------------------------------------------------------------//
  // Constructor: Register which fields your material model needs.
  //---------------------------------------------------------------------------//
  PDHE::PDHE(const Teuchos::ParameterList& params): Material(params),
  m_Youngs_Modulus(0.0), m_GB_Diff_Coeff(0.0), m_Sat_Val_Hyd_Conc(0.0), m_Critic_Energy_Rel_Rate(0.0), m_density(0.0), m_poissons_ratio(0.0),
  m_horizon(0.0), N_t(0), N_h(0), 
  m_h(0.0), m_min_grid_spacing(0.0)

  { 
    //---------------------------------------------------------------------------//
    // In this section, all the required parameters are read from the ParameterList
    //---------------------------------------------------------------------------//

    // Read material parameters
    m_Youngs_Modulus = params.get<double>("Young's Modulus");
    m_GB_Diff_Coeff = params.get<double>("Grain boundary diffusion coefficient");
    m_Sat_Val_Hyd_Conc = params.get<double>("Saturated value of hydrogen concentration");
    m_Critic_Energy_Rel_Rate = params.get<double>("Critical energy release rate");
    m_density = params.get<double>("Density");
    m_poissons_ratio = params.get<double>("Poisson's ratio");

    // Read Input Parameters
    m_horizon = params.get<double>("Horizon");
    N_t = params.get<int>("No. of load steps");
    N_h = params.get<int>("No. of steps for hydrogen concentration");

    // Read Geometrical parameters
    m_h = params.get<double>("Thickness");
    m_min_grid_spacing = params.get<double>("Minimum grid spacing");

    //---------------------------------------------------------------------------//
    // In this section, all the fields are read from the FieldManager where all 
    // data fields are tracked in Peridigm
    //---------------------------------------------------------------------------//

    PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
    m_modelCoordinatesFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
    m_coordinatesFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
    m_volumeFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
    m_concentrationFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Temperature");
    m_damageFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Damage");
    //m_displacementFieldID = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Prescribed_Displacement");
    m_bodyForceFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Body_Force");


    //---------------------------------------------------------------------------//
    // In this section, all the fields are registered so Peridigm knows that 
    // we use these data fields.
    //---------------------------------------------------------------------------//

    m_fieldIds.push_back(m_modelCoordinatesFieldId);
    m_fieldIds.push_back(m_coordinatesFieldId);
    m_fieldIds.push_back(m_volumeFieldId);
    m_fieldIds.push_back(m_concentrationFieldId);
    m_fieldIds.push_back(m_damageFieldId);
    //m_fieldIds.push_back(m_displacementFieldID);
    m_fieldIds.push_back(m_bodyForceFieldId);
  }

//---------------------------------------------------------------------------//
// Implementation of pure virtuals:
//---------------------------------------------------------------------------//

std::string PDHE::Name() const 
{return "PDHE";}                  // Returns anything that identifies your material in logs

double PDHE::YoungsModulus() const 
{return m_Youngs_Modulus;}

double PDHE::GBdiffCoeff() const 
{return m_GB_Diff_Coeff;}

double PDHE::SatValHydConc() const 
{return m_Sat_Val_Hyd_Conc;}

double PDHE::CriticEnergyRelRate() const 
{return m_Critic_Energy_Rel_Rate;}

double PDHE::Density() const 
{return m_density;}

double PDHE::PoissonsRatio() const 
{return m_poissons_ratio;}

double PDHE::Horizon() const 
{return m_horizon;}

int PDHE::Nt() const 
{return N_t;}

int PDHE::Nh() const 
{return N_h;}

double PDHE::h() const 
{return m_h;}

double PDHE::minGridSpacing() const 
{return m_min_grid_spacing;}

std::vector<int> PDHE::FieldIds() const 
{return m_fieldIds;} // Return the vector of field IDs that is used

// Added below to fulfill the pure virtual BulkModulus() and ShearModulus() from Material.hpp:  
double PDHE::BulkModulus() const {
  // Return a value or formula; here, just return 0.0 for demonstration.
  return 0.0;
}

double PDHE::ShearModulus() const {
  // Return a value or formula; here, just return 0.0 for demonstration.
  return 0.0;
}

/*std::vector<int> readNodeSet(const std::string& fileName) 
{
  std::vector<int> nodeSet;
  std::ifstream file(fileName.c_str());
  if (!file.is_open())
  {
    std::cerr << "Error: Could not open node set file: " << fileName << std::endl;
    return nodeSet;
  }
  
  std::string line;
  while (std::getline(file, line)) 
  {
    std::istringstream iss(line);
    int nodeID;
    if (iss >> nodeID) {
      nodeSet.push_back(nodeID);
    }
  }
  return nodeSet;
}*/

//---------------------------------------------------------------------------//
// computeForce(): This is called by Peridigm to compute internal forces.
//---------------------------------------------------------------------------//
void PDHE::computeForce(const double dt,
                      const int numOwnedPoints,                                              
                      const int* ownedIDs,
                      const int* neighborhoodList,
                      PeridigmNS::DataManager& dataManager) const
  {
    double *modelCoord = nullptr; 
    double *currentCoord = nullptr;
    double *oldCoord = nullptr;
    double *volume = nullptr; 
    double *concentration = nullptr;
    double *damage = nullptr;
    //double *displacement = nullptr;
    double *body_force = nullptr;

    // Access Field IDs with the IDs you grabbed in constructor
    dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&modelCoord);
    dataManager.getData(m_coordinatesFieldId,      PeridigmField::STEP_NP1)->ExtractView(&currentCoord);
    dataManager.getData(m_volumeFieldId,           PeridigmField::STEP_NONE)->ExtractView(&volume);
    dataManager.getData(m_concentrationFieldId,    PeridigmField::STEP_NP1)->ExtractView(&concentration);
    dataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->ExtractView(&damage);
    dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_N)->ExtractView(&oldCoord);
    //dataManager.getData(m_displacementFieldID, PeridigmField::STEP_NP1)->ExtractView(&displacement);


    string outputFileName = "output.txt";
    fs::path currentDir = fs::current_path();
    fs::path outputPath = currentDir / outputFileName;

    ofstream outFile(outputPath, std::ios::out | std::ios::app);
    //ofstream outFile(outputPath);
    /*
    cout << "Youngs Modulus: " << m_Youngs_Modulus << endl ;
    cout << "Grain boundary diffusion coefficient: " << m_GB_Diff_Coeff << endl;
    cout << "Saturated value of hydrogen concentration: " << m_Sat_Val_Hyd_Conc << endl;
    cout << "Critical energy release rate: " << m_Critic_Energy_Rel_Rate << endl;
    cout << "Density: " << m_density << endl;
    cout << "Horizon: " << m_horizon << endl;
    cout << "No. of load steps: " << N_t << endl;
    cout << "No. of steps for hydrogen concentration: " << N_h << endl;
    cout << "Thickness: " << m_h << endl;
    cout << "Minimum grid spacing: " << m_min_grid_spacing << endl << endl; */

    double dh = (6 * m_GB_Diff_Coeff)/(M_PI * m_h * (pow(m_horizon,3))); // PD bond constant
    double c = (6 * m_Youngs_Modulus)/(M_PI*pow(m_horizon,4)*(1 - (2*m_poissons_ratio))); // PD parameter

    // Time steps size decleration
    double time_step_size_CDM = (sqrt(m_density * m_Youngs_Modulus)) * 0.7 * m_min_grid_spacing;
    double time_step_size_EFM = (sqrt(m_density * m_Youngs_Modulus)) * 0.5 * m_min_grid_spacing;
    
    // Other variables decleration
    double m_ii = (M_PI * m_horizon * m_horizon * m_h * c)*(time_step_size_CDM * time_step_size_CDM)/(20 * m_min_grid_spacing);
    std::vector<double> U_dot_half(3*numOwnedPoints);
    std::vector<double> U_dot_n_plus_half(3*numOwnedPoints);
    std::vector<double> M(3*numOwnedPoints);
    std::vector<double> P(3*numOwnedPoints);
    std::vector<double> M_inverse(3*numOwnedPoints);
    std::vector<double> K(3*numOwnedPoints);
    std::vector<double> displacement(3*numOwnedPoints);
    double numerator, denominator;
    double k_n = (6*m_Youngs_Modulus)/(M_PI*m_h*pow(m_horizon,3)*(1-m_poissons_ratio));
    double k_t = (6*m_Youngs_Modulus*(1-m_poissons_ratio))/(M_PI*m_h*pow(m_horizon,3)*(1-m_poissons_ratio));

   // std::vector<int> myBoundaryNodes = readNodeSet("nodeset_top.txt");

    if(outFile.is_open())
    {
      cout << "PDHE simulation started... "<< endl << endl;
      /*for(int i=0; i < N_t ; i++)
      {
        cout << "Load steps: "<< i+1 << endl << endl;
        // Displacement BC
        for(auto nodeID : myBoundaryNodes)
        {
          if(i <= 1000)
          {currentCoord[3*nodeID + 1] = currentCoord[3*nodeID + 1] + (0.02/1000)*i;}
          else
          {currentCoord[3*nodeID + 1] = currentCoord[3*nodeID + 1] + 0.02;}
        }*/

        for(int j=0; j < N_h ; j++)
        {
          cout << "Hydrogen steps: "<< j+1 << endl;
          int neighIndex = 0; // index into neighborhoodList
          for(int iID=0; iID < numOwnedPoints; ++iID)
          {
            int nodeID = ownedIDs[iID];

            currentCoord[3*nodeID] = oldCoord[3*nodeID];
            currentCoord[3*nodeID + 1] = oldCoord[3*nodeID + 1];
            currentCoord[3*nodeID + 2] = oldCoord[3*nodeID + 2];

            double x = currentCoord[3*nodeID]; // x
            double y = currentCoord[3*nodeID + 1]; // y
            double z = currentCoord[3*nodeID + 2]; // z
            double Volume_i = volume[nodeID];
            elementroutinehydrogen output;

            displacement[3*nodeID] = x - modelCoord[3*nodeID];
            displacement[3*nodeID + 1] = y - modelCoord[3*nodeID + 1];
            displacement[3*nodeID + 2] = z - modelCoord[3*nodeID + 2];

            int numNeighbors = neighborhoodList[neighIndex++];
            double concentration_nodeID = concentration[nodeID];
            output = element_routine_hydrogen(nodeID, currentCoord, x, y, z, neighborhoodList, neighIndex, numNeighbors, m_horizon, concentration, concentration_nodeID, time_step_size_EFM, dh, Volume_i, volume);
            concentration[nodeID] = output.conc;
            neighIndex = output.neighindex;
          }
        }
        numerator = 0.0; denominator = 0.0; // Variables used for simplication of calculation
        
        int neighIndex = 0; // index into neighborhoodList
        for(int iID=0; iID < numOwnedPoints; ++iID)
        {
          int nodeID = ownedIDs[iID];

          currentCoord[3*nodeID] = oldCoord[3*nodeID];
          currentCoord[3*nodeID + 1] = oldCoord[3*nodeID + 1];
          currentCoord[3*nodeID + 2] = oldCoord[3*nodeID + 2];
          
          double Px; double Py; double Pz;
          double x = currentCoord[3*nodeID]; // x
          double y = currentCoord[3*nodeID + 1]; // y
          double z = currentCoord[3*nodeID + 2]; // z
          double Volume_i = volume[nodeID];
          
          int numNeighbors = neighborhoodList[neighIndex++];
          double concenctration_nodeID = concentration[nodeID];

          PDResult pdResult = element_routine_PD(Volume_i, volume, c, m_horizon, k_n, k_t, m_Sat_Val_Hyd_Conc, m_Critic_Energy_Rel_Rate, currentCoord, x, y, z, nodeID, neighborhoodList, neighIndex, numNeighbors, displacement, concentration, concenctration_nodeID, m_min_grid_spacing);
          Px = pdResult.Px; Py = pdResult.Py; Pz = pdResult.Pz; damage[nodeID] = pdResult.damage; neighIndex = pdResult.neighindex;
          P[3*nodeID] = Px /*+ body_force[nodeID]*/;
          P[3*nodeID + 1] = Py /*+ body_force[nodeID + 1]*/;
          P[3*nodeID + 2] = Pz /*+ body_force[nodeID + 2]*/;

          // Adavptive dynamic relaxtaion method
          M_inverse[3*nodeID] = (1/m_ii); // inverse matrix computation
          M_inverse[3*nodeID + 1] = (1/m_ii);
          M_inverse[3*nodeID + 2] = (1/m_ii);

          U_dot_half[3*nodeID] = 0.5 * M_inverse[3*nodeID] * P[3*nodeID]*time_step_size_CDM;
          U_dot_half[3*nodeID + 1] = 0.5 * M_inverse[3*nodeID + 1] * P[3*nodeID + 1]*time_step_size_CDM;
          U_dot_half[3*nodeID + 2] = 0.5 * M_inverse[3*nodeID + 2] * P[3*nodeID + 2]*time_step_size_CDM;

          if(nodeID == 0) // At initial material point F_n-1 = 1
            {
              if(U_dot_half[3*nodeID] == 0.0)
                {K[3*nodeID] = 0;}
              else if(U_dot_half[3*nodeID + 1] == 0.0)
                {K[3*nodeID + 1] = 0;}
              else if(U_dot_half[3*nodeID + 2] == 0.0)
                {K[3*nodeID + 2] = 0;}
              else
              {
                K[3*nodeID] = -((P[3*nodeID]/m_ii) - (1/m_ii))/(time_step_size_CDM * U_dot_half[3*nodeID]);
                K[3*nodeID + 1] = -((P[3*nodeID + 1]/m_ii) - (1/m_ii))/(time_step_size_CDM * U_dot_half[3*nodeID + 1]);
                K[3*nodeID + 2] = -((P[3*nodeID + 2]/m_ii) - (1/m_ii))/(time_step_size_CDM * U_dot_half[3*nodeID + 2]);
              }
            }

          else 
            {
              if(U_dot_half[3*nodeID] == 0.0)
                {K[3*nodeID] = 0;}
              else if(U_dot_half[3*nodeID + 1] == 0.0)
                {K[3*nodeID + 1] = 0;}
              else if(U_dot_half[3*nodeID + 2] == 0.0)
                {K[3*nodeID + 2] = 0;}
              else
              {
                K[3*nodeID] = -((P[3*nodeID]/m_ii) - (P[3*nodeID-1]/m_ii))/(time_step_size_CDM * U_dot_half[3*nodeID]);
                K[3*nodeID + 1] = -((P[3*nodeID + 1]/m_ii) - (P[3*nodeID + 1-1]/m_ii))/(time_step_size_CDM * U_dot_half[3*nodeID + 1]);
                K[3*nodeID + 2] = -((P[3*nodeID + 2]/m_ii) - (P[3*nodeID + 2-1]/m_ii))/(time_step_size_CDM * U_dot_half[3*nodeID + 2]);
              }
            }
          numerator = numerator + (displacement[3*nodeID] * K[3*nodeID] * displacement[3*nodeID])
                                + (displacement[3*nodeID + 1] * K[3*nodeID + 1] * displacement[3*nodeID + 1])
                                + (displacement[3*nodeID + 2] * K[3*nodeID + 2] * displacement[3*nodeID + 2]);
          denominator = denominator + (displacement[3*nodeID] * displacement[3*nodeID])
                                    + (displacement[3*nodeID + 1] * displacement[3*nodeID + 1])
                                    + (displacement[3*nodeID + 2] * displacement[3*nodeID + 2]);

        }

        double c_n;
        c_n = 2 * sqrt(numerator/denominator);
        if(numerator <= 0.0 || c_n > 2.0)
        {c_n = 1.9;}

        for(int iID=0; iID < numOwnedPoints; ++iID)
        {
          int nodeID = ownedIDs[iID];
          U_dot_n_plus_half[3*nodeID] = (((2 - (time_step_size_CDM * c_n)) * U_dot_half[3*nodeID]) +
                           (2*time_step_size_CDM * M_inverse[3*nodeID] * P[3*nodeID]))/(2 + (time_step_size_CDM*c_n));
          U_dot_n_plus_half[3*nodeID + 1] = (((2 - (time_step_size_CDM * c_n)) * U_dot_half[3*nodeID + 1]) +
                           (2*time_step_size_CDM * M_inverse[3*nodeID + 1] * P[3*nodeID + 1]))/(2 + (time_step_size_CDM*c_n));
          U_dot_n_plus_half[3*nodeID + 2] = (((2 - (time_step_size_CDM * c_n)) * U_dot_half[3*nodeID + 2]) +
                           (2*time_step_size_CDM * M_inverse[3*nodeID + 2] * P[3*nodeID + 2]))/(2 + (time_step_size_CDM*c_n));

          displacement[3*nodeID] = displacement[3*nodeID] + (time_step_size_CDM * U_dot_n_plus_half[3*nodeID]);
          displacement[3*nodeID + 1] = displacement[3*nodeID + 1] + (time_step_size_CDM * U_dot_n_plus_half[3*nodeID + 1]);
          displacement[3*nodeID + 2] = displacement[3*nodeID + 2] + (time_step_size_CDM * U_dot_n_plus_half[3*nodeID + 2]);

          currentCoord[3*nodeID] = currentCoord[3*nodeID] + displacement[3*nodeID];
          currentCoord[3*nodeID + 1] = currentCoord[3*nodeID + 1] + displacement[3*nodeID + 1];
          currentCoord[3*nodeID + 2] = currentCoord[3*nodeID + 2] + displacement[3*nodeID + 2];

          oldCoord[3*nodeID] = currentCoord[3*nodeID];
          oldCoord[3*nodeID + 1] = currentCoord[3*nodeID + 1];
          oldCoord[3*nodeID + 2] = currentCoord[3*nodeID + 2];

          double x = currentCoord[3*nodeID]; // x
          double y = currentCoord[3*nodeID + 1]; // y
          outFile << x << " " << y << " " << displacement[3*nodeID + 1] << endl; // if required damage value
          //cout << x << " " << y << " " << z << " " << y/*damage[nodeID]*/ << endl; // if required damage value
        }
        outFile << endl << endl;
      }
      outFile.close();
      cout << "Data exported to " << outputPath << endl;  
    }
  }
//}