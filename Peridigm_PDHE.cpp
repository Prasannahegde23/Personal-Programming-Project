#include <iostream>
#include <vector>
#include <cmath>  // or <math.h>
#include <fstream>         // For file I/O
#include <filesystem>      // For current_path()
#include <Teuchos_ParameterList.hpp>
#include <algorithm> // For std::max_element
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
  m_h(0.0), m_min_grid_spacing(0.0),
  m_modelCoordinatesFieldId(-1), m_coordinatesFieldId(-1), m_volumeFieldId(-1), m_concentrationFieldId(-1), m_damageFieldId(-1), m_bodyForceFieldId(-1)

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
    //m_displacementFieldID = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Displacement");
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

std::vector<int> readNodeSet(const std::string& fileName) 
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
}

//---------------------------------------------------------------------------//
// computeForce(): This is called by Peridigm to compute internal forces.
//---------------------------------------------------------------------------//
void PDHE::computeForce(const double dt,
                      const int numOwnedPoints,                                              
                      const int* ownedIDs,
                      const int* neighborhoodList,
                      PeridigmNS::DataManager& dataManager) const
  {
    double *modelCoord; 
    double *currentCoord;
    double *oldCoord;
    double *volume; 
    double *concentration;
    double *damage;
    //double *displacement;
    double *body_force;

    // Access Field IDs with the IDs you grabbed in constructor
    dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&modelCoord);
    dataManager.getData(m_coordinatesFieldId,      PeridigmField::STEP_NP1)->ExtractView(&currentCoord);
    dataManager.getData(m_volumeFieldId,           PeridigmField::STEP_NONE)->ExtractView(&volume);
    dataManager.getData(m_concentrationFieldId,    PeridigmField::STEP_NP1)->ExtractView(&concentration);
    dataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->ExtractView(&damage);
    dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_N)->ExtractView(&oldCoord);
    //dataManager.getData(m_displacementFieldID, PeridigmField::STEP_NP1)->ExtractView(&displacement);


    string outputFileName = "Output.txt";
    fs::path currentDir = fs::current_path();
    fs::path outputPath = currentDir / outputFileName;

    ofstream outFile(outputPath, std::ios::out | std::ios::app);
    //ofstream outFile(outputPath);
    //double converted_Youngs_Modulus = m_Youngs_Modulus * 1000;     //  GPa to MPa
    //double converted_GB_Diff_Coeff = m_GB_Diff_Coeff * 1e6;       // m2/s to mm2/s
    //double m_Sat_Val_Hyd_Conc = m_Sat_Val_Hyd_Conc * 1e-6;
    //double converted_density = m_density * 1e-6;

    double dh = (6 * m_GB_Diff_Coeff)/(M_PI * m_h * (pow(m_horizon,3))); // PD bond constant
    double c = (6 * m_Youngs_Modulus)/(M_PI*m_h*pow(m_horizon,3)*(1 - m_poissons_ratio)); // PD parameter
    double lambda_ii = 1.0;

    // Time steps size decleration
    //double time_step_size_CDM = (sqrt(m_density / (m_Youngs_Modulus))) * 0.7 * (m_min_grid_spacing);
    //double time_step_size_CDM = (sqrt(lambda_ii * 0.0001 / (m_Youngs_Modulus))) * 0.7 * (m_min_grid_spacing);
    double time_step_size_CDM = (sqrt(lambda_ii * 0.0001 / (m_Youngs_Modulus))) * 0.7 * (m_min_grid_spacing);
    //double time_step_size_EFM = (m_min_grid_spacing * m_min_grid_spacing) /(2 * converted_GB_Diff_Coeff);
    double time_step_size_EFM = (sqrt(m_density / (m_Youngs_Modulus))) * (m_min_grid_spacing);
    //double time_step_size_ADR = 1.0;

    
    // Other variables decleration
    std::vector<std::vector<double>> m_bondFactor(numOwnedPoints); 
    //double m_ii = (M_PI * m_horizon * m_horizon * m_h * c)*(time_step_size_ADR * time_step_size_ADR)/(20 * m_min_grid_spacing);

    //std::vector<double> U_dot_n_minus_half(2*numOwnedPoints);
    //std::vector<double> U_dot_n_plus_half(2*numOwnedPoints);

    std::vector<double> old_force(2*numOwnedPoints);
    std::vector<double> old_displacement(2*numOwnedPoints);
    std::vector<double> displacement_n_minus_one(2*numOwnedPoints);

    //std::vector<double> M(3*numOwnedPoints);
    std::vector<double> P(2*numOwnedPoints);
    std::vector<double> K(2*numOwnedPoints);
    //std::vector<double> M_inverse(2*numOwnedPoints);
    //std::vector<double> coordinates(2*numOwnedPoints);
    std::vector<double> displacement(2*numOwnedPoints, 0.0);
    //std::vector<double> concentration(numOwnedPoints, 0.0);
    std::vector<double> old_concentration(numOwnedPoints, 0.0);
    //std::vector<double> new_concentration(numOwnedPoints, 0.0);
    

    double numerator, denominator;
    double k_n = (6*m_Youngs_Modulus)/(M_PI*m_h*pow(m_horizon,3)*(1-m_poissons_ratio));
    double k_t = (6*m_Youngs_Modulus*(1-(3*m_poissons_ratio)))/(M_PI*m_h*pow(m_horizon,3)*(1-m_poissons_ratio));

    std::vector<int> myBoundaryNodes1 = readNodeSet("nodeset_top.txt");
    std::vector<int> myBoundaryNodes2 = readNodeSet("nodeset_bottom.txt");
    std::vector<int> myBoundaryNodes3 = readNodeSet("nodeset_concentration.txt");

    for(int iID=0; iID < numOwnedPoints; ++iID)
    {
      int nodeID = ownedIDs[iID];
      old_concentration[nodeID] = concentration[nodeID];
    }

    std::vector<int> numNeighbors(numOwnedPoints);
    {
      int idx = 0;
      for(int iID=0; iID<numOwnedPoints; ++iID){
        numNeighbors[iID] = neighborhoodList[idx++];
        idx += numNeighbors[iID];
      }
    }

    for(int iID=0; iID<numOwnedPoints; ++iID)
    {m_bondFactor[iID].assign(numNeighbors[iID], 1.0);}

    if(outFile.is_open())
    {
      cout << "PDHE simulation started... "<< endl << endl;
      for(int i=0; i < N_t ; i++)
      {
        cout << "Load steps: "<< i+1 << endl << endl;
        // Displacement BC
        //outFile << "myBoundaryNodes1:" << endl;
        for(auto nodeID : myBoundaryNodes1)
        {
          if(i <= 1000)
          { 
            displacement[2*nodeID] = 0.0;
            displacement[2*nodeID + 1] = (2*1e-8)*i;
          }
          else
          { 
            displacement[2*nodeID] = 0.0;
            displacement[2*nodeID + 1] = (2*1e-5);
          }
          /*displacement[2*nodeID] = 0.0;
          displacement[2*nodeID + 1] = (2*1e-8)*i;*/
          //outFile << "After oldCoord[3*nodeID + 1] " << oldCoord[3*nodeID + 1] << endl;
        }

        for(auto nodeID : myBoundaryNodes2)
        {
          //outFile << "Before oldCoord[3*nodeID + 1] " << oldCoord[3*nodeID + 1] << endl;
          if(i <= 1000)
          {
            displacement[2*nodeID] = 0.0;
            displacement[2*nodeID + 1] = -(2*1e-8)*i; 
          }
          else
          {
            displacement[2*nodeID] = 0.0;
            displacement[2*nodeID + 1] = -(2*1e-5);
          }
          /*displacement[2*nodeID] = 0.0;
          displacement[2*nodeID + 1] = -(2*1e-8)*i;*/
          //outFile << "After oldCoord[3*nodeID + 1] " << oldCoord[3*nodeID + 1] << endl;
        }

        //old_concentration = concentration;

        for(int j=0; j < N_h ; j++)
        {
          //cout << "Hydrogen steps: "<< j+1 << endl;

          for(auto nodeID : myBoundaryNodes3) 
          {old_concentration[nodeID] = m_Sat_Val_Hyd_Conc;}

          //outFile << "Hydrogen_step: " << j << endl;
          int neighIndex = 0; // index into neighborhoodList
          for(int iID=0; iID < numOwnedPoints; ++iID)
          {
            int nodeID = ownedIDs[iID];

            //outFile << displacement[3*nodeID + 1] << endl;
            double x = modelCoord[3*nodeID];
            double y = modelCoord[3*nodeID + 1];
            double Volume_i = volume[nodeID];
            elementroutinehydrogen output;
            //double time_step_size_EFM = (sqrt(converted_density* volume[nodeID]/ (converted_Youngs_Modulus))) * (m_min_grid_spacing);

            int numNeighbors = neighborhoodList[neighIndex++];
            double concentration_nodeID = old_concentration[nodeID];
            output = element_routine_hydrogen(nodeID, modelCoord, x, y, neighborhoodList, neighIndex, numNeighbors, m_horizon, old_concentration, concentration_nodeID, time_step_size_EFM, dh, Volume_i, volume);
            old_concentration[nodeID] = output.conc;
            neighIndex = output.neighindex;
            
          }
          for(auto nodeID : myBoundaryNodes3) 
          {old_concentration[nodeID] = m_Sat_Val_Hyd_Conc;}

          //old_concentration.swap(new_concentration);
          //outFile << endl << endl;
        }

        numerator = 0.0; denominator = 0.0; // Variables used for simplication of calculation
        //outFile << "Iteration: " << iter << endl;
        int neighIndex = 0; // index into neighborhoodList
        for(int iID=0; iID < numOwnedPoints; ++iID)
        {
          int nodeID = ownedIDs[iID];

          //currentCoord[3*nodeID] = currentCoord[3*nodeID] + displacement[3*nodeID];
          //currentCoord[3*nodeID + 1] = currentCoord[3*nodeID + 1] + displacement[3*nodeID + 1];
          //currentCoord[3*nodeID + 2] = currentCoord[3*nodeID + 2] + displacement[3*nodeID + 2];
          //if(i >= 900)
          //{outFile << currentCoord[3*nodeID] << " " << currentCoord[3*nodeID + 1] << " " << displacement[3*nodeID + 1] << endl;} // if required damage value}
          double Px; double Py;
          double x = modelCoord[3*nodeID]; // x
          double y = modelCoord[3*nodeID + 1]; // y
          double Volume_i = volume[nodeID];

          if(old_concentration[nodeID] > 0.0 && damage[nodeID] >= 0.36)
            {old_concentration[nodeID] = m_Sat_Val_Hyd_Conc;}

          int numNeighbors = neighborhoodList[neighIndex++];
          double concenctration_nodeID = old_concentration[nodeID];

          PDResult pdResult = element_routine_PD(Volume_i, volume, c, m_h, m_horizon, k_n, k_t, m_Sat_Val_Hyd_Conc, m_Critic_Energy_Rel_Rate, modelCoord, x, y, nodeID, neighborhoodList, neighIndex, numNeighbors, displacement, old_concentration, concenctration_nodeID, m_min_grid_spacing, m_bondFactor[iID]/*, outFile*/);
          Px = pdResult.Px; Py = pdResult.Py; damage[nodeID] = pdResult.damage; neighIndex = pdResult.neighindex;



          P[2*nodeID] = Px /*+ body_force[nodeID]*/;
          P[2*nodeID + 1] = Py /*+ body_force[nodeID + 1]*/;

          //if(i >= 0)
          //{
           // outFile << "P[2*nodeID]: " << P[2*nodeID] << endl;
           // outFile << "P[2*nodeID + 1]: " << P[2*nodeID + 1] << endl;
          //}
         // outFile << endl;

          // Adavptive dynamic relaxtaion method
         /* M_inverse[2*nodeID] = (1/m_ii); // inverse matrix computation
          M_inverse[2*nodeID + 1] = (1/m_ii);

          if(i == 0)
          {
            U_dot_n_minus_half[2*nodeID] = 0.5 * M_inverse[2*nodeID] * P[2*nodeID]*time_step_size_CDM;
            U_dot_n_minus_half[2*nodeID + 1] = 0.5 * M_inverse[2*nodeID + 1] * P[2*nodeID + 1]*time_step_size_CDM;

            //coordinates[3*nodeID] = modelCoord[3*nodeID];
            //coordinates[3*nodeID + 1] = modelCoord[3*nodeID + 1];
            //coordinates[3*nodeID + 2] = modelCoord[3*nodeID + 2];
          }

          else
          {
            U_dot_n_minus_half[2*nodeID] = U_dot_n_plus_half[2*nodeID];
            U_dot_n_minus_half[2*nodeID + 1] = U_dot_n_plus_half[2*nodeID + 1];
          }*/


          if(i == 0)
            {
              displacement_n_minus_one[2*nodeID] = 0.0;
              displacement_n_minus_one[2*nodeID + 1] = 0.0;
            }

          if(((displacement[2*nodeID] - displacement_n_minus_one[2*nodeID]) == 0.0))
            {K[2*nodeID] = 0.0;}
          else if(((displacement[2*nodeID + 1] - displacement_n_minus_one[2*nodeID + 1]) == 0.0))
            {K[2*nodeID + 1] = 0.0;}
          else
              {
                K[2*nodeID] = -((P[2*nodeID] - old_force[2*nodeID])/lambda_ii)/(displacement[2*nodeID] - displacement_n_minus_one[2*nodeID]);
                K[2*nodeID + 1] = -((P[2*nodeID + 1] - old_force[2*nodeID + 1])/lambda_ii)/(displacement[2*nodeID + 1] - displacement_n_minus_one[2*nodeID + 1]);
              }

          numerator = numerator + (displacement[2*nodeID] * K[2*nodeID] * displacement[2*nodeID])
                                + (displacement[2*nodeID + 1] * K[2*nodeID + 1] * displacement[2*nodeID + 1]);
          denominator = denominator + (displacement[2*nodeID] * displacement[2*nodeID])
                                    + (displacement[2*nodeID + 1] * displacement[2*nodeID + 1]);

        }
        //outFile << "numerator " << numerator << endl;
        //outFile << "denominator " << denominator << endl;
        double c_n = 2 * sqrt(numerator/denominator);
        if(c_n >= 2.0 || numerator < 0.0 || denominator <= 0.0)
          {c_n = 1.9;}
        //outFile << "c_n : " << c_n << endl;


        for(int iID=0; iID < numOwnedPoints; ++iID)
        {
          int nodeID = ownedIDs[iID];
          /*U_dot_n_plus_half[2*nodeID] = (((2 - (time_step_size_CDM * c_n)) * U_dot_n_minus_half[2*nodeID]) +
                           (2*time_step_size_CDM * M_inverse[2*nodeID] * P[2*nodeID]))/(2 + (time_step_size_CDM*c_n));
          U_dot_n_plus_half[2*nodeID + 1] = (((2 - (time_step_size_CDM * c_n)) * U_dot_n_minus_half[3*nodeID + 1]) +
                           (2*time_step_size_CDM * M_inverse[2*nodeID + 1] * P[2*nodeID + 1]))/(2 + (time_step_size_CDM*c_n));

          */

         old_displacement[2*nodeID] = displacement[2*nodeID];
         old_displacement[2*nodeID + 1] = displacement[2*nodeID + 1];

          int flag = 0;
          for(int nodecheckID = 0; nodecheckID < myBoundaryNodes1.size(); nodecheckID++)
          {
            if(nodeID == myBoundaryNodes1[nodecheckID])
            {
              flag = 1;
              break;
            }

            else
              continue;
          }

          for(int nodecheckID = 0; nodecheckID < myBoundaryNodes2.size(); nodecheckID++)
          {
            if(nodeID == myBoundaryNodes2[nodecheckID])
            {
              flag = 1;
              break;
            }

            else
              continue;
          }

          if(flag == 1)
          {
            currentCoord[3*nodeID] = currentCoord[3*nodeID] + displacement[2*nodeID];
            currentCoord[3*nodeID + 1] = currentCoord[3*nodeID + 1] + displacement[2*nodeID + 1];
          }
          
          /*else
          {
            displacement[2*nodeID] = displacement[2*nodeID] + (time_step_size_CDM * U_dot_n_plus_half[2*nodeID]);
            displacement[2*nodeID + 1] = displacement[2*nodeID + 1] + (time_step_size_CDM * U_dot_n_plus_half[2*nodeID + 1]);

            currentCoord[3*nodeID] = currentCoord[3*nodeID] + displacement[2*nodeID];
            currentCoord[3*nodeID + 1] = currentCoord[3*nodeID + 1] + displacement[2*nodeID + 1];
          }*/

          else
          {
            //double time_step_size_CDM = (sqrt(converted_density* volume[nodeID]/ (converted_Youngs_Modulus))) * (m_min_grid_spacing);

            displacement[2*nodeID] = ((2* time_step_size_CDM*time_step_size_CDM * P[2*nodeID]) + (4 * displacement[2*nodeID]) +
                                    (((c_n*time_step_size_CDM) -2) * displacement_n_minus_one[2*nodeID]))/(2 + (time_step_size_CDM*c_n));
            displacement[2*nodeID + 1] = ((2* time_step_size_CDM*time_step_size_CDM * P[2*nodeID + 1]) + (4 * displacement[2*nodeID + 1]) +
                                    (((c_n*time_step_size_CDM) -2) * displacement_n_minus_one[2*nodeID + 1]))/(2 + (time_step_size_CDM*c_n));

            currentCoord[3*nodeID] = currentCoord[3*nodeID] + displacement[2*nodeID];
            currentCoord[3*nodeID + 1] = currentCoord[3*nodeID + 1] + displacement[2*nodeID + 1];
          }
          old_force[2*nodeID] = P[2*nodeID];
          old_force[2*nodeID + 1] = P[2*nodeID + 1];

          displacement_n_minus_one[2*nodeID] = old_displacement[2*nodeID];
          displacement_n_minus_one[2*nodeID + 1] = old_displacement[2*nodeID + 1];

          double x = currentCoord[3*nodeID]; // x
          double y = currentCoord[3*nodeID + 1]; // y
          //double z = coordinates[3*nodeID + 2]; // z
          //displacement[2*nodeID + 1] damage[nodeID]
          if(i >= 2980)
          {
            /*for(auto b : myBoundaryNodes3) 
            {
              if(i <= 1000)
            {
              displacement[2*b + 1] = (2*1e-8)*i; 
            }
            else
            {
              displacement[2*nodeID + 1] = (2*1e-5);
            }
            }*/

            outFile << x << " " << y << " " << displacement[2*nodeID + 1] << " " << old_concentration[nodeID] << " " << damage[nodeID] << endl;} // if required damage value}
          //cout << x << " " << y << " " << z << " " << y/*damage[nodeID]*/ << endl; // if required damage value
        }
        if(i >= 2980)
        {outFile << endl << endl;}
      }
      outFile.close();
      cout << "Data exported to " << outputPath << endl;
    }
  }
}