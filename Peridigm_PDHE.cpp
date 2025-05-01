#include <iostream>
#include <vector>
#include <cmath> 
#include <fstream> 
#include <filesystem> 
#include <Teuchos_ParameterList.hpp>
#include <algorithm> 
#include <unordered_set>
#include <numeric> 
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
  // Constructor: Register required fields and read parameters from input
  //---------------------------------------------------------------------------//
  PDHE::PDHE(const Teuchos::ParameterList& params): 
  Material(params), // Base class constructor

  // Initialize all scalar members to zero or invalid ID
  m_Youngs_Modulus(0.0), 
  m_GB_Diff_Coeff(0.0), 
  m_Sat_Val_Hyd_Conc(0.0), 
  m_Critic_Energy_Rel_Rate(0.0), 
  m_density(0.0), 
  m_poissons_ratio(0.0),
  m_horizon(0.0), 
  N_t(0), N_h(0), N(0),
  m_h(0.0), m_min_grid_spacing(0.0),
  m_modelCoordinatesFieldId(-1), 
  m_coordinatesFieldId(-1), 
  m_volumeFieldId(-1), 
  m_concentrationFieldId(-1), 
  m_damageFieldId(-1), 
  m_bodyForceFieldId(-1)

  { 
    //---------------------------------------------------------------------------//
    // In this section, all the required parameters are read from the ParameterList
    //---------------------------------------------------------------------------//

    // Read material parameters from the user-specified ParameterList
    m_Youngs_Modulus = params.get<double>("Young's Modulus");
    m_GB_Diff_Coeff = params.get<double>("Grain boundary diffusion coefficient");
    m_Sat_Val_Hyd_Conc = params.get<double>("Saturated value of hydrogen concentration");
    m_Critic_Energy_Rel_Rate = params.get<double>("Critical energy release rate");
    m_density = params.get<double>("Density");
    m_poissons_ratio = params.get<double>("Poisson's ratio");

    // Read simulation control parameters
    m_horizon = params.get<double>("Horizon");
    N_t = params.get<int>("No. of load steps");
    N_h = params.get<int>("No. of steps for hydrogen concentration");
    N = params.get<int>("Capture and save the simulation frame from N load steps");

    // Read geometric parameters
    m_h = params.get<double>("Thickness");
    m_min_grid_spacing = params.get<double>("Minimum grid spacing");

    //Read File names
    m_top_displacement = params.get<std::string>("File name for displacement in +ve x/y-direction");
    m_bottom_displacement = params.get<std::string>("File name for displacement in -ve x/y-direction");
    m_crack_top = params.get<std::string>("File name for crack top lip");
    m_crack_bottom = params.get<std::string>("File name for crack bottom lip");
    m_concentration_file = params.get<std::string>("File name for concentration");

    BC_status = params.get<bool>("Boundary condition test"); // Used to capture status of boundary condition test

    //---------------------------------------------------------------------------//
    // In this section, all the fields are read from the FieldManager where all 
    // data fields are tracked in Peridigm
    //---------------------------------------------------------------------------//

    // Obtain field IDs from Peridigm's FieldManager
    PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
    m_modelCoordinatesFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
    m_coordinatesFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
    m_volumeFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
    m_concentrationFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Temperature");
    m_damageFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Damage");
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
    m_fieldIds.push_back(m_bodyForceFieldId);
  }

//---------------------------------------------------------------------------//
// Material metadata: unique name and property accessors
//---------------------------------------------------------------------------//

std::string PDHE::Name() const 
{return "PDHE";}

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

int PDHE::Captureloadsteps() const 
{return N;}

double PDHE::h() const 
{return m_h;}

double PDHE::minGridSpacing() const 
{return m_min_grid_spacing;}

std::vector<int> PDHE::FieldIds() const 
{return m_fieldIds;}

double PDHE::BulkModulus() const 
{return 0.0;}

double PDHE::ShearModulus() const 
{return 0.0;}

bool PDHE::BCtest() const 
{return BC_status;}

//---------------------------------------------------------------------------//
// readNodeSet: Utility to load a list of node IDs from a text file       //
//---------------------------------------------------------------------------//

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
    if (iss >> nodeID) 
      {nodeSet.push_back(nodeID);}
  }
  return nodeSet;
}

//---------------------------------------------------------------------------//
// UnionFind: Simple disjoint-set data structure for connectivity tests     //
//---------------------------------------------------------------------------//

struct UnionFind 
{
  std::vector<int> parent;
  // Initialize each node to be its own parent
  UnionFind(int N) : parent(N) { std::iota(parent.begin(), parent.end(), 0); }
   // Path compression
  int find(int i){ return parent[i]==i ? i : parent[i]=find(parent[i]); }
  void unite(int i,int j)
  {
    i = find(i); j = find(j);
    if(i!=j) parent[j]=i;
  }
};

//---------------------------------------------------------------------------//
// computeForce(): This is called by Peridigm
//---------------------------------------------------------------------------//
void PDHE::computeForce(const double dt,
                      const int numOwnedPoints,                                              
                      const int* ownedIDs,
                      const int* neighborhoodList,
                      PeridigmNS::DataManager& dataManager) const
  {
    // Extract pointers to all required field data
    double *modelCoord; 
    double *currentCoord;
    double *oldCoord;
    double *volume; 
    double *concentration;
    double *body_force;

    // Access Field IDs with the IDs that were grabbed in constructor
    dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&modelCoord);
    dataManager.getData(m_coordinatesFieldId,      PeridigmField::STEP_NP1)->ExtractView(&currentCoord);
    dataManager.getData(m_volumeFieldId,           PeridigmField::STEP_NONE)->ExtractView(&volume);
    dataManager.getData(m_concentrationFieldId,    PeridigmField::STEP_NP1)->ExtractView(&concentration);
    dataManager.getData(m_coordinatesFieldId,      PeridigmField::STEP_N)->ExtractView(&oldCoord);

    // Prepare output file for logging results
    string outputFileName = "Output.txt";
    fs::path currentDir = fs::current_path();
    fs::path outputPath = currentDir / outputFileName;

    ofstream outFile(outputPath, std::ios::out | std::ios::app);

    // Compute peridynamic bond constants
    double dh = (6 * m_GB_Diff_Coeff)/(M_PI * m_h * (pow(m_horizon,3))); // PD bond constant
    double c = (6 * m_Youngs_Modulus)/(M_PI*m_h*pow(m_horizon,3)*(1 - m_poissons_ratio)); // PD parameter
    double lambda_ii = 1.0; // lambda scaling

    // Estimate stable time step sizes
    double time_step_size_CDM = (sqrt(m_density * 0.0001 / (m_Youngs_Modulus))) * (m_min_grid_spacing);
    double time_step_size_EFM = (sqrt(m_density / (m_Youngs_Modulus))) * (m_min_grid_spacing);

    
    // Other variables decleration
    std::vector<double> old_force(2*numOwnedPoints);
    std::vector<double> old_displacement(2*numOwnedPoints);
    std::vector<double> displacement_n_minus_one(2*numOwnedPoints);
    std::vector<double> displacement(2*numOwnedPoints, 0.0);
    std::vector<double> old_concentration(numOwnedPoints, 0.0);
    std::vector<double> P(2*numOwnedPoints);
    std::vector<double> Reac_Force_bnd_nodes(2*numOwnedPoints, 0.0);
    std::vector<double> K(2*numOwnedPoints);
    double numerator, denominator;

    //Precompute PD parameters
    double k_n = (6*m_Youngs_Modulus)/(M_PI*m_h*pow(m_horizon,3)*(1-m_poissons_ratio));
    double k_t = (6*m_Youngs_Modulus*(1-(3*m_poissons_ratio)))/(M_PI*m_h*pow(m_horizon,3)*(1-m_poissons_ratio));

    // Load boundary node sets: top displacement, bottom displacement, hydrogen BCs, crack nodeIds
    
    std::vector<int> myBoundaryNodes1 = readNodeSet(m_top_displacement);
    std::vector<int> myBoundaryNodes2 = readNodeSet(m_bottom_displacement);
    std::vector<int> myBoundaryNodes3 = readNodeSet(m_concentration_file);
    std::vector<int> crackTop = readNodeSet(m_crack_top);
    std::vector<int> crackBottom = readNodeSet(m_crack_bottom);


    // Initialize old_concentration from the field
    for(int iID=0; iID < numOwnedPoints; ++iID)
    {
      int nodeID = ownedIDs[iID];
      old_concentration[nodeID] = concentration[nodeID];
    }

    // Build global-to-local ID map for quick neighbor lookups
    std::unordered_map<int,int> globalToLocal;
    globalToLocal.reserve(numOwnedPoints);
    for(int iID = 0; iID < numOwnedPoints; ++iID)
      {globalToLocal[ ownedIDs[iID] ] = iID;}

    // Read neighbor counts per node to reconstruct adjacency lists, will be used for storing bond factor
    std::vector<int> numNeighbors(numOwnedPoints);
    {
      int idx = 0;
      for(int iID=0; iID<numOwnedPoints; ++iID)
      {
        numNeighbors[iID] = neighborhoodList[idx++];
        idx += numNeighbors[iID];
      }
    }

    // Initialize bond factors (1.0 = intact, 0.0 = broken)
    m_bondFactor.resize(numOwnedPoints);
    damage.assign(numOwnedPoints, 0.0);
    for(int iID=0; iID<numOwnedPoints; ++iID)
    {m_bondFactor[iID].assign(numNeighbors[iID], 1.0);}

    // Convert the raw neighbor list into a vector-of-vectors for easier access
    vector<vector<int>> neighborListVec(numOwnedPoints);
    int idx=0;
    for(int i=0;i<numOwnedPoints;++i)
    {
      int n = neighborhoodList[idx++];
      neighborListVec[i].assign(&neighborhoodList[idx], &neighborhoodList[idx+n]);
      idx += n;
    }

    //Sever initial crack bonds symmetrically, avoiding boundary nodes
    std::unordered_set<int> crackTopSet(crackTop.begin(), crackTop.end());
    std::unordered_set<int> crackBotSet(crackBottom.begin(), crackBottom.end());
    std::unordered_set<int> displacementBC;
    for(auto id : myBoundaryNodes1) 
      {displacementBC.insert(id);}
    for(auto id : myBoundaryNodes2) 
      {displacementBC.insert(id);}

// break all bonds along the initial crack face, _symmetrically_
  for(int i=0; i<numOwnedPoints; ++i)
  {
    int g = ownedIDs[i];
    // Not severing bonds attached to displacement BC nodes
    if(displacementBC.count(g))
      continue;

    auto &nbrs  = neighborListVec[i];
    auto &bonds = m_bondFactor[i];
    for(int n=0; n<(int)nbrs.size(); ++n)
    {
      int h = nbrs[n];
    // Skip if neighbor is a BC node
      if(displacementBC.count(h))
        continue;

      bool top2bot = crackTopSet.count(g) && crackBotSet.count(h);
      bool bot2top = crackBotSet.count(g) && crackTopSet.count(h);
      if(top2bot || bot2top)
      {
        // Sever bond on both sides
        bonds[n] = 0.0;

        // mirror-sever on the other side
        auto it = globalToLocal.find(h);
        if(it != globalToLocal.end()){
          int j = it->second;
          auto &nbrsJ  = neighborListVec[j];
          auto &bondsJ = m_bondFactor[j];
          for(int m=0; m<(int)nbrsJ.size(); ++m)
          {
            if(nbrsJ[m] == g)
            {
              bondsJ[m] = 0.0;
              break;
            }
          }
        }
      }
    }
  }

  // Begin main load-step loop for mechanical and diffusion updates
  if(outFile.is_open())
  { 
    // Header for output file
    if(BC_status == true)
    {
      outFile << "#Header"<< endl;
      outFile << "nodeID value"<< endl << endl << endl;
    }

    else
    {
      outFile << "#Header"<< endl;
      outFile << "x y Displacement(m) Hydrogen_concentration(mol/m^2) Damage"<< endl << endl;
    }


    cout << "PDHE simulation started... "<< endl << endl;
    for(int i=0; i < N_t+1 ; i++)
    {
      cout << "Load steps: "<< i << endl << endl;
        
      // Apply prescribed displacement BCs on top and bottom nodes
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
      }

      for(auto nodeID : myBoundaryNodes2)
      {
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
      }

      //Writing Load step number into text file
      if(i % 100 == 0)
        {outFile << "Load step: " << i << endl;}

      // Hydrogen diffusion subcycling: N_h sub-steps per mechanical step
      for(int j=0; j < N_h ; j++)
      {
        // Enforce saturated concentration at BC nodes
        for(auto nodeID : myBoundaryNodes3) 
        {old_concentration[nodeID] = m_Sat_Val_Hyd_Conc;}

        // Loop over all points for element_routine_hydrogen
        int neighIndex = 0; // index into neighborhoodList
        for(int iID=0; iID < numOwnedPoints; ++iID)
        {
          int nodeID = ownedIDs[iID];
          double x = modelCoord[3*nodeID];
          double y = modelCoord[3*nodeID + 1];
          double Volume_i = volume[nodeID];
          elementroutinehydrogen output;

          // Extract neighbor count and run hydrogen update
          int numNeighbors = neighborhoodList[neighIndex++];
          double concentration_nodeID = old_concentration[nodeID];
          output = element_routine_hydrogen(nodeID, modelCoord, x, y, neighborhoodList, neighIndex, numNeighbors, m_horizon, old_concentration, concentration_nodeID, time_step_size_EFM, dh, Volume_i, volume);
          old_concentration[nodeID] = output.conc;
          neighIndex = output.neighindex;
            
        }

        // Reapply Dirichlet BC for hydrogen in order to have saturated hydrogen concentration at crack nodes
        for(auto b : myBoundaryNodes3) 
          {old_concentration[b] = m_Sat_Val_Hyd_Conc;}
      }

      // Connectivity tracking: use UnionFind to cluster intact bonds
      UnionFind uf(numOwnedPoints);
      int idx = 0;
      for(int iID=0; iID<numOwnedPoints; ++iID)
      {
        int nbors = neighborhoodList[idx++];
        for(int n=0; n<nbors; ++n)
        {
          int nbrGlobal = neighborhoodList[idx + n];
          if(m_bondFactor[iID][n] > 0.0)
          {
            auto it = globalToLocal.find(nbrGlobal);
            if(it != globalToLocal.end())
            uf.unite(iID, it->second);
          }
        }
        idx += nbors;
      }

      // Sever bonds between different UF clusters (crack growth)
      idx = 0;
      for(int iID=0; iID<numOwnedPoints; ++iID)
      {
        int g = ownedIDs[iID];
        int nbors = neighborhoodList[idx++];

        if(displacementBC.count(g))
        {
          idx += nbors;
          continue;
        }

        for(int n=0; n<nbors; ++n)
        {
          int h = neighborhoodList[idx + n];
          if(displacementBC.count(h))
            {continue;}

          auto it = globalToLocal.find(h);
          if(it != globalToLocal.end() && uf.find(iID) != uf.find(it->second))
          {
            // Break bond on both sides
            m_bondFactor[iID][n] = 0.0;
            auto &nbrsJ = neighborListVec[it->second];
            auto &bJ     = m_bondFactor[it->second];
            for(int m=0; m<(int)nbrsJ.size(); ++m)
            {
              if(nbrsJ[m] == g)
              {
                bJ[m] = 0.0;
                break;
              }
            }
          }
        }
        idx += nbors;
      }


    // Compute internal forces, damage, and effective stiffness K
    numerator = 0.0; denominator = 0.0; // Variables used for simplication of calculation
    int neighIndex = 0; // index into neighborhoodList
    for(int iID=0; iID < numOwnedPoints; ++iID)
    {
      int nodeID = ownedIDs[iID];
      double Px; double Py;
      double x = modelCoord[3*nodeID];
      double y = modelCoord[3*nodeID + 1];
      double Volume_i = volume[nodeID];

      // If damage is high, clamp conc to increase hydrogen embrittlement
      if(old_concentration[nodeID] > 0.0 && damage[nodeID] >= 0.36)
        {old_concentration[nodeID] = m_Sat_Val_Hyd_Conc;}

      int numNeighbors = neighborhoodList[neighIndex++];
      double concenctration_nodeID = old_concentration[nodeID];

      // Zero out force and restore bonds for BC rows
      if(displacementBC.count(nodeID)) 
      {
        damage[iID] = 0.0;
        std::fill(m_bondFactor[iID].begin(), m_bondFactor[iID].end(), 1.0);
        P[2*nodeID]     = 0.0;
        P[2*nodeID + 1] = 0.0;
        neighIndex += numNeighbors;
        continue;
      }

      // Call PD element routine to compute forces and updated bonds
      PDResult pdResult = element_routine_PD(Volume_i, volume, c, m_h, m_horizon, k_n, k_t, m_Sat_Val_Hyd_Conc, m_Critic_Energy_Rel_Rate, modelCoord, x, y, nodeID, neighborhoodList, neighIndex, numNeighbors, displacement, old_concentration, concenctration_nodeID, m_min_grid_spacing, m_bondFactor[iID]/*, outFile*/);
      Px = pdResult.Px; Py = pdResult.Py; neighIndex = pdResult.neighindex; m_bondFactor[iID] = pdResult.bondFac;
      P[2*nodeID] = Px /*+ body_force[nodeID]*/;
      P[2*nodeID + 1] = Py /*+ body_force[nodeID + 1]*/;

      // Update damage based on broken bonds
      double sum = std::accumulate(m_bondFactor[iID].begin(), m_bondFactor[iID].end(), 0.0);
      double dNew = 1.0 - sum/m_bondFactor[iID].size();
      damage[iID] = std::max(damage[iID], dNew);


       // Compute local stiffness K for computing nodal displacements
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
      
      // Accumulate for computing damping co-efficient
      numerator = numerator + (displacement[2*nodeID] * K[2*nodeID] * displacement[2*nodeID])
                            + (displacement[2*nodeID + 1] * K[2*nodeID + 1] * displacement[2*nodeID + 1]);
      denominator = denominator + (displacement[2*nodeID] * displacement[2*nodeID])
                                + (displacement[2*nodeID + 1] * displacement[2*nodeID + 1]);

    }

    // Calculate damping co-efficient, force the damping co-efficient if it is met with the conditions
    double c_n = 2 * sqrt(numerator/denominator);
    if(c_n >= 2.0 || numerator < 0.0 || denominator <= 0.0)
      {c_n = 1.9;}

    // Update displacements and coordinates using central-difference scheme
    for(int iID=0; iID < numOwnedPoints; ++iID)
    {
      int nodeID = ownedIDs[iID];

      old_displacement[2*nodeID] = displacement[2*nodeID];
      old_displacement[2*nodeID + 1] = displacement[2*nodeID + 1];

      int flag = 0;
      for(std::size_t nodecheckID = 0; nodecheckID < myBoundaryNodes1.size(); nodecheckID++)
      {
        if(nodeID == myBoundaryNodes1[nodecheckID])
        {
          flag = 1;
          break;
        }

        else
          {continue;}
      }

      for(std::size_t nodecheckID = 0; nodecheckID < myBoundaryNodes2.size(); nodecheckID++)
      {
        if(nodeID == myBoundaryNodes2[nodecheckID])
        {
          flag = 1;
          break;
        }

        else
          {continue;}
      }

      if(flag == 1)
      {
        // If at all wanted for deformed configuration
        currentCoord[3*nodeID] = currentCoord[3*nodeID] + displacement[2*nodeID];
        currentCoord[3*nodeID + 1] = currentCoord[3*nodeID + 1] + displacement[2*nodeID + 1];
      }

      else
      {
        displacement[2*nodeID] = ((2* time_step_size_CDM*time_step_size_CDM * P[2*nodeID]) + (4 * displacement[2*nodeID]) +
                                (((c_n*time_step_size_CDM) - 2) * displacement_n_minus_one[2*nodeID]))/(2 + (time_step_size_CDM*c_n));
        displacement[2*nodeID + 1] = ((2* time_step_size_CDM*time_step_size_CDM * P[2*nodeID + 1]) + (4 * displacement[2*nodeID + 1]) +
                                    (((c_n*time_step_size_CDM) - 2) * displacement_n_minus_one[2*nodeID + 1]))/(2 + (time_step_size_CDM*c_n));

        currentCoord[3*nodeID] = currentCoord[3*nodeID] + displacement[2*nodeID];
        currentCoord[3*nodeID + 1] = currentCoord[3*nodeID + 1] + displacement[2*nodeID + 1];
      }
      
      old_force[2*nodeID] = P[2*nodeID];
      old_force[2*nodeID + 1] = P[2*nodeID + 1];

      displacement_n_minus_one[2*nodeID] = old_displacement[2*nodeID];
      displacement_n_minus_one[2*nodeID + 1] = old_displacement[2*nodeID + 1];

      double x = modelCoord[3*nodeID];
      double y = modelCoord[3*nodeID + 1];
      
      // Log results after capture threshold

      // Will be used for boundary condition test
      if(BC_status == true && i >= N)
      {
        for(std::size_t nodecheckID = 0; nodecheckID < myBoundaryNodes1.size(); nodecheckID++)
        {

          if(nodeID == myBoundaryNodes1[nodecheckID])
          {
            outFile << "Check for Displacement boundary condition in m" << endl;
            outFile << nodeID << " " << displacement[2*nodeID + 1] << endl;
            break;
          }
          else
            {continue;}
        }

        for(std::size_t nodecheckID = 0; nodecheckID < myBoundaryNodes3.size(); nodecheckID++)
        {
          if(nodeID == myBoundaryNodes3[nodecheckID])
          {
            outFile << "Check for Concentration boundary condition in mol/m^2" << endl;
            outFile << nodeID << " " << old_concentration[nodeID] << endl << endl;
            break;
          }
          else
            {continue;}
        }
      }

      else
      {
        if(i % 100 == 0)
        {
          outFile << x << " " << y << " " << 
          displacement[2*nodeID + 1] << " " << 
          old_concentration[nodeID] << " " << 
          damage[nodeID] << endl;
        }
      }
    }
    
    if(i % 100 == 0)
      {outFile << endl << endl;}
  }
  // End load-step loop and indicate the path of the exported simulation
  outFile.close();
  cout << "Data exported to " << outputPath << endl;
    
    }
  }
}
