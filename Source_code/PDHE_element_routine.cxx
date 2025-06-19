#include <iostream>
#include <vector>
#include <cmath>
#include "PDHE_element_routine.h" 
#include "PDHE_material_routine.h"
#include "PDHE_boundary_effects.h"
   
using namespace std;

// Computes relative displacement vector between two nodes
std::vector<double> eta_vec(std::vector<double> &disp_of_center, std::vector<double> &disp_of_j) 
{
    std::vector<double> result(2);
    for(int i=0; i<2; i++)
    {
        result[i] = disp_of_j[i] - disp_of_center[i];
    }

    return result;
}

// Computes the initial bond vector between two nodes in the undeformed configuration
std::vector<double> xi_vec(double x, double y, double neigh_x, double neigh_y) 
{   
    std::vector<double> result(2);
    result[0] = neigh_x - x;
    result[1] = neigh_y - y;

    return result;
}

// Computes magnitude of 2D vector
double mod_xi(const std::vector<double> &result) // Used for calculating magnitude of "xi" vector
{
    double res = sqrt((result[0]*result[0]) + (result[1]*result[1]));
    return res;
}

//----------------------------------------------------------------------------//
// element_routine_hydrogen: Performs a single explicit diffusion substep for
//                            hydrogen concentration at one node
// Inputs:
//   nodeID               - global ID of the current node
//   modelCoord           - flattened array of undeformed coordinates
//   x, y                 - undeformed position of this node
//   neighborhoodList     - neighborhood array
//   neighIndex           - index into neighborhoodList where neighbor IDs start
//   numNeighbors         - number of neighbors for this node
//   m_horizon            - peridynamic horizon length
//   old_concentration    - reference to vector holding concentrations
//   concentration_nodeID - concentration at current node
//   time_step_size_EFM   - diffusion substep size
//   dh                   - hydrogen diffusion bond constant
//   Volume_i             - volume of this node
//   volume               - pointer to volumes of all nodes
// Outputs:
//   finalresult.conc     - updated concentration
//   finalresult.neighindex - updated index after iterating neighbors
//----------------------------------------------------------------------------//

elementroutinehydrogen element_routine_hydrogen(int nodeID, double *modelCoord,double x, double y,const int *neighborhoodList, int neighIndex, int numNeighbors, double m_horizon, std::vector<double> &old_concentration, double concentration_nodeID, double time_step_size_EFM, double dh, double Volume_i, double *volume)
{
    double C_dot = 0.0;
    double individual_C_dot;
    elementroutinehydrogen finalresult;

    // Loop over each neighbor to compute pairwise diffusion contributions
    for(int n=0; n < numNeighbors; ++n)
    {
        int neighborID = neighborhoodList[neighIndex++];
        double neighborVolume = volume[neighborID];
        double neigh_x = modelCoord[3*neighborID];
        double neigh_y = modelCoord[3*neighborID + 1];
        
        // Compute initial bond vector magnitude
        std::vector<double> xi;
        xi =  xi_vec(x, y, neigh_x, neigh_y);
        double mag_xi = mod_xi(xi);
        double concentration_neighborID = old_concentration[neighborID];

        // Calculate diffusion rate from neighbor j into i
        individual_C_dot = material_routine_Hydrogen_conc(mag_xi, neighborVolume, dh, concentration_nodeID, concentration_neighborID);
        C_dot = C_dot + individual_C_dot; // Accumulated diffusion flux
    }

    // Update concentration by Explicit Euler scheme
    concentration_nodeID = concentration_nodeID + (C_dot * time_step_size_EFM); //Concentration of point
    finalresult.conc = concentration_nodeID;
    finalresult.neighindex = neighIndex;
    return finalresult;
}

//----------------------------------------------------------------------------//
// element_routine_PD: Computes peridynamic internal forces and updates bond
//                    damage for one node
// Inputs:
//   Volume_i               - volume of this node
//   volume                 - pointer to volumes of all nodes
//   c                      - dilatational stiffness constant
//   m_h                    - material thickness
//   m_horizon              - peridynamic horizon length
//   k_n, k_t               - normal and tangential stiffness constants
//   m_Sat_Val_Hyd_Conc     - saturation concentration (for hydrogen embrittlement)
//   m_Critic_Energy_Rel_Rate - critical energy release rate for bond failure
//   modelCoord             - flattened array of undeformed coordinates
//   x, y                   - undeformed position of this node
//   nodeID                 - global index of this node
//   neighborhoodList       - array listing neighbor counts and IDs
//   neighIndex             - starting index into neighborhoodList
//   numNeighbors           - number of neighbors
//   displacement           - current displacements vector (size 2*N)
//   old_concentration      - concentration vector at current time
//   concentration_nodeID   - concentration at this node
//   m_min_grid_spacing     - smallest grid spacing
//   bondFactor             - vector of bond factors (0 to 1)
// Outputs:
//   result.Px, result.Py   - summed internal force components on node
//   result.bondFac         - updated bondFactor for current node
//   result.neighindex      - updated neighbor index after loop
//----------------------------------------------------------------------------//

PDResult element_routine_PD(double Volume_i, double *volume, double c, double m_h, double m_horizon, double k_n, double k_t, double m_Sat_Val_Hyd_Conc, double m_Critic_Energy_Rel_Rate, double *modelCoord, double x, double y, int nodeID,const int *neighborhoodList, int neighIndex, int numNeighbors, /*double *displacement*/std::vector<double> &displacement, std::vector<double> &old_concentration, double concentration_nodeID, double m_min_grid_spacing, std::vector<double> &bondFactor/*, std::ostream &logStream*/)
{
    // Preallocate small vectors for neighbor computations
    std::vector<double> PD_force_i(2, 0.0);
    PDResult result;
    std::vector<double> PD_force_i_j(2);
    std::vector<double> disp_center(2);
    std::vector<double> disp_j(2);
    std::vector<double> eta(2); 
    std::vector<double> xi(2);

    // Loop over neighbors to compute pairwise PD forces
    for(int n=0; n < numNeighbors; ++n)
    {
        int neighborID = neighborhoodList[neighIndex++];
        double neighborVolume = volume[neighborID];
        double neigh_x = modelCoord[3*neighborID];
        double neigh_y = modelCoord[3*neighborID + 1];

        // Relative displacement vectors
        disp_center[0] = displacement[2*nodeID];
        disp_center[1] = displacement[2*nodeID + 1];
        disp_j[0] = displacement[2*neighborID];
        disp_j[1] = displacement[2*neighborID + 1];

        // Check if bond is intact and within horizon
        double dx = (modelCoord[3*neighborID] + disp_j[0]) - (x + disp_center[0]);
        double dy = (modelCoord[3*neighborID+1] + disp_j[1]) - (y + disp_center[1]);
        if(std::sqrt(dx*dx + dy*dy) > m_horizon || bondFactor[n]==0.0)
          continue; // skip broken or out-of-range bonds
        
        // Compute relative vectors
        eta = eta_vec(disp_center, disp_j);
        xi = xi_vec(x, y, neigh_x, neigh_y);

        double mag_xi = mod_xi(xi);
        double concentration_neighborID = old_concentration[neighborID];
        PDOutputs out;   

        // Calculate PD force and bond factor between point i and j
        out = material_routine_PD(c, m_h, m_horizon, k_n, k_t, mag_xi, m_Sat_Val_Hyd_Conc, m_Critic_Energy_Rel_Rate, nodeID, neighborID, n, concentration_nodeID, concentration_neighborID, m_min_grid_spacing, Volume_i, neighborVolume, bondFactor, eta, xi);
        PD_force_i_j[0] = out.PDforce_x;
        PD_force_i_j[1] = out.PDforce_y;
        bondFactor[n] = out.bondVal;

        // accumulate forces in x and y
        PD_force_i[0] = PD_force_i[0] + PD_force_i_j[0];
        PD_force_i[1] = PD_force_i[1] + PD_force_i_j[1];
    }

    // Package results for the caller
    result.Px = PD_force_i[0];
    result.Py = PD_force_i[1];
    result.neighindex = neighIndex;
    result.bondFac = bondFactor;

    return result;
}