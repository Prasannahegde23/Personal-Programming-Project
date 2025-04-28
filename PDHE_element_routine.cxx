#include <iostream>
#include <vector>
#include <cmath>
#include "PDHE_element_routine.h" 
#include "PDHE_material_routine.h"
#include "PDHE_boundary_effects.h"
  
using namespace std;

std::vector<double> eta_vec(std::vector<double> &disp_of_center, std::vector<double> &disp_of_j) 
{
    std::vector<double> result(2);
    for(int i=0; i<2; i++)
    {
        result[i] = disp_of_j[i] - disp_of_center[i];
    }
    return result;
}

std::vector<double> xi_vec(double x, double y, double neigh_x, double neigh_y) 
{   
    std::vector<double> result(2);
    result[0] = neigh_x - x;
    result[1] = neigh_y - y;
    return result;
}

double mod_xi(const std::vector<double> &result)
{
    return sqrt((result[0]*result[0]) + (result[1]*result[1]));
}

elementroutinehydrogen element_routine_hydrogen(int nodeID, double *modelCoord,double x, double y,const int *neighborhoodList, int neighIndex, int numNeighbors, double m_horizon, std::vector<double> &old_concentration, double concentration_nodeID, double time_step_size_EFM, double dh, double Volume_i, double *volume)
{
    double C_dot = 0.0;
    double individual_C_dot;
    elementroutinehydrogen finalresult;
    for(int n=0; n < numNeighbors; ++n)
    {
        int neighborID = neighborhoodList[neighIndex++];
        double neighborVolume = volume[neighborID];
        double neigh_x = modelCoord[3*neighborID];
        double neigh_y = modelCoord[3*neighborID + 1];
        

        std::vector<double> xi;
        xi =  xi_vec(x, y, neigh_x, neigh_y);
        double mag_xi = mod_xi(xi);
        double concentration_neighborID = old_concentration[neighborID];
        individual_C_dot = material_routine_Hydrogen_conc(mag_xi, neighborVolume, dh, concentration_nodeID, concentration_neighborID);
        C_dot = C_dot + individual_C_dot;
    }

    concentration_nodeID = concentration_nodeID + (C_dot * time_step_size_EFM); //Concentration of point
    finalresult.conc = concentration_nodeID;
    finalresult.neighindex = neighIndex;
    return finalresult;
}

PDResult element_routine_PD(double Volume_i, double *volume, double c, double m_h, double m_horizon, double k_n, double k_t, double m_Sat_Val_Hyd_Conc, double m_Critic_Energy_Rel_Rate, double *modelCoord, double x, double y, int nodeID,const int *neighborhoodList, int neighIndex, int numNeighbors, /*double *displacement*/std::vector<double> &displacement, std::vector<double> &old_concentration, double concentration_nodeID, double m_min_grid_spacing, std::vector<double> &bondFactor/*, std::ostream &logStream*/)
{
    std::vector<double> PD_force_i(2, 0.0);
    PDResult result;
    std::vector<double> PD_force_i_j(2);
    std::vector<double> disp_center(2);
    std::vector<double> disp_j(2);
    std::vector<double> eta(2); 
    std::vector<double> xi(2);
    //std::vector<double> b_d(numNeighbors);
    double damage_val = 0.0;
    double Volume_j = 0.0;
    for(int n=0; n < numNeighbors; ++n)
    {
        int neighborID = neighborhoodList[neighIndex++];
        double neighborVolume = volume[neighborID];
        double neigh_x = modelCoord[3*neighborID];
        double neigh_y = modelCoord[3*neighborID + 1];

        disp_center[0] = displacement[2*nodeID];
        disp_center[1] = displacement[2*nodeID + 1];
        disp_j[0] = displacement[2*neighborID];
        disp_j[1] = displacement[2*neighborID + 1];

        double dx = (modelCoord[3*neighborID] + disp_j[0]) - (x + disp_center[0]);
        double dy = (modelCoord[3*neighborID+1] + disp_j[1]) - (y + disp_center[1]);
        if(std::sqrt(dx*dx + dy*dy) > m_horizon || bondFactor[n]==0.0)
          continue;  // permanently skip any broken or over-stretched pair
          
        eta = eta_vec(disp_center, disp_j);
        xi = xi_vec(x, y, neigh_x, neigh_y);

        double mag_xi = mod_xi(xi);
        double concentration_neighborID = old_concentration[neighborID];
        PDOutputs out;   

        out = material_routine_PD(c, m_h, m_horizon, k_n, k_t, mag_xi, m_Sat_Val_Hyd_Conc, m_Critic_Energy_Rel_Rate, nodeID, neighborID, n, concentration_nodeID, concentration_neighborID, m_min_grid_spacing, Volume_i, neighborVolume, bondFactor, eta, xi);
        PD_force_i_j[0] = out.PDforce_x;
        PD_force_i_j[1] = out.PDforce_y;
        bondFactor[n] = out.bondVal;
        PD_force_i[0] = PD_force_i[0] + PD_force_i_j[0];
        PD_force_i[1] = PD_force_i[1] + PD_force_i_j[1];
        
        damage_val = damage_val + (bondFactor[n] * neighborVolume);
        Volume_j = Volume_j + neighborVolume;
    }

    result.Px = PD_force_i[0];
    result.Py = PD_force_i[1];
    result.damage = 1 - (damage_val)/(Volume_j);
    result.neighindex = neighIndex;
    result.bondFac = bondFactor;

    return result;
}