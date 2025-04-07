#include <iostream>
#include <vector>
#include <cmath>
#include "PDHE_element_routine.h" 
#include "PDHE_material_routine.h"
#include "PDHE_boundary_effects.h"

using namespace std;

double computeDistance(double point_i[3], double point_j[3]) 
{
    double distance = sqrt((point_i[0] - point_j[0]) * (point_i[0] - point_j[0]) +
                    (point_i[1] - point_j[1]) * (point_i[1] - point_j[1]) +
                    (point_i[2] - point_j[2]) * (point_i[2] - point_j[2]));

    return distance;
}

std::vector<double> eta_vec(std::vector<double> disp_of_center, std::vector<double> disp_of_j) 
{
    std::vector<double> result(3);
    for(int i=0; i<3; i++)
    {
        result[i] = disp_of_j[i] - disp_of_center[i];
    }
    return result;
}

std::vector<double> xi_vec(double x, double y, double z, double neigh_x, double neigh_y, double neigh_z) 
{   
    std::vector<double> result(3);
    result[0] = x - neigh_x;
    result[1] = y - neigh_y;
    result[2] = z - neigh_z;
    return result;
}

double mod_xi(const std::vector<double>& result)
{
    return sqrt((result[0]*result[0]) + (result[1]*result[1]) + (result[2]*result[2]));
}

elementroutinehydrogen element_routine_hydrogen(int nodeID, double *currentCoord,double x, double y, double z,const int *neighborhoodList, int neighIndex, int numNeighbors, double m_horizon, double *concentration, double concentration_nodeID, double time_step_size_EFM, double dh, double Volume_i, double *volume)
{
    double C_dot = 0.0;
    double individual_C_dot;
    elementroutinehydrogen finalresult;
    for(int n=0; n < numNeighbors; ++n)
    {
        int neighborID = neighborhoodList[neighIndex++];
        double neighborVolume = volume[neighborID];
        double neigh_x = currentCoord[3*neighborID];
        double neigh_y = currentCoord[3*neighborID + 1];
        double neigh_z = currentCoord[3*neighborID + 2];

        std::vector<double> xi;
        xi =  xi_vec(x, y, z, neigh_x, neigh_y, neigh_z);
        double mag_xi = mod_xi(xi);
        double concentration_neighborID = concentration[neighborID];
        individual_C_dot = material_routine_Hydrogen_conc(mag_xi, neighborVolume, dh, concentration_nodeID, concentration_neighborID);
        C_dot = C_dot + individual_C_dot;
    }

    concentration[nodeID] = concentration[nodeID] + (C_dot * time_step_size_EFM); //Concentration of point
    finalresult.conc = concentration[nodeID];
    finalresult.neighindex = neighIndex;
    return finalresult;
}

PDResult element_routine_PD(double Volume_i, double *volume, double c, double m_horizon, double k_n, double k_t, double m_Sat_Val_Hyd_Conc, double m_Critic_Energy_Rel_Rate, double *currentCoord, double x, double y, double z, int nodeID,const int *neighborhoodList, int neighIndex, int numNeighbors, double *displacement, double *concentration, double concentration_nodeID, double m_min_grid_spacing)
{
    std::vector<double> PD_force_i(3, 0.0);
    PDResult result;
    std::vector<double> PD_force_i_j(3);
    std::vector<double> b_d(numNeighbors);
    double damage_val = 0.0;
    double Volume_j = 0.0;
    for(int n=0; n < numNeighbors; ++n)
    {
        int neighborID = neighborhoodList[neighIndex++];
        double neighborVolume = volume[neighborID];
        double neigh_x = currentCoord[3*neighborID];
        double neigh_y = currentCoord[3*neighborID + 1];
        double neigh_z = currentCoord[3*neighborID + 2];

        std::vector<double> disp_center(3);
        std::vector<double> disp_j(3);
        disp_center[0] = displacement[3*nodeID];
        disp_center[1] = displacement[3*nodeID + 1];
        disp_center[2] = displacement[3*nodeID + 2];
        disp_j[0] = displacement[3*neighborID];
        disp_j[1] = displacement[3*neighborID + 1];
        disp_j[2] = displacement[3*neighborID + 2];

        std::vector<double> eta(3); std::vector<double> xi(3);
        eta = eta_vec(disp_center, disp_j);
        xi = xi_vec(x, y, z, neigh_x, neigh_y, neigh_z);
        double mag_xi = mod_xi(xi);
        double concentration_neighborID = concentration[neighborID];
        PDOutputs out;

        out = material_routine_PD(c, m_horizon, k_n, k_t, mag_xi, m_Sat_Val_Hyd_Conc, m_Critic_Energy_Rel_Rate, nodeID, neighborID, n, concentration_nodeID, concentration_neighborID, m_min_grid_spacing, Volume_i, neighborVolume, b_d, eta, xi);
        PD_force_i_j[0] = out.PDforce_x;
        PD_force_i_j[1] = out.PDforce_y;
        PD_force_i_j[2] = out.PDforce_z;
        b_d[n] = out.bondVal;

        for(int k=0; k<3; k++)
        {
            PD_force_i[k] = PD_force_i[k] + PD_force_i_j[k];
        }

        damage_val = damage_val + (b_d[n] * neighborVolume);
        Volume_j = Volume_j + neighborVolume;
    }
    cout << "PD_force[0] in element routine returning value: " << PD_force_i[0] << endl;
    cout << "PD_force[1] in element routine returning value: " << PD_force_i[0] << endl;
    cout << "PD_force[2] in element routine returning value: " << PD_force_i[2] << endl;
    result.Px = PD_force_i[0];
    result.Py = PD_force_i[1];
    result.Pz = PD_force_i[2];
    result.damage = 1 - (damage_val)/(Volume_j);

    return result;
}