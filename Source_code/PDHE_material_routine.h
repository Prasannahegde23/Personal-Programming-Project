#ifndef PDHE_MATERIAL_ROUTINE_H
#define PDHE_MATERIAL_ROUTINE_H

using namespace std;
 
struct PDOutputs 
{
    double PDforce_x;
    double PDforce_y;
    double bondVal; 
};

double material_routine_Hydrogen_conc(double mag_xi, double neighborVolume, 
                            double dh, double concentration_nodeID, double concentration_neighborID);

PDOutputs material_routine_PD(double c, double m_h, double m_horizon, double k_n, double k_t, double mag_xi, double m_Sat_Val_Hyd_Conc, double m_Critic_Energy_Rel_Rate, 
    int nodeID, int neighborID, int m, double concentration_nodeID, double concentration_neighborID, double m_min_grid_spacing,
    double Volume_i, double neighborVolume, std::vector<double> &b_d, std::vector<double> &eta, std::vector<double> &xi);

double vector_magnitude(std::vector<double> &v);

std::vector<double> addVectors(std::vector<double> &a, std::vector<double> &b);

#endif 