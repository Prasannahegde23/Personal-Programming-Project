#ifndef PDHE_BOUNDARY_EFFECTS_H
#define PDHE_BOUNDARY_EFFECTS_H
 
struct PDParameter
{
    double K_n;
    double K_t;
}; 
 
PDParameter surface_correction(double k_n, double k_t,
                          double m_horizon, double V_i, double V_j, double m_h);

std::vector<double> volume_correction(std::vector<double> &PD_force, double m_horizon, double mag_xi, 
                            double m_min_grid_spacing, double nodeID, double neighborID);



#endif
