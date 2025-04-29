#include <iostream>
#include <vector>
#include <cmath> 
#include "PDHE_boundary_effects.h" 
 
using namespace std;

//----------------------------------------------------------------------------//
// surface_correction:
// Inputs:
//   k_n        - nominal normal bond stiffness
//   k_t        - nominal tangential bond stiffness
//   m_horizon  - peridynamic horizon length
//   V_i, V_j   - volumes of the interacting nodes i and j
//   m_h        - thickness of the material (used to compute reference volume)
// Returns:
//   PDParameter struct with corrected k_n and k_t values
//----------------------------------------------------------------------------//

PDParameter surface_correction(double k_n, double k_t, double m_horizon, double V_i, double V_j, double m_h)
{
    double surf_correc_factor;
    double V_0 = m_h*M_PI*pow(m_horizon,2);
    PDParameter result;
    surf_correc_factor = (2 * V_0)/(V_i + V_j);
    result.K_n = surf_correc_factor * k_n;
    result.K_t = surf_correc_factor * k_t;
    return result; // returning corrected PD parameter
}

//----------------------------------------------------------------------------//
// volume_correction:
// Inputs:
//   PD_force         - original pairwise force vector
//   m_horizon        - peridynamic horizon
//   mag_xi           - undeformed bond length
//   m_min_grid_spacing - grid spacing
//   nodeID, neighborID - IDs
// Returns:
//   Modified PD_force vector after volume correction
//----------------------------------------------------------------------------//

std::vector<double> volume_correction(std::vector<double> &PD_force, double m_horizon, double mag_xi, double m_min_grid_spacing, double nodeID, double neighborID)
{   
    double vol_corr_factor;
    double check = m_horizon - (m_min_grid_spacing/2);
    if (mag_xi < check)
    {
        vol_corr_factor = 1;
        for(int k = 0; k<2; k++)
        {
            PD_force[k] = vol_corr_factor * PD_force[k];
        }
    }

    if (mag_xi > check)
    {
        vol_corr_factor = (m_horizon - mag_xi + (m_min_grid_spacing/2))/m_min_grid_spacing;
        for(int k = 0; k<2; k++)
        {
            PD_force[k] = vol_corr_factor * PD_force[k];
        }
    }
    return PD_force;
}


