#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "PDHE_material_routine.h"
#include "PDHE_boundary_effects.h"
 
using namespace std;

// Computes magnitude of 2D vector
double vector_magnitude(std::vector<double> &v)
{
    double res = sqrt((v[0]*v[0]) + (v[1]*v[1]));
    
    // Used for Unit test:
    /* 
    Used for Unit test:
        cout << "Magnitude of xi: " << res << endl;
    */
    return res;
}

// Computes addition of 2D vector
std::vector<double> addVectors(std::vector<double> &a, std::vector<double> &b)  
{
    std::vector<double> result(2);
    for(int i=0; i<2; i++) 
    {result[i] = a[i] + b[i];}
    
    // Used for Unit test:
    /* 
        cout << "Result of vector: " << endl;
        cout << "x " << result[0] << endl;
        cout << "y " << result[1] << endl;
    */
    return result;
}

//----------------------------------------------------------------------------//
// material_routine_Hydrogen_conc: Compute hydrogen diffusion rate between a pair
// of nodes using Fick-like formulation
// Inputs:
//   mag_xi                   - initial bond length
//   neighborVolume           - volume of neighbor node
//   dh                       - diffusion coefficient constant for grain boundary
//   concentration_nodeID     - hydrogen concentration at current node
//   concentration_neighborID - hydrogen concentration at neighbor node
// Returns:
//   C_dot: net concentration rate
//----------------------------------------------------------------------------//

double material_routine_Hydrogen_conc(double mag_xi, double neighborVolume, double dh, double concentration_nodeID, double concentration_neighborID)
{
    double C_dot;

    if(mag_xi < 1e-3) 
    {C_dot = 0.0;}

    else
    {C_dot = dh*(concentration_neighborID - concentration_nodeID)*neighborVolume/(mag_xi*mag_xi);}
    return C_dot;
}

//----------------------------------------------------------------------------//
// material_routine_PD: Compute peridynamic pair force and update bond factor
// considering hydrogen embrittlement and boundary effects
// Inputs:
//   c, m_h, m_horizon        - peridynamic constants and material thickness
//   k_n, k_t                 - nominal normal/tangential bond stiffness
//   mag_xi                   - undeformed bond length
//   m_Sat_Val_Hyd_Conc       - saturation concentration for embrittlement
//   m_Critic_Energy_Rel_Rate - critical energy release rate for bond failure
//   nodeID, neighborID       - global IDs of the interacting nodes
//   n                        - index in bondFactor vector for this neighbor
//   concentration_nodeID     - current node concentration
//   concentration_neighborID - neighbor concentration
//   m_min_grid_spacing       - smallest cell spacing
//   Volume_i, neighborVolume - volumes for surface correction
//   b_d                      - reference to bondFactor array for this node
//   eta, xi                  - relative displacement and initial bond vectors
// Outputs in PDOutputs:
//   PDforce_x, PDforce_y     - force components on this node from this bond
//   bondVal                  - updated bond factor (0 or 1)
//----------------------------------------------------------------------------//

PDOutputs material_routine_PD(double c, double m_h, double m_horizon, double k_n, double k_t, double mag_xi, double m_Sat_Val_Hyd_Conc, double m_Critic_Energy_Rel_Rate, int nodeID, int neighborID, int n, double concentration_nodeID, double concentration_neighborID, double m_min_grid_spacing, double Volume_i, double neighborVolume, std::vector<double> &b_d, std::vector<double> &eta, std::vector<double> &xi)
{
    PDOutputs out;

    std::vector<double> eta_n(2);
    std::vector<double> eta_t(2);
    double kn; double kt;
    double bond_factor;

    PDParameter pdparameter;
    // Apply boundary correction to stiffness using surface_correction
    pdparameter = surface_correction(k_n, k_t, m_horizon, Volume_i, neighborVolume, m_h); // boundary effects correction
    kn = pdparameter.K_n;
    kt = pdparameter.K_t;

    // Decompose relative displacement into normal and tangential components
    double dot = (eta[0]*xi[0]) + (eta[1]*xi[1]);
    for(int i=0; i<2; i++)
    {
        eta_n[i] = (dot * xi[i])/(mag_xi*mag_xi); //eta_n
        eta_t[i] = eta[i] - eta_n[i]; //eta_t
    }
   
   // Calculation of components of PD forces
   std::vector<double> PD_force(2);
    for(int i=0; i<2; i++)
    {PD_force[i] = ((kn * eta_n[i]) + (kt * eta_t[i]))/mag_xi;}

    std::vector<double> result(2);
    result = addVectors(xi, eta);
    
    // Compute current stretch and critical stretch
    double s_c0 = sqrt((4 * m_Critic_Energy_Rel_Rate)/(c * m_h * pow(m_horizon,4))); // critical stretch
    double s = (vector_magnitude(result) - mag_xi)/mag_xi ; // current stretch
    
    // Coupling of mechanics-chemo including hydrogen coverage
    if(m_Sat_Val_Hyd_Conc > 1e-12)
    {
        double Phi_i =  concentration_nodeID/m_Sat_Val_Hyd_Conc; // Hydrogen coverage of i
        double Phi_j =  concentration_neighborID/m_Sat_Val_Hyd_Conc; // Hydrogen coverage of j
        double Phi = 0.5 * (Phi_i + Phi_j);
        s_c0 = s_c0 * (1 - (1.0467*Phi) + (0.1687*pow(Phi,2)));
        Phi = std::clamp(Phi, 0.0, 1.0);
    }


    // Scalar factor b_d which determines bond breakage
    if (s >= s_c0 || b_d[n] == 0.0)
    {
        bond_factor = 0.0;
        for(int i=0; i<2; i++)
        { PD_force[i] = 0.0;}
    }
    else 
    {
        bond_factor = 1.0;
        PD_force = volume_correction(PD_force, m_horizon, mag_xi, m_min_grid_spacing, nodeID, neighborID);
    }

    out.PDforce_x = PD_force[0];
    out.PDforce_y = PD_force[1];
    out.bondVal = bond_factor;

    return out; 
}