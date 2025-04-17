#include <iostream>
#include <cmath>
#include <vector>
#include "PDHE_material_routine.h"
#include "PDHE_boundary_effects.h"
 
using namespace std;
 
double vector_magnitude(std::vector<double> &v)
{
    return sqrt((v[0]*v[0]) + (v[1]*v[1]));
}

std::vector<double> addVectors(std::vector<double> &a, std::vector<double> &b)  
{
    std::vector<double> result(2);
    for(int i=0; i<2; i++) 
    {result[i] = a[i] + b[i];}
    return result;
}

double material_routine_Hydrogen_conc(double mod_xi, double neighborVolume, double dh, double concentration_nodeID, double concentration_neighborID)
{
    double C_dot = dh*(concentration_neighborID - concentration_nodeID)*neighborVolume/mod_xi;
    return C_dot;
}

PDOutputs material_routine_PD(double c, double m_h, double m_horizon, double k_n, double k_t, double mag_xi, double m_Sat_Val_Hyd_Conc, double m_Critic_Energy_Rel_Rate, int nodeID, int neighborID, int n, double concentration_nodeID, double concentration_neighborID, double m_min_grid_spacing, double Volume_i, double neighborVolume, std::vector<double> &b_d, std::vector<double> &eta, std::vector<double> &xi)
{
    std::vector<double> eta_n(2);
    std::vector<double> eta_t(2);
    double kn; double kt;
    
    PDParameter pdparameter;
    pdparameter = surface_correction(k_n, k_t, m_horizon, Volume_i, neighborVolume, m_h); // boundary effects correction
    kn = pdparameter.K_n;
    kt = pdparameter.K_t;
    // Calculation of eta_n and eta_t

    double dot = (eta[0]*xi[0]) + (eta[1]*xi[1]) + eta[2]*xi[2];

    for(int i=0; i<2; i++)
    {
        //eta_n[i] = ((pow(eta[0],2) + pow(eta[1],2) + pow(eta[2],2))/pow(mag_xi,2))*xi[i]; //eta_n
        eta_n[i] = (dot * xi[i])/(mag_xi*mag_xi); //eta_n
        eta_t[i] = eta[i] - eta_n[i]; //eta_t
    }
   
   // Calculation of components of PD forces
   std::vector<double> PD_force(2);

    for(int i=0; i<2; i++)
    {PD_force[i] = ((kn * eta_n[i]) + (kt * eta_t[i]))/mag_xi;}

    std::vector<double> result(2);
    result = addVectors(xi, eta);
    
    double s_c0 = sqrt((4 * m_Critic_Energy_Rel_Rate)/(c * m_h * pow(m_horizon,4))); // critical stretch
    double s = (vector_magnitude(result) - vector_magnitude(xi))/vector_magnitude(xi) ;
    
    // Coupling of mechanics-chemo
    if(m_Sat_Val_Hyd_Conc > 0.0)
    {
        double Phi_i =  concentration_nodeID/m_Sat_Val_Hyd_Conc; // Hydrogen coverage of i
        double Phi_j =  concentration_neighborID/m_Sat_Val_Hyd_Conc; // Hydrogen coverage of j
        double Phi = 0.5 * (Phi_i + Phi_j);
        s_c0 = s_c0 * (1 - (1.0467*Phi) + (0.1687*pow(Phi,2)));
    }


    // Scalar factor b_d which determines bond breakage
    if (s < s_c0)
        {b_d[n] = 1.0;}
    else 
        {b_d[n] = 0.0;}

    if(b_d[n] == 0.0)
    {
        for(int i=0; i<2; i++)
        { PD_force[i] = 0.0;}
    }

    else
    {PD_force = volume_correction(PD_force, m_horizon, mag_xi, m_min_grid_spacing, nodeID, neighborID);}

    PDOutputs out;

    out.PDforce_x = PD_force[0];
    out.PDforce_y = PD_force[1];
    out.bondVal = b_d[n];

    return out; 
}