#ifndef PDHE_ELEMENT_ROUTINE_H
#define PDHE_ELEMENT_ROUTINE_H
 
struct PDResult 
{
    double Px;
    double Py;
    double damage;
    int neighindex;
    std::vector<double> bondFac;
};
  
struct elementroutinehydrogen 
{ 
    double conc;
    int neighindex;
};

elementroutinehydrogen element_routine_hydrogen(int nodeID, double *modelCoord,double x, double y, 
    const int *neighborhoodList, int neighIndex, int numNeighbors, double m_horizon, std::vector<double> &old_concentration, 
    double concentration_nodeID, double time_step_size_EFM, double dh, double Volume_i, double *volume);


PDResult element_routine_PD(double Volume_i, double *volume, double c, double m_h, double m_horizon, double k_n, double k_t, 
    double m_Sat_Val_Hyd_Conc, double m_Critic_Energy_Rel_Rate, double *modelCoord, double x, double y, 
    int nodeID,const int *neighborhoodList, int neighIndex, int numNeighbors, /*double *displacement*/std::vector<double> &displacement, std::vector<double> &old_concentration, 
    double concentration_nodeID, double m_min_grid_spacing, std::vector<double> &bondFactor/*, std::ostream &logStream*/);

#endif