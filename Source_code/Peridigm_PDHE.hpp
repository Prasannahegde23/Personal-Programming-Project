#ifndef PERIDIGM_PDHE_HPP
#define PERIDIGM_PDHE_HPP

#include "Peridigm_Material.hpp"
#include "Peridigm_Field.hpp"
 
namespace PeridigmNS {

  class PDHE : public PeridigmNS::Material 
  {
  public: 
    // Constructor s
    PDHE(const Teuchos::ParameterList& params);

    // Destructor
    virtual ~PDHE(){} 

    // These override the pure virtual methods in Material:
    virtual std::string Name() const;
    virtual double YoungsModulus() const;
    virtual double GBdiffCoeff() const;
    virtual double SatValHydConc() const;
    virtual double CriticEnergyRelRate() const;
    virtual double Density() const;
    virtual double PoissonsRatio() const;
    virtual double Horizon() const;
    virtual int Nt() const;
    virtual int Nh() const;
    virtual int Captureloadsteps() const;
    virtual double h() const;
    virtual double minGridSpacing() const;
    virtual bool BCtest() const;
    

    virtual std::vector<int> FieldIds() const;

    // Compute force
    virtual void computeForce(const double dt,
                                const int numOwnedPoints,
                                const int* ownedIDs,
                                const int* neighborhoodList,
                                PeridigmNS::DataManager& dataManager)const;
    
    virtual double BulkModulus() const;
    virtual double ShearModulus() const;

  private:
    // Material Parameters
    double m_Youngs_Modulus;
    double m_GB_Diff_Coeff;             // Grain boundary diffusion coefficient
    double m_Sat_Val_Hyd_Conc;          // Saturated value of hydrogen concentration
    double m_Critic_Energy_Rel_Rate;    // Critical energy release rate
    double m_density;
    double m_poissons_ratio;

    // Input Parameters
    double m_horizon;
    int N_t;                          // No. of load steps
    int N_h;                          // No. of steps for hydrogen concentration
    int N;                            // To Capture and save the simulation frame from N load steps

    // Geometrical Parameters
    double m_h;                         // thickness
    double m_min_grid_spacing;

    // Field IDs to access from the DataManager
    int m_modelCoordinatesFieldId;   // original/reference coordinates
    int m_coordinatesFieldId;        // current coordinates
    int m_volumeFieldId;             // per-node volume
    int m_concentrationFieldId;      // Used for hydrogen concentration
    int m_damageFieldId;             // Used to store Damage value
    int m_bodyForceFieldId;       // Used to store Body force value

    // List of all field IDs
    std::vector<int> m_fieldIds;

    // List of file names
    std::string m_top_displacement;
    std::string m_bottom_displacement;
    std::string m_crack_top;
    std::string m_crack_bottom;
    std::string m_concentration_file;

    bool BC_status; // Used to access for Boundary condition test

    mutable std::vector<std::vector<double>> m_bondFactor;
    mutable std::vector<double> damage;
  };

}
extern "C" {
  PeridigmNS::Material* create(const Teuchos::ParameterList& params);
}

#endif
