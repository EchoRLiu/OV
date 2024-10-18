import tellurium as te
import libsbml

"""
the fourth model includes the temporal infection rate and delayed lysis of infected cells
"""

delayed_infection_model = """
model delayed_infection

  // ODE functions
  C_u' = rho * C_u * (1 - kappa * (C_u + C_i + C_i1 + C_i2 + C_i3 + C_i4 + C_i5 + C_l)) - psi * V * C_u;
  C_i' = psi * V * C_u - phi * C_i;
  C_i1' = phi * C_i - phi * C_i1;
  C_i2' = phi * C_i1 - phi * C_i2;
  C_i3' = phi * C_i2 - phi * C_i3;
  C_i4' = phi * C_i3 - phi * C_i4;
  C_i5' = phi * C_i4 - phi * C_i5;
  C_l' = phi * C_i5 - alpha * C_l;
  V' = beta * alpha * C_l - psi * V * C_u - delta * V;

  // initial conditions
  // t=0 is when the virus is injected
  C_u = 1 / (kappa + exp(-rho) * (1/400 - kappa)); // number of uninfected tumor cells at t=0
  // 400 cells/nL is the number of tumor injected at t=-1
  C_i = 0; // number of infected tumor cells at t=0
  C_i1 = 0; // number of infected tumor cells at t=0
  C_i2 = 0; // number of infected tumor cells at t=0
  C_i3 = 0; // number of infected tumor cells at t=0
  C_i4 = 0; // number of infected tumor cells at t=0
  C_i5 = 0; // number of infected tumor cells at t=0
  C_l = 0; // number of infected tumor cells that are going to be lysed at t=0
  V = virus_injection; // number of virus injected at t=0

  // condition dependent parameters
  virus_injection = 3 * 1E9; // pfu

  // temporal infection rate
  psi := psi_0 * exp(- tau * time);
  
  // parameters to be estimated
  rho = 0.01;
  kappa = 0.01;
  psi_0 = 0.01;
  tau = 0.01;
  phi = 0.01;
  alpha = 0.01;
  beta = 0.01;
  delta = 0.01;
  // the scaling parameter should handle the scaling from cell number to volume of tumor

  // it seems like lambda is reserved and hence cannot be used as a parameter name

end
"""

# Load the model
r = te.loada(delayed_infection_model)

# convert model back to Antimony / SBML
ant_str_before = r.getAntimony()
sbml_str_before = r.getSBML()

# write xml file
# with open('petab_files/delayed_infection.xml', 'w') as f:
#     f.write(sbml_str_before)

r.exportToSBML('petab_files/delayed_infection.xml', current=False)