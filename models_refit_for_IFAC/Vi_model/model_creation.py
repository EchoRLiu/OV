import tellurium as te
import libsbml

Vi_model_event = """

model Vi_model_event

  // ODE functions
  C_u' = rho * C_u * (1 - (C_u + C_i) / kappa) - psi_1 * V_e * C_u;
  C_i' = rho * C_i * (1 - (C_u + C_i) / kappa) + psi_1 * V_e * C_u - alpha * V_i;
  V_i' = psi_1 * V_e * C_u + psi_2 * V_e * C_i + zeta * C_i - alpha * V_i * V_i / C_i;
  V_e' = alpha * V_i * V_i / C_i - psi_1 * V_e * C_u - psi_2 * V_e * C_i - delta * V_e;

  // initial conditions
  // t=0 is when the virus is injected
  C_u = 1 / (1/kappa + exp(-rho) * (1/400 - 1/kappa)); // number of uninfected tumor cells at t=0
  // 400 cells/nL is the number of tumor injected at t=-1
  C_i = 1e-12; // number of infected tumor cells at t=0
  V_i = 1e-12; // chance of lysis at t=0
  V_e = 1e-12 + virus_injection; // number of virus injected at t=0

  // condition dependent parameters
  virus_injection = 1 * 1E9; // pfu
  
  // parameters to be estimated
  rho = 0.01;
  kappa = 0.01;
  psi_1 = 0.01;
  psi_2 = 0.01;
  zeta = 0.01;
  alpha = 0.01;
  delta = 0.01;
  // the scaling parameter should handle the scaling from cell number to volume of tumor

  // it seems like lambda is reserved and hence cannot be used as a parameter name

  // Event to set C_u to zero when C_u falls below 1e-8
  at (C_u < 1e-2): C_u = 0;
end
"""

# Load the model
r = te.loada(Vi_model_event)

# convert model back to Antimony / SBML
ant_str_before = r.getAntimony()
sbml_str_before = r.getSBML()

# write xml file
# with open('petab_files/Vi_model_event.xml', 'w') as f:
#     f.write(sbml_str_before)

r.exportToSBML('petab_files/Vi_model_event.xml', current=False)