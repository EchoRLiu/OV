import tellurium as te
import libsbml

init_model = """
model vvDD_effect

  // ODE functions
  C_u' = rho * C_u * (1 - kappa * (C_u + C_i)) - psi_1 * V * C_u;
  C_i' = psi_1 * V * C_u - alpha * C_i;
  V' = beta * alpha * C_i - psi_2 * V * C_u - delta * V;

  // initial conditions
  // t=0 is when the virus is injected
  C_u = 1/(kappa + exp(-rho + ln(1/(4 * 1E5) - kappa))); // 400 cells/nL but the data is in volume of tumor, so we need to convert it to volume of tumor
  C_i = 0;
  V = virus_injection;

  // condition dependent parameters
  virus_injection = 3 * 1E9; // pfu
  
  // parameters to be estimated
  rho = 0.01;
  kappa = 0.01;
  psi_1 = 0.01;
  psi_2 = 0.01;
  beta = 0.01;
  alpha = 0.01;
  delta = 0.01;

  // it seems like lambda is reserved and hence cannot be used as a parameter name

end
"""

# Load the model
r = te.loada(init_model)

# convert model back to Antimony / SBML
ant_str_before = r.getAntimony()
sbml_str_before = r.getSBML()

# write xml file
# with open('petab_files/init_model.xml', 'w') as f:
#     f.write(sbml_str_before)

r.exportToSBML('petab_files/init_model.xml', current=False)