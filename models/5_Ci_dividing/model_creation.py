import tellurium as te
import libsbml

"""
the fifth model includes the delayed lysis model, where the infected cells can still divie into more infected cells
with only one infected state
without the following:
  # // temporal infection rate
  # psi := psi_0 * exp(- tau * time);
"""

dividing_infected_cells_model = """
model dividing_infected_cells

  // ODE functions
  C_u' = rho * C_u * (1 - kappa * (C_u + C_i)) - psi * V * C_u;
  C_i' = rho * C_i * (1 - kappa * (C_u + C_i)) + psi * V * C_u - alpha * C_i;
  V' = beta * alpha * C_i - psi * V * C_u - delta * V;

  // initial conditions
  // t=0 is when the virus is injected
  C_u = 1 / (kappa + exp(-rho) * (1/400 - kappa)); // number of uninfected tumor cells at t=0
  // 400 cells/nL is the number of tumor injected at t=-1
  C_i = 0; // number of infected tumor cells at t=0
  V = virus_injection; // number of virus injected at t=0

  // condition dependent parameters
  virus_injection = 3 * 1E9; // pfu
  
  // parameters to be estimated
  rho = 0.01;
  kappa = 0.01;
  psi = 0.01;
  alpha = 0.01;
  beta = 0.01;
  delta = 0.01;
  // the scaling parameter should handle the scaling from cell number to volume of tumor

  // it seems like lambda is reserved and hence cannot be used as a parameter name

end
"""

# Load the model
r = te.loada(dividing_infected_cells_model)

# convert model back to Antimony / SBML
ant_str_before = r.getAntimony()
sbml_str_before = r.getSBML()

# write xml file
# with open('petab_files/dividing_infected_cells.xml', 'w') as f:
#     f.write(sbml_str_before)

r.exportToSBML('petab_files/dividing_infected_cells.xml', current=False)