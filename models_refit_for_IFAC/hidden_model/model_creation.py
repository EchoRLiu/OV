import tellurium as te
import libsbml

"""
the fifth model includes the delayed lysis model, where the infected cells can still divie into more infected cells
with five hidden state
without the following:
  # // temporal infection rate
  # psi := psi_0 * exp(- tau * time);
"""

# Generate ODE functions and initial conditions
ode_functions = []
initial_conditions = []
n_hidden = 5

# ODE for C_u
ode_functions.append("C_u' = rho * C_u * (1 - (C_u + C_i + " + " + ".join([f"C_i{j}" for j in range(1, n_hidden + 1)]) + " + C_l)/kappa) - psi * V * C_u;")
initial_conditions.append("C_u = 1 / (1/kappa + exp(-rho) * (1/400 - 1/kappa));")

# ODE for C_i
ode_functions.append("C_i' = rho * C_i * (1 - (C_u + C_i + " + " + ".join([f"C_i{j}" for j in range(1, n_hidden + 1)]) + " + C_l)/kappa) + psi * V * C_u - phi * C_i;")
initial_conditions.append("C_i = 0;")

# ODE and initial conditions for C_i1 to C_i{n_hidden}
for j in range(1, n_hidden + 1):
  if j == 1:
    ode_functions.append(f"C_i{j}' = rho * C_i{j} * (1 - (C_u + C_i + " + " + ".join([f"C_i{k}" for k in range(1, n_hidden + 1)]) + f" + C_l)/kappa) + phi * C_i - phi * C_i{j};")
  else:
    ode_functions.append(f"C_i{j}' = rho * C_i{j} * (1 - (C_u + C_i + " + " + ".join([f"C_i{k}" for k in range(1, n_hidden + 1)]) + f" + C_l)/kappa) + phi * C_i{j-1} - phi * C_i{j};")
  initial_conditions.append(f"C_i{j} = 0;")

# ODE for C_l
ode_functions.append(f"C_l' = rho * C_l * (1 - (C_u + C_i + " + " + ".join([f"C_i{j}" for j in range(1, n_hidden + 1)]) + f" + C_l)/kappa) + phi * C_i{n_hidden} - alpha * C_l;")
initial_conditions.append("C_l = 0;")

# ODE for V
ode_functions.append("V' = beta * alpha * C_l - psi * V * C_u - delta * V;")
initial_conditions.append("V = virus_injection;")

dividing_infected_cells_model = f"""
model dividing_infected_cells

  // ODE functions
  {' '.join(ode_functions)}

  // initial conditions
  // t=0 is when the virus is injected
  {' '.join(initial_conditions)}

  // condition dependent parameters
  virus_injection = 1 * 1E9; // pfu
  
  // parameters to be estimated
  rho = 0.01;
  kappa = 0.01;
  psi = 0.01;
  phi = 0.01;
  alpha = 0.01;
  beta = 0.01;
  delta = 0.01;
  // the scaling parameter should handle the scaling from cell number to volume of tumor

  // it seems like lambda is reserved and hence cannot be used as a parameter name

  // Event to set C_u to zero when C_u falls below 1e-8
  at (C_u < 1e-2): C_u = 0;
end
"""

# Load the model
r = te.loada(dividing_infected_cells_model)

# convert model back to Antimony / SBML
ant_str_before = r.getAntimony()
sbml_str_before = r.getSBML()

# write xml file
# with open('petab_files/dividing_infected_cells_v3.xml', 'w') as f:
#     f.write(sbml_str_before)

r.exportToSBML('petab_files/dividing_infected_cells_event.xml', current=False)