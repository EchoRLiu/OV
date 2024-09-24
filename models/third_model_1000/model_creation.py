import tellurium as te
import libsbml

# Generate ODE functions and initial conditions
ode_functions = []
initial_conditions = []

# ODE for C_u
ode_functions.append("C_u' = rho * C_u * (1 - kappa * (C_u + C_i + " + " + ".join([f"C_i{j}" for j in range(1, 1001)]) + " + C_l)) - psi * V * C_u;")
initial_conditions.append("C_u = 1 / (kappa + exp(-rho) * (1/400 - kappa));")

# ODE for C_i
ode_functions.append("C_i' = psi * V * C_u - phi * C_i;")
initial_conditions.append("C_i = 0;")

# ODE and initial conditions for C_i1 to C_i1000
for j in range(1, 1001):
  if j == 1:
    ode_functions.append(f"C_i{j}' = phi * C_i - phi * C_i{j};")
  else:
    ode_functions.append(f"C_i{j}' = phi * C_i{j-1} - phi * C_i{j};")
  initial_conditions.append(f"C_i{j} = 0;")

# ODE for C_l
ode_functions.append("C_l' = phi * C_i1000 - alpha * C_l;")
initial_conditions.append("C_l = 0;")

# ODE for V
ode_functions.append("V' = beta * alpha * C_l - psi * V * C_u - delta * V;")
initial_conditions.append("V = virus_injection;")

# Combine all parts into the model string
with_delay_approx_1000_model = f"""
model vvDD_effect_with_delay_approx_1000

  // ODE functions
  {' '.join(ode_functions)}

  // initial conditions
  // t=0 is when the virus is injected
  {' '.join(initial_conditions)}

  // condition dependent parameters
  virus_injection = 3 * 1E9; // pfu
  
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

end
"""

# Load the model
r = te.loada(with_delay_approx_1000_model)

# convert model back to Antimony / SBML
ant_str_before = r.getAntimony()
sbml_str_before = r.getSBML()

# write xml file
# with open('petab_files/with_delay_approx_1000.xml', 'w') as f:
#     f.write(sbml_str_before)

r.exportToSBML('petab_files/with_delay_approx_1000.xml', current=False)