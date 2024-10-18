import numpy as np
import pandas as pd
import scipy
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import os
import h5py
import copy

import amici
from petab.C import *
import petab
import petab.C
import pypesto
import pypesto.petab
from pypesto.optimize import minimize
from pypesto.startpoint import uniform
from pypesto.engine import MultiProcessEngine, MultiThreadEngine
from pypesto.optimize.optimizer import FidesOptimizer
import pypesto.optimize as optimize

from pypesto.visualize import waterfall
from pypesto.visualize import parameters
from pypesto.visualize.model_fit import visualize_optimized_model_fit
from pypesto.visualize import profiles

import pypesto.profile as profile
from pypesto.optimize import ScipyOptimizer
from pypesto.profile.options import ProfileOptions

from pypesto.visualize.model_fit import visualize_optimized_model_fit, _get_simulation_rdatas

def hex_to_rgba_gradient(color1, color2, n):
    '''
    Create a gradient in rgba between two hex colors
    '''
    # Convert to rgba
    c1 = matplotlib.colors.to_rgba(matplotlib.colors.hex2color(color1))
    c2 = matplotlib.colors.to_rgba(matplotlib.colors.hex2color(color2))

    return [[(c1[i]*(n-j-1) + c2[i]*j)/(n-1) for i in range(4)] for j in range(n)]

# find the index for cut off based on Chi square distribution CI 95%
def find_cut_off_index(result, ci = 0.95):
    '''
    Find the cut off index for the data based on the Chi square distribution
    '''

    # calculate the chi square distribution
    cut_off_value = scipy.stats.chi2.ppf(ci, 1)

    # find the index
    best_fval = result.optimize_result.list[0].fval

    for i in range(len(result.optimize_result.list)):
        if result.optimize_result.list[i].fval > best_fval + cut_off_value:
            break
    
    return i - 1

def find_cut_off_x_trace(result, ci = 0.95, flatten = True):

    cut_off_value = scipy.stats.chi2.ppf(ci, 1)
    best_fval = result.optimize_result.list[0].fval

    # store the optimized x trace that are below the cut off value
    x_trace_within_cut_off = []
    if flatten:
        for i in range(find_cut_off_index(result, ci)):
            
            fval_trace = result.optimize_result.list[i].history.get_fval_trace()
            x_trace = result.optimize_result.list[i].history.get_x_trace()

            for j in range(len(fval_trace)):
                if fval_trace[j] < best_fval + cut_off_value:
                    x_trace_within_cut_off.append(x_trace[j])
    else:
        for i in range(find_cut_off_index(result, ci)):
            
            fval_trace = result.optimize_result.list[i].history.get_fval_trace()
            x_trace = result.optimize_result.list[i].history.get_x_trace()

            x_trace_within_cut_off_i = []
            for j in range(len(fval_trace)):
                if fval_trace[j] < best_fval + cut_off_value:
                    x_trace_within_cut_off_i.append(x_trace[j])
            x_trace_within_cut_off.append(x_trace_within_cut_off_i)

    return x_trace_within_cut_off


# Plot setting
plt.rcParams['font.size'] = 30

dpi = 100
wid = int(2560/dpi)
hei = int(1600/dpi)

# number of optimization runs
n_runs = 5000

# Define the folder where you want to save the figures
folder_path = "/Users/yuhongliu/Documents/OV/figures/third_model_1000/n"+str(n_runs)+"/"

# If the folder does not exist, create it
if not os.path.exists(folder_path):
    os.makedirs(folder_path)

# optimization
hierarchical = True

petab_yaml = 'petab_files/with_delay_approx_1000.yaml'
petab.validate(petab_yaml)
petab_problem = petab.Problem.from_yaml(petab_yaml)

np.random.seed(500)

problem = pypesto.petab.PetabImporter(
        petab_problem,
        hierarchical=hierarchical,
        model_name=f"WITH_DELAY_APPROX_1000_Model",
    ).create_problem(force_compile=True)

problem.objective.amici_model.setAllStatesNonNegative()

loading = True

if loading:
    # load result history from file
    result = pypesto.store.read_result('optimization_history/n'+ str(n_runs) +'.hdf5')

else:
    # optimize the model
    result = minimize(
        problem=problem,
        optimizer=FidesOptimizer(verbose=False, options={'maxiter': 5000}),
        n_starts=n_runs,
        engine=MultiProcessEngine(),
        # startpoint_method=uniform,
        history_options = pypesto.HistoryOptions(trace_record=True, storage_file='optimization_history/n'+ str(n_runs) +'.hdf5'),
        filename='optimization_history/n'+ str(n_runs) +'.hdf5',
    )

# print result summary
print(result.summary())

parameters_from_result = dict(zip(problem.x_names, result.optimize_result.list[0]['x']))

# Scale all parameters and put them into a dictionary
scaled_parameters = {key: 10**value for key, value in parameters_from_result.items()}

# Print the scaled parameters
print("Scaled parameters:")
for key, value in scaled_parameters.items():
    print(f"{key}: {value}")

return_dict = problem.objective(result.optimize_result.list[0].x, return_dict=True)
rdatas = return_dict['rdatas']
edatas = problem.objective.edatas
x_axis = [edata.id for edata in edatas]
simulation = [rdata.y.reshape(5, -1)[:,0] for rdata in rdatas]
data = [np.array(edata.getObservedData()).reshape(5, -1) for edata in edatas]

# visualize the fitting result by comparing the simulation and the data
# separate into vvDD and ctrl with two subplots

# plot the ctrl data and simulation, simulation[0] (1d array) and data[0] 
# (2d array as it has 5 timepoints [3,4,5,6,7] and 10 samples)
plt.figure(figsize=(10, 8))
for i in range(10):
    plt.plot(np.array([3,4,5,6,7]), data[0][:, i], 
             marker='o', linestyle='-', color='red', alpha=0.3, label='data' if i == 0 else "")
plt.plot(np.array([3,4,5,6,7]), simulation[0], 
         marker='*', linestyle='-', color='red', alpha=1, label='simulation')
plt.xlabel('Time (days)')
plt.ylabel('Tumor Volume ($\mu m^3$)')
plt.title('Tumor Volume vs Time for PBS Conditions')
plt.legend()
plt.grid(False)
plt.xticks([3, 4, 5, 6, 7])
plt.tight_layout()
plt.savefig(folder_path + 'ctrl_fit.pdf', dpi=300, bbox_inches='tight')
plt.show()

# plot the vvDD data and simulation, simulation[1] (1d array) and data[1] 
# (2d array as it has 5 timepoints [3,4,5,6,7] and 10 samples)
plt.figure(figsize=(10, 8))
for i in range(10):
    plt.plot(np.array([3,4,5,6,7]), data[1][:, i], 
             marker='o', linestyle='-', color='blue', alpha=0.3, label='data' if i == 0 else "")
plt.plot(np.array([3,4,5,6,7]), simulation[1], 
         marker='*', linestyle='-', color='blue', alpha=1, label='simulation')
plt.xlabel('Time (days)')
plt.ylabel('Tumor Volume ($\mu m^3$)')
plt.title('Tumor Volume vs Time for vvDD Conditions')
plt.legend()
plt.grid(False)
plt.xticks([3, 4, 5, 6, 7])
plt.tight_layout()
plt.savefig(folder_path + 'vvDD_fit.pdf', dpi=300, bbox_inches='tight')
plt.show()

plt.figure(figsize=(14, 8))

for i in range(10):
    plt.plot(np.array([3,4,5,6,7]), data[0][:, i], 
             marker='o', linestyle='-', color='red', alpha=0.3, label='data - ctrl' if i == 0 else "")
plt.plot(np.array([3,4,5,6,7]), simulation[0], 
         marker='*', linestyle='-', color='red', alpha=1, label='simulation - ctrl')

for i in range(10):
    plt.plot(np.array([3,4,5,6,7]), data[1][:, i], 
             marker='o', linestyle='-', color='blue', alpha=0.3, label='data - vvDD' if i == 0 else "")
plt.plot(np.array([3,4,5,6,7]), simulation[1], 
         marker='*', linestyle='-', color='blue', alpha=1, label='simulation - vvDD')

plt.xlabel('Time (days)')
plt.ylabel('Tumor Volume ($\mu m^3$)')
plt.title('Tumor Volume vs Time for vvDD Conditions')
plt.legend()
plt.grid(False)
plt.xticks([3, 4, 5, 6, 7])
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')  
plt.tight_layout()
plt.savefig(folder_path + 'ctrl_vvDD_fit.pdf', dpi=300, bbox_inches='tight')
plt.show()

# calculate the mean and std of the data within ctrl and vvDD groups and compare again with the simulation

# Calculate mean and standard error for vvDD
vvDD_mean = data[1].mean(axis=1)
vvDD_se = data[1].std(axis=1)
# Calculate mean and standard error for PBS
pbs_mean = data[0].mean(axis=1)
pbs_se = data[0].std(axis=1)

# Plot both provided and calculated data
plt.figure(figsize=(12, 8))
# Plot provided data
plt.errorbar(np.array([3,4,5,6,7]), vvDD_mean, 
             yerr=[vvDD_mean - vvDD_se],
             fmt='o--', color='blue', alpha=0.3, ecolor='lightblue', capsize=5, label='data - vvDD')
plt.errorbar(np.array([3,4,5,6,7]), pbs_mean, 
             yerr=[pbs_mean - pbs_se],
             fmt='o--', color='red', alpha=0.3, ecolor='lightcoral', capsize=5, label='data - ctrl')
# Plot simulated data
plt.plot(np.array([3,4,5,6,7]), simulation[0], 
         marker='*', linestyle='-', color='red', alpha=1, label='simulation - ctrl')
plt.plot(np.array([3,4,5,6,7]), simulation[1], 
         marker='*', linestyle='-', color='blue', alpha=1, label='simulation - vvDD')

plt.xlabel('Time (days)')
plt.ylabel('Tumor Volume ($\mu$ $m^3$)') 
plt.title('Mean Tumor Volume with Standard Deviation vs Time')
plt.legend()
plt.grid(False)
plt.xticks([3, 4, 5, 6, 7])
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')  
plt.tight_layout()
plt.savefig(folder_path + 'mean_tumor_volume_with_std.pdf', dpi=300, bbox_inches='tight')
plt.show()

"""
the C_u, C_i, V dynamics before the injection of virus
for both ctrl and vvDD conditions
"""

import numpy as np

def C_u(t, 
       kappa = scaled_parameters['kappa'],
       rho = scaled_parameters['rho'],
       Cu_0 = 400.):
    # return the value of C_u at time t
    return 1/(kappa + np.exp(-rho*t)*(1/Cu_0 - kappa))

def C_i(t):
    # t has to be between 0 and 1
    # raise error if t is not in the range
    if t.any() < 0 or t.any() > 1:
        raise ValueError("t has to be between 0 and 1")
    # return the value of C_i at time t
    return np.zeros_like(t)

def C_l(t):
    # t has to be between 0 and 1
    if t.any() < 0 or t.any() > 1:
        raise ValueError("t has to be between 0 and 1")
    # return the value of C_l at time t
    return np.zeros_like(t)

def V(t):
    # t has to be between 0 and 1
    if t.any() < 0 or t.any() > 1:
        raise ValueError("t has to be between 0 and 1")
    # return the value of V at time t
    return np.zeros_like(t)

"""
visualize the temporal dynamics of the virus, uninfected and infected tumor cells using the fitted model from the result
from day 3 to day 13
get the simulation results for the optimized parameters
"""

from pypesto.visualize.model_fit import visualize_optimized_model_fit, _get_simulation_rdatas

amici_model = problem.objective.amici_model

species_to_plot = ['C_u', 'C_i', 'C_i1', 'C_i2', 'C_i3', 'C_i4', 'C_i5', 'C_l', 'V']

# simulate from day 3 to day 12
stop_day = 10
timepoints = np.linspace(start=0, stop=stop_day, num=50)

simulation_rdatas = _get_simulation_rdatas(
    result=result,
    problem=problem,
    start_index = 0,
    simulation_timepoints=timepoints,
)

# Plot all state trajectories
for c_ in range(len(problem.objective.edatas)):
    for species in species_to_plot:
        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        label = f"{species} - {problem.objective.edatas[c_].id}"
        if problem.objective.edatas[c_].id == 'ctrl':
            ax.plot(timepoints, simulation_rdatas[c_]['x'][:, amici_model.getStateIds().index(species)], color='red', label=label, lw = 3)
        else:
            ax.plot(timepoints, simulation_rdatas[c_]['x'][:, amici_model.getStateIds().index(species)], color='blue', label=label, lw = 3)
        ax.set_ylabel(species)
        ax.set_xlabel('Time (days)')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        # Put a legend to the right of the current axis
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.title(f"State {species} over time")
        # Change x-axis range and labels
        ax.set_xticks(np.arange(0, stop_day, 1))
        ax.set_xticklabels(np.arange(3, stop_day+3, 1))
        plt.xticks(rotation=60)
        plt.savefig(folder_path + f"{species}_over_time_{problem.objective.edatas[c_].id}.pdf", dpi=300, bbox_inches='tight')
        plt.show()

"""
visualize the temporal dynamics of the virus, uninfected and infected tumor cells using the fitted model from the result
from day 2 to day 13
"""

# obtain the dynamics of C_u, C_i, V over time
timepoints_previous = np.linspace(start=0, stop=1, num=100)
C_u_previous = C_u(timepoints_previous)
C_i_previous = C_i(timepoints_previous)
C_i1_previous = C_i(timepoints_previous)
C_i2_previous = C_i(timepoints_previous)
C_i3_previous = C_i(timepoints_previous)
C_i4_previous = C_i(timepoints_previous)
C_i5_previous = C_i(timepoints_previous)
C_l_previous = C_l(timepoints_previous)
V_previous = V(timepoints_previous)

C_u_ctrl_whole = C_u(np.concatenate([timepoints_previous, timepoints+1]))

# adjust all state trajectories with previous dynamics
for c_ in range(len(problem.objective.edatas)):
    for species in species_to_plot:
        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        label = f"{species} - {problem.objective.edatas[c_].id}"
        if problem.objective.edatas[c_].id == 'ctrl':
            ax.plot(timepoints + 1, simulation_rdatas[c_]['x'][:, amici_model.getStateIds().index(species)], color='red', label=label, lw = 3)
            if species == 'C_u':
                ax.plot(timepoints_previous, C_u_previous, color='red', lw = 3)
                ax.plot(np.concatenate([timepoints_previous, timepoints+1]), C_u_ctrl_whole, '--', color='purple', label='C_u - ctrl (analytical)', lw = 3)
            elif species == 'C_i':
                ax.plot(timepoints_previous, C_i_previous, color='red', lw = 3)
            elif species == 'C_i1':
                ax.plot(timepoints_previous, C_i1_previous, color='red', lw = 3)
            elif species == 'C_i2':
                ax.plot(timepoints_previous, C_i2_previous, color='red', lw = 3)
            elif species == 'C_i3':
                ax.plot(timepoints_previous, C_i3_previous, color='red', lw = 3)
            elif species == 'C_i4':
                ax.plot(timepoints_previous, C_i4_previous, color='red', lw = 3)
            elif species == 'C_i5':
                ax.plot(timepoints_previous, C_i5_previous, color='red', lw = 3)
            elif species == 'C_l':
                ax.plot(timepoints_previous, C_l_previous, color='red', lw = 3)
            elif species == 'V':
                ax.plot(timepoints_previous, V_previous, color='red', lw = 3)
        else:
            ax.plot(timepoints + 1, simulation_rdatas[c_]['x'][:, amici_model.getStateIds().index(species)], color='blue', label=label, lw = 3)
            if species == 'C_u':
                ax.plot(timepoints_previous, C_u_previous, color='blue', lw = 3)
            elif species == 'C_i':
                ax.plot(timepoints_previous, C_i_previous, color='blue', lw = 3)
            elif species == 'C_i1':
                ax.plot(timepoints_previous, C_i1_previous, color='blue', lw = 3)
            elif species == 'C_i2':
                ax.plot(timepoints_previous, C_i2_previous, color='blue', lw = 3)
            elif species == 'C_i3':
                ax.plot(timepoints_previous, C_i3_previous, color='blue', lw = 3)
            elif species == 'C_i4':
                ax.plot(timepoints_previous, C_i4_previous, color='blue', lw = 3)
            elif species == 'C_i5':
                ax.plot(timepoints_previous, C_i5_previous, color='blue', lw = 3)
            elif species == 'C_l':
                ax.plot(timepoints_previous, C_l_previous, color='blue', lw = 3)
            elif species == 'V':
                ax.plot(timepoints_previous, V_previous, color='blue', lw = 3)
        
        ax.set_ylabel(species)
        ax.set_xlabel('Time (days)')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.title(f"State {species} over time")
        # Change x-axis range and labels
        ax.set_xticks(np.arange(0, stop_day+1, 1))
        ax.set_xticklabels(np.arange(2, stop_day+3, 1))
        # rotate the x-axis labels
        plt.xticks(rotation=60)
        
        plt.savefig(folder_path + f"{species}_over_time_{problem.objective.edatas[c_].id}_adjusted.pdf", dpi=300, bbox_inches='tight')
        plt.show()

# get the statistics for 95% CI
cut_off_index = find_cut_off_index(result)
x_trace_within_cut_off = np.array(find_cut_off_x_trace(result, ci = 0.95))
x_trace_within_cut_off_df = pd.DataFrame(x_trace_within_cut_off, columns=problem.x_names)

x_trace_within_cut_off_NF = find_cut_off_x_trace(result, ci = 0.95, flatten=False)
x_trace_within_cut_off_NF_top100 = x_trace_within_cut_off_NF[:100]

# unpack and reshape x_trace_within_cut_off_NF_top100 into an array with the same column number as x_trace_within_cut
x_trace_within_cut_off_NF_top100_unpacked = []
for i in range(len(x_trace_within_cut_off_NF_top100)):
    for j in range(len(x_trace_within_cut_off_NF_top100[i])):
        x_trace_within_cut_off_NF_top100_unpacked.append(x_trace_within_cut_off_NF_top100[i][j])

x_trace_within_cut_off_NF_top100_unpacked_df = pd.DataFrame(x_trace_within_cut_off_NF_top100_unpacked, columns=problem.x_names)

plt.rcParams.update({'font.size': 30})

waterfall(result, size=(wid, hei))
plt.savefig(os.path.join(folder_path, 'waterfall_plot.pdf'), dpi=dpi, bbox_inches="tight")
plt.show()

pypesto.visualize.optimizer_history(result)
plt.savefig(os.path.join(folder_path, 'optimizer_history_plot.pdf'), dpi=dpi, bbox_inches="tight")
plt.show()

plt.rcParams.update({'font.size': 30})

fig, axs = plt.subplots(1, 3, figsize=(wid, hei), sharey=False, )

pypesto.visualize.parameters(result, ax = axs[0], plot_inner_parameters=False, start_indices=cut_off_index,  colors=hex_to_rgba_gradient('#A7C9F8', '#28518B', cut_off_index))
pypesto.visualize.parameters(result, ax = axs[1], plot_inner_parameters=False, start_indices=300,  colors=hex_to_rgba_gradient('#A7C9F8', '#28518B', 300))
pypesto.visualize.parameters(result, ax = axs[2], plot_inner_parameters=False, start_indices=100,  colors=hex_to_rgba_gradient('#A7C9F8', '#28518B', 100))

axs[0].set_title('95% CI', fontsize=30)
axs[1].set_title('Top 300', fontsize=30)
axs[2].set_title('Top 100', fontsize=30)

# set all the x axis, x and y labels to have fontsize 30
for j in range(3):
    axs[j].set_xticklabels(axs[j].get_xticklabels(), fontsize=20)
    axs[j].set_xlabel('Parameter Value', fontsize=30)
    axs[j].set_ylabel('Parameter', fontsize=30)

plt.suptitle('Parameter Estimation', fontsize=30)
plt.tight_layout()
plt.savefig(os.path.join(folder_path, 'parameters_plot.pdf'), dpi=dpi, bbox_inches="tight")
plt.show()

pypesto.visualize.optimization_scatter(result, start_indices=cut_off_index, size=(wid, hei))
plt.savefig(os.path.join(folder_path, 'parameters_scatter_plot_cut_off_'+str(cut_off_index)+'.pdf'), dpi=dpi, bbox_inches="tight")
plt.show()