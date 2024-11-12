import numpy as np
import pandas as pd
import scipy
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import os
import h5py
import copy
import sys

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

plt.rcParams.update({'font.size': 30})
dpi = 100
wid = int(2560/dpi)
hei = int(1600/dpi)

def main():    

    n_runs, maxiter, load = 5000, 5000, False

    n_procs = int(os.environ.get("SLURM_CPUS_ON_NODE", os.cpu_count()))

    # optimization
    hierarchical = False

    petab_yaml = 'petab_files/dividing_infected_cells.yaml'
    petab.validate(petab_yaml)
    petab_problem = petab.Problem.from_yaml(petab_yaml)

    np.random.seed(500)

    problem = pypesto.petab.PetabImporter(
            petab_problem,
            hierarchical=hierarchical,
            model_name=f"DIVIDING_INFECTED_CELLS_Model",
        ).create_problem(force_compile=True)
    
    problem.objective.amici_model.setAllStatesNonNegative()

    # some model properties
    print("Model parameters:", list(problem.objective.amici_model.getParameterIds()), "\n")
    print("Model const parameters:", list(problem.objective.amici_model.getFixedParameterIds()), "\n")
    print("Model outputs:   ", list(problem.objective.amici_model.getObservableIds()), "\n")
    print("Model states:    ", list(problem.objective.amici_model.getStateIds()), "\n")

    if load:
        result = pypesto.store.read_result('optimization_history/n'+ str(n_runs) +'_diff_scale_v2.hdf5')
    else:
        result = minimize(
            problem=problem,
            optimizer=FidesOptimizer(verbose=False, options={'maxiter': maxiter}),
            n_starts=n_runs,
            engine=MultiProcessEngine(n_procs=n_procs),
            # startpoint_method=uniform,
            history_options = pypesto.HistoryOptions(trace_record=True, storage_file='optimization_history/n'+ str(n_runs) +'_diff_scale_v2.hdf5'),
            filename='optimization_history/n'+ str(n_runs) +'_diff_scale_v2.hdf5',
        )

    # print result summary
    print(result.summary())

if __name__ == "__main__":
    main()
