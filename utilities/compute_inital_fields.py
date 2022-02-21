"""
Usage: Compute 1D laminar premixed flame
"""

import os
import numpy as np
import cantera as ct

if __name__ == '__main__':
    print(ct.__version__)

    ##################################
    # User-defined input - Fuel is CH4
    ##################################
    pgas       =   1*ct.one_atm 
    Tin        =   900.0
    phi        =   0.5
    mixture    =   'O2:0.21, N2:0.79'
    dir_out = 'out_states'
    ###################################

    if not os.path.exists(dir_out):
        os.makedirs(dir_out)

    #Import gas phases with mixture transport model
    gas = ct.Solution('gri30.xml')
    #Set gas state to that of the unburned gas
    gas.TPX = Tin, pgas, mixture
    gas.set_equivalence_ratio(phi=phi, fuel='CH4:1', oxidizer=mixture)
    gas()
    print( "Phi = " + str(gas.get_equivalence_ratio()) )

    ##########################################
    # Create the free laminar premixed flame #
    ##########################################

    # Solver settings
    initial_grid = 2*np.array([0.0, 0.001, 0.01, 0.02, 0.029, 0.03],'d')/3
    tol_ss    = [1.0e-5, 1.0e-8] # [rtol atol] for steady-state problem
    tol_ts    = [1.0e-5, 1.0e-8] # [rtol atol] for time stepping
    loglevel  = 1
    refine_grid = True

    f = ct.FreeFlame(gas, initial_grid)
    f.flame.set_steady_tolerances(default=tol_ss)
    f.flame.set_transient_tolerances(default=tol_ts)
    f.inlet.X = gas.X
    f.inlet.T = Tin

    # First solve - energy off
    f.energy_enabled = False
    f.set_refine_criteria(ratio=7.0, slope=1, curve=1)
    # Max number of times the Jacobian will be used before re-evaluation
    f.set_max_jac_age(50, 50)
    #Set time steps whenever Newton convergence fails
    f.set_time_step(5.e-06, [10, 20, 80]) #s
    f.transport_model = 'UnityLewis' # Mix, Multi, UnityLewis
    #Calculation
    f.solve(loglevel, refine_grid)

    # Second solve - energy on
    f.energy_enabled = True
    # New refinement criteria
    f.set_refine_criteria(ratio=5.0, slope=0.5, curve = 0.5)
    # Transport model
    # f.transport_model = 'Multi'
    f.solve(loglevel, refine_grid)

    #Third solve:
    f.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.05)
    f.solve(loglevel, refine_grid)

    #Fourth solve:
    f.set_refine_criteria(ratio=2.0, slope=0.02, curve=0.02, prune=0.01)
    f.solve(loglevel, refine_grid)

    grid = f.flame.grid
    n_points = f.flame.n_points
    T = np.zeros(n_points)
    for n in np.arange(n_points):
        T[n]= f.T[n]
    np.savetxt(dir_out+'/T', T)
    np.savetxt(dir_out+'/grid', grid)

    n_species = gas.n_species
    for i in np.arange(n_species):
        species_name = gas.species_name(i)
        if species_name in ['CH4', 'O2', 'N2', 'CO2', 'H2O']:
            print('species: ' + str(species_name))
            Yi = np.zeros(n_points)
            for n in np.arange(n_points):
                f.set_gas_state(n)
                Yi[n]= gas.Y[i]
            np.savetxt(dir_out + '/' + str(species_name), Yi)
