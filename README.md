# Implications of a single case of paralysis for cryptic poliovirus circulation in Rockland County, New York

## Requirements
The `tidyverse` suite including `dplyr`, `tidyr` and `ggplot2`. `patchwork`, `lubridate`, `LaplacesDemon` and `ExtraDistr` are also used.

## General overview
All models fits and simulations are run using the `simulation_cluster_twoimmune.R` script. For the full suite of results, we run this script multiple times on the Harvard cluster (<https://www.rc.fas.harvard.edu/>) to generate many thousands of simulated trajectories. Cluster job submission is managed by `cluster_submit.sh`.

Briefly, the script runs the following actions:

1. Read in a table of assumed model parameter priors from `pars/priors.csv`.
2. For the high immunity and low immunity scenario analyses, run `nsims` simulations with a single call to `random_simulation_twoimmune`, each time drawing a new set of parameters from the prior (using the `simulate_priors` function).
3. Extract the simulation results and filter down to only those draws which were consistent with the observed data.
4. Using the successful draws, continue the simulation runs to generate complete trajectories. This part uses the final conditions from the first call to `random_simulation_twoimmune`, and calls `restart_simulations_table_twoimmune` to restart all the trajectories where `random_simulation_twoimmune` left off.
5. Simulate runs from the same posterior distribution using population parameters for New York City. Again, this just uses `random_simulation_twoimmune` but with different input parameters.

The simulation can instead be run locally by changing the working directory on L6/7, and by setting the starting index `i` on L17 rather.

## Detail of main simulation flow
The simulation is run with the following cascades:

1. `simulation_cluster_twoimmune.R` is the top level master script. This calls...
2. `random_simulation_twoimmune`, which manages the multiple prior draws and simulation runs and returns a large list of all the simulated trajectories and parameter values. This functions calls two other functions:
    a. `simulate_priors`, which simply returns a large table of parameter values with `n` draws from the priors. It takes many arguments, but these all just correspond to a parameter in the model as described in the methods.
    b. `run_simulation_twoimmune`, the main simulation code. This is called one time for each draw from the prior, and takes single value for each model parameter. It returns data frames for incidence and the final state of the simulation.
3. `restart_simulations_table_twoimmune` takes in the final conditions (which are just a bunch of data frames) from `random_simulation_twoimmune` and restarts the each trajectory where they left off. In the restart, there is no need to stay consistent with any observed data, and the simulation is simply run for a set duration.

## Top level simulation functions
### Initial run
`random_simulation_twoimmune` manages the top level of the simulation. 

KEY ARGUMENTS: It takes in all of the arguments corresponding to the prior distributions, the vector of observed paralysis incidence, the population size and simulation duration, and the number of simulations to be run. Each parameter with a prior must have its two parameters (mean/var or shape/scale or rate) specified. There are some optional flags to fix many of the model parameters rather than generating prior draws, but these can be ignored.

KEY OUTPUTS: a list with:

1. `simulation_results` is a merged data frame giving the incidence trajectories for _all_ of the simulated trajectories.
2. `final_conditions` is a merged data frame with the final conditions (time step, population compartment sizes etc) of each simulation. These are used to restart the simulations later.
3. `par_table` is a data frame with the prior draws.
4. `infectious_periods` is a merged data frame giving the generation interval distribution of the fully susceptible group for each simulation
5. `infectious_ps_periods` as in 4) but for the partially susceptible group
6. `incubation_periods` merged data frame giving the incubation period distribution for each simulation
7. `tmax_vector` vector, where each entry gives the number of time steps each simulation was run for
8. `n_paralysis` data frame giving the total number of paralysis cases simulated before termination for each simulation
9. `data_are_consistent` vector with an entry for each simulation. This gives the vector of flags showing if each trajectory remained consistent with the observed data or not.

### Restart run
`restart_simulations_table_twoimmune` is very similar to `random_simulation_twoimmune`, but rather than generating new prior draws, it simply uses the model parameters and final conditions of a previous call to `random_simulation_twoimmune` and continues running the trajectories for `tmax` time steps, regardless of their consistency with any observed data.

KEY ARGUMENTS: 

1. `use_sims` is a vector of the simulation IDs for which we'd like to extend the trajectories. These should all have corresponding entries in `pars` and `final_conditions`.
2. `pars` is a large data frame giving the model parameters used to generate the trajectories.
3. `final_conditions` a large data frame with all of the final conditions of the previous simulation runs.
4. `t_starts` is a vector with one entry for each simulation. This gives the time step at which the simulation should be restarted ie. the final time step from the previous call to `random_simulation_twoimmune`.
5. `vaccinate_proportion` is an optional matrix with one row for each vaccination strategy to be trialled. By default, this is a matrix with a single row dictating no vaccination.
6. `tmax` how long to run the extended simulations for in time steps.
7. `nruns` optional, but allows each extended simulation to be run multiple times (ie. with different random course)
8. `P` is the population size and should correspond to whatever was used in the initial call to `random_simulation_twoimmune`.

KEY OUTPUTS: one combined data farme giving the daily incidence of paralysis/infections for each immune class and Rt, with a simulation index.

### Main body of simulation code
`run_simulation_twoimmune` is the main simulation function. It does not call any other user-defined functions. 

KEY ARGUMENTS: It takes arguments for each of the key model parameters, including R0, the generation interval distribution parameters, the infection to paralysis fraction etc. A key argument is `observed_data`, which is simply a vector of consecutive days giving the number of observed paralysis cases. For a normal run, the `continue_run` and `restart_simulation` flags should be set to `FALSE`. These are accompanied by `NULL` for `final_conditions`, `t_start` and `vaccinate_proportion`. These are only used for special cases described later.

KEY OUTPUTS: a list with:

1. `inc_dat` A data frame of infection incidence and paralysis incidence of each of the susceptible groups, as well as Rt.
2. `final_conditions` A data frame giving the entire state of the simulation (population sizes in each immune state, incidence curves etc) at the end.
3. `t_end` the final time point of the simulation
4. `n_paralysis` the total number of paralysis cases simulated
5. `data_are_consistent` flag showing whether the simulation ended entirely consistently with `observed_data` or not.

The main body of the simulation starts at L178. This is a while loop that steps the simulation forward in daily increments. The while loop ensures that the simulation terminates after `tmax` time steps. Otherwise, there are two conditions for the simulation to continue running:

1. If `continue_run==FALSE`, the simulation starts a timer after the first case of paralysis is simulated. The simulation then continues for as many days as there are observations in `observed_data`, as long as the sequence of of observed data and the sequence of simulated paralysis cases continue to be consistent.
2. If `continue_run==TRUE`, then the simulation does NOT need to be consistent with the observed data, and the simulation will simply run for `tmax` time steps.

New infections and cases of paralysis are handled from L196 to L285. These lines are simply draws from various distributions and then bookkeeping around population sizes.

The block at L294 acts as a sort of timer -- once the first case of paralysis has had onset, each day the simulated trajectory is checked for consistency with the observed data on L298. This checks that the vector of observed data and the incidence of paralysis to date are ENTIRELY consistent.

L307 onwards simply compiles the simulation state, incidence curves etc into data frames and values to be returned.


 
