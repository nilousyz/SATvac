# SATvac
Stochastic Agent-based T cell Vaccine

# Description:
This is the MATLAB source code for the SATvac (Stochastic Agent-based T-cell Vaccination) model, a stochastic agent-based framework for simulating CD8+ T cell dynamics following vaccination. This model captures main immune response phases including activation, expansion, and contraction. It also tracks T cell differentiation into effector and memory cell types and explains the variability observed in immune responses by modeling stochasticity at the single cell level.

# Instructions:
1. Open the MouseDriver.m file
2. Set the experiment number by changing the EXP_N variable:
EXP_N = 1 Training set
EXP_N = 2 Testing set 1
EXP_N = 3 Testing set 2
3. The default number of simulation runs is 200. To change this, open the plot_simulation_results_cloud.m file (located in the Functions folder) and update the limit in the parfor loop.
4. Run the MouseDriver.m file to start the simulation.
