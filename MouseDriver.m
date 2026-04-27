clear all;
close all;
clc

%% Setup
addpath("Functions")
steptime = 0.10; % 6 minute step times
stoptime = 24 * 27; 
Exp_N = 1; % Can be 1,2,3, referring to training == 1 and testing sets == 2,3
%X= [NA0/NP0 PAPbind    PE-Ileave TArrival PBPleave NA0   PJ2Jmean  PJ2JSD   PJLeaving  PC-ISplitmean PC-ISplitSD
X = [0.96	 0.00275	0.00021	  14.47	   0.104	2800  15.78	    4.808	 0.0100	    13	          9];

%% Simulation
tic
[dat,total] = MouseSim(round(X(6)),round(X(6)/X(1)),X(2),X(4),X(5),X(3), stoptime, steptime , X(7) , X(8) , X(9) , X(10) , X(11));

%% Results Demonstration
plot_simulation_results_cloud(X,Exp_N,stoptime,steptime)
rt=toc;
fprintf('Completed in %1.2f seconds.\n', rt);