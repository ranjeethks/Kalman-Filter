%Ranjeeth KS, University of Calgary
clc;
clear all;
close all;

[InpDataStatic InpDataKinematic InpDataSimulator epochSt epochKi epochSim]= ReadData();

%%
%settings
settings.no_iterationsinLS = 10;
settings.method = 'curvi'; % Curvilinear System is implemented

%zenith error of pseudorange 0.25 m
settings.sigma02 = (0.25).^2; 

%zenith error of pseudorange rate 3cm/s as given in the lab handout
settings.sigma02prr = (0.03).^2; 

% R Matrix is weighed based on elevation
settings.weightedDOP = 1;  

%%%%%%%%%% TYPE THE FOLLOWING ID FOR THE DATA TO BE PROCESSED %%%%%%%%%%%
%       settings.dataID = 'STA' if STATIC DATA TO BE PROCESSED          %
%       settings.dataID = 'KIN' if KINEMATIC DATA TO BE PROCESSED       %
%       settings.dataID = 'SIM' if SIMULATED DATA TO BE PROCESSED       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
settings.dataID = 'SIM'; 

%%%%%%%%%% TYPE THE FOLLOWING ID FOR THE MEASUREMENT TO BE USED %%%%%%%%%%%
%   settings.measurement = 'PRO' if PSEUDO-RANGE MEASUREMENT ONLY         %
%   settings.measurement = 'PRR' if PSEUDO-RANGE + PSEUDO-RANGE RATE      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
settings.measurement = 'PRO';

%%
%computes LS with weighted DOP
%runLS=singlepoint(InpDataStatic,settings);

if (settings.dataID == 'STA')
     %Set the following variables depending on the task

     %States are posiion and clock only: X = [dxe, dxn, dxu, clk bias, clk drift]
     settings.states = 'PCO' 

     %Only pseudorange measurement is considered
     settings.measurement = 'PRO'; 

     settings.epochs = epochSt;
     settings.num_states = 5; % 5 states
     runKF_Static = KalmanFilter(InpDataStatic, settings); 

elseif (settings.dataID == 'KIN')
    %Set the following variables depending on the task

    %States are posiion, clock, and only X = [dxe, dxn, dxu, clk bias, clk drift]
    settings.states = 'PVC'; %position, velocity, and clock: [dxe, dxn, dxu, clk bias, dve, dvn, dvu, clk drift]

    settings.epochs = epochKi;
    settings.num_states = 8; % 8 states

    runKF_Kinematic = KalmanFilter(InpDataKinematic, settings); %change the input data as per the task

elseif (settings.dataID == 'SIM')
    %Set the following variables depending on the task
    %[dxe, dxn, clk bias, dve, dvn, clk drift, angular rate of rotation]
    settings.states = 'PVC'; 

    %pseudo-range with rate is used for simulated data
    settings.measurement = 'PRR'; 
    settings.epochs = epochSim;
    settings.num_states = 7; %2D position + 2D velocity + clock bias + clock drift + angular rate of rotation

    runKF_Simulation = KalmanFilterSimulation(InpDataSimulator, settings); 
 
end





analysis(run);
%%
