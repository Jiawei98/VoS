clear
tic

%**************************************************************************
% Parameters Setting
%**************************************************************************

% gamma in Random Fourier Features
gamma = 2;

% training window list
trnwin_list = [12, 60, 120];

% number of simulations for rff_function
nSim = 100;

% number of simulations for rff_best_function
nSim_best = 10;

% variance standardization = True
stdize = 1;

% demean for X and Y
demeanX = 0;
demeanY = 0;

%**************************************************************************
% Predictions for rff_function with 100 simulations
% Note: Parallelization or HPC can be used to expedite the for loop
%**************************************************************************

for trnwin = trnwin_list
    for random_seed = 1:nSim
        rff_function(gamma, trnwin, random_seed, stdize, demeanX, demeanY);
    end
end

%**************************************************************************
% Predictions for rff_best_function with 10 simulations
% Note: Parallelization or HPC can be used to expedite the for loop
%**************************************************************************

for trnwin = trnwin_list
    for random_seed = 1:nSim_best
        rff_best_function(gamma, trnwin, random_seed, stdize, demeanX, demeanY);
    end
end
