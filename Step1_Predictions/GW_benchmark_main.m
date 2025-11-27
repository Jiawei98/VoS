clear
tic

%**************************************************************************
% Parameters Setting
%**************************************************************************

% training window list
trnwin_list = [12, 60, 120, 360];

% demean for X and Y
demeanX = 0;
demeanY = 0;

%**************************************************************************
% Get GW benchmark
%**************************************************************************

for trnwin = trnwin_list
    GW_benchmark_function(trnwin, demeanX, demeanY);
    GW_best_benchmark_function(trnwin, demeanX, demeanY);
end