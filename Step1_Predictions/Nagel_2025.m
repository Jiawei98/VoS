%**************************************************************************
% Load Data
%**************************************************************************
clear;
stdizex = 1;
stdizey = 1;
demean = 0;
load GYdata

% Add lag return Y_{t-1} as predictor
X = [X lagmatrix(Y,1)];

% Volatility standardize
if stdizex == 1
    X = volstdbwd(X,[]);
end
if stdizey == 1
    Y2 = 0;
    for j = 1:12
        Y2 = Y2 + lagmatrix(Y.^2, j);
    end
    Y2 = Y2 / 12;
    Y  = Y ./ sqrt(Y2);
end

Y     = Y(37:end);
X     = X(37:end,:);
dates = dates(37:end,:);

Tfull = length(Y);
X     = X';
Y     = Y';
d     = size(X,1);

%**************************************************************************
% PARAMETERS
%**************************************************************************
Tlag = 12;    
gamma = 2;    
weights_EW = ones(1, Tlag) / Tlag;               % Equal weights
weights_DW = (12 - (0:11)) / 78;                 % Fixed decreasing weights

%**************************************************************************
% STORAGE
%**************************************************************************
yhat_kernel = nan(1, Tfull);
yhat_DW     = nan(1, Tfull);
yhat_EW     = nan(1, Tfull);
yhat_PV     = nan(1, Tfull);
yhat_EWPV   = nan(1, Tfull);
yhat_DWPV   = nan(1, Tfull);
yhat_EWSD   = nan(1, Tfull);   

save_path = './rff_SeparateSims/';
if ~exist(save_path, 'dir'), mkdir(save_path); end

%**************************************************************************
% MAIN LOOP
%**************************************************************************
for t = (Tlag+1):Tfull
    
    %---------------------------------------------------
    % Rolling window: demean and standardize
    %---------------------------------------------------
    trnloc  = (t-Tlag):t-1;
    Ztrn    = X(:, trnloc);
    Ytrn    = Y(trnloc);
    Ztst    = X(:, t);
    
    if demean == 1
        Ymn = nanmean(Ytrn);
        Zmn = nanmean(Ztrn, 2);
    else
        Ymn = 0;
        Zmn = 0;
    end
    
    Ytrn = Ytrn - Ymn;
    Ztrn = Ztrn - Zmn;
    Ztst = Ztst - Zmn;
    
    %---------------------------------------------------
    % Kernel vector and matrix
    %---------------------------------------------------
    kvec = zeros(Tlag,1);
    Kmat = zeros(Tlag, Tlag);
    Ztrn_cols = Ztrn(:, end:-1:end-Tlag+1); 
    
    for k = 1:Tlag
        dx = Ztst - Ztrn_cols(:, k);
        kvec(k) = exp( -0.5 * gamma^2 * (dx' * dx) );
        
        for j = 1:Tlag
            dx2 = Ztrn_cols(:, k) - Ztrn_cols(:, j);
            Kmat(k,j) = exp( -0.5 * gamma^2 * (dx2' * dx2) );
        end
    end
    
    %---------------------------------------------------
    % Build Yvec from Ytrn
    %---------------------------------------------------
    Yvec = Ytrn(end:-1:end-Tlag+1);
    
    %---------------------------------------------------
    % Kernel OLS forecast
    %---------------------------------------------------
    weights_kernel = (kvec' / Kmat);
    yhat_kernel(t) =  weights_kernel * Yvec' + Ymn;
    
    %---------------------------------------------------
    % DW forecast (decreasing weights)
    %---------------------------------------------------
    yhat_DW(t) = weights_DW * Yvec' + Ymn;
    
    %---------------------------------------------------
    % EW forecast (equal weights)
    %---------------------------------------------------
    yhat_EW(t) = weights_EW * Yvec' + Ymn;
    
    %---------------------------------------------------
    % PV forecast (predictor volatility only)
    %---------------------------------------------------
    sigma2_x = mean(var(Ztrn_cols, 0, 2));   % avg cross-predictor variance
    yhat_PV(t) = 1 / sigma2_x;
    
    %---------------------------------------------------
    % EWPV forecast (EW scaled by 1/σ²_x,t)
    %---------------------------------------------------
    yhat_EWPV(t) = (1 / sigma2_x) * (weights_EW * Yvec') + Ymn;

    %---------------------------------------------------
    % DWPV forecast (DW scaled by 1/σ²_x,t)
    %---------------------------------------------------
    yhat_DWPV(t) = (1 / sigma2_x) * (weights_DW * Yvec') + Ymn;

    %---------------------------------------------------
    % EWSD forecast: EW scaled by 1/Var(Ytrn)
    %---------------------------------------------------
    sigma_y = nanstd(Ytrn, 0);              % rolling variance of Y (demeaned window)
    yhat_EWSD(t) = (1 / sigma_y) * (weights_EW * Yvec') + Ymn;
end

%**************************************************************************
% PREPARE PLOTTING
%**************************************************************************
plot_dates = dates(Tlag+1:end);
years  = floor(plot_dates / 100);
months = mod(plot_dates, 100);
plot_datenum = datenum(years, months, 1);

plot_yhat_kernel = yhat_kernel(Tlag+1:end);
plot_yhat_EW     = yhat_EW(Tlag+1:end);
plot_yhat_DW     = yhat_DW(Tlag+1:end);
plot_yhat_PV     = yhat_PV(Tlag+1:end);
plot_yhat_EWPV   = yhat_EWPV(Tlag+1:end);
plot_yhat_DWPV   = yhat_DWPV(Tlag+1:end);
plot_yhat_EWSD   = yhat_EWSD(Tlag+1:end);

%**************************************************************************
% PLOT: Forecast Comparison
%**************************************************************************
figure;
plot(plot_datenum, plot_yhat_kernel, 'LineWidth', 2); hold on;
plot(plot_datenum, plot_yhat_EW,     'LineWidth', 2);
plot(plot_datenum, plot_yhat_DW,     'LineWidth', 2);
plot(plot_datenum, plot_yhat_PV,     'LineWidth', 2);
plot(plot_datenum, plot_yhat_EWPV,   'LineWidth', 2);
plot(plot_datenum, plot_yhat_DWPV,   'LineWidth', 2);
plot(plot_datenum, plot_yhat_EWSD,   'LineWidth', 2);
datetick('x','yyyy');
xlabel('Year');
ylabel('Forecasted Return');
legend('Kernel OLS','EW','DW','PV','EW/PV','DW/PV','EWSD','Location','best');
title('Forecast Comparison');
grid on;

%**************************************************************************
% Correlation matrix of forecasts
%**************************************************************************
forecast_names = {'Kernel', 'EW', 'DW', 'PV', 'EWPV', 'DWPV', 'EWSD'};

forecast_mat = [plot_yhat_kernel(:), plot_yhat_EW(:), plot_yhat_DW(:), ...
                plot_yhat_PV(:), plot_yhat_EWPV(:), plot_yhat_DWPV(:), ...
                plot_yhat_EWSD(:)];
corr_mat = corr(forecast_mat);

corr_table = array2table(corr_mat, ...
    'VariableNames', forecast_names, ...
    'RowNames', forecast_names);

disp('Correlation matrix of forecasts:');
disp(corr_table);

%**************************************************************************
% Correlation matrix of timing strategy
%**************************************************************************
% Ensure vectors are column-aligned
Y_plot = Y(Tlag+1:end)';  % 1 × T
Y_plot = Y_plot(:);       % T × 1

% Create timing-adjusted return series
timing_Kernel = Y_plot .* plot_yhat_kernel(:);
timing_EW     = Y_plot .* plot_yhat_EW(:);
timing_DW     = Y_plot .* plot_yhat_DW(:);
timing_PV     = Y_plot .* plot_yhat_PV(:);
timing_EWPV   = Y_plot .* plot_yhat_EWPV(:);
timing_DWPV   = Y_plot .* plot_yhat_DWPV(:);
timing_EWSD   = Y_plot .* plot_yhat_EWSD(:);

% Stack into matrix
timing_mat = [timing_Kernel, timing_EW, timing_DW, timing_PV, timing_EWPV, timing_DWPV, timing_EWSD];

% Compute correlation while ignoring rows with NaNs
corr_timing = corr(timing_mat, 'rows', 'pairwise');

% Label the matrix
forecast_names = {'Kernel', 'EW', 'DW', 'PV', 'EWPV', 'DWPV', 'EWSD'};
corr_timing_table = array2table(corr_timing, ...
    'VariableNames', forecast_names, ...
    'RowNames', forecast_names);

% Display result
disp('Correlation matrix of timing-adjusted returns (pairwise NaN removal):');
disp(corr_timing_table);

%**************************************************************************
% SAVE
%**************************************************************************
save([save_path '/Nagel_' num2str(Tlag) '.mat'], ...
     'yhat_kernel', 'yhat_EW', 'yhat_DW', 'yhat_PV', 'yhat_EWPV', 'yhat_DWPV', 'yhat_EWSD');