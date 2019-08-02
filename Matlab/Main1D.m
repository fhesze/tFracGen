function Main1D

%% This function generates a 1D Gaussian random field with specified
% correlation function and N(0,1). This field is subsequentely transformed
% into a Gaussian field with N(muY, si2Y) representing the log hydraulic
% conductivity Y and finally a log-nomal field with LogN(myK, si2K)
% representing the hydraulic conductivity K.
%
% Further details on the used methods can be found in 'Generating random 
% fields with a truncated power-law variogram. A comparison of several
% numerical methods with respect to accuracy and reproduction of structural
% features'.
%
% Copyright (C) 2013  Falk He??e, Vladyslav Prykhodko.
% The terms of the license agreement can be found in the README.txt.
%

%% Defining Parameters ----------------------------------------------------
% 
% In this part the parameters of the random field are defined. There are 
% grouped according to their use.
%

%%% Geostatistical parameters
%
% Here the following parameters are defined:
% 
% * mu_K, Expectation value of K,
% * si^2_K, Variance of K (values according to Sudicky1986, MacKay1986),
% * mu_Y, Expectation value of Y,
% * si^2_Y, Variance of Y,
% * L_geo, maximum geological length scale (correlation or cut-off length),
% * l_geo, minimum geological length scale,
% * H, Hurst coeffient (used only for fractal model functions) as well as
% * func, type of the model function (Gauss, Exp, tFracGauss, tFracExp) 
%

muK         = 0.416666;                                                   
si2K        = 0.29;                                                        
muY         = log(muK) - 0.5.*log(1 + si2K./muK.^2);                       
si2Y        = log(1 + si2K./muK.^2);                                           
L_geo       = 1;                                                          
l_geo       = 0.0;                                                                                      
H           = 0.35;  

func        = 'Gauss';

%%% Geometrical parameters
%
% These are the geometrical parameters, which might change according to 
% the implementation.
%
% * L_num, are the maximum numerical lenght scales (domain size)
% * nGrid, are the number of numerical grid points
% * l_num, are the according minimum numerical lenght scales (resolution)
%

nGrid       = 1000;
L_num       = 10.*L_geo;
l_num       = L_num./(nGrid - 1);

x           = ((1:nGrid(1)) - 1)*l_num(1);

%%% Numerical parameters
%
% These are the numerical parameters
%
% * nReal, is the number of realizations
% * nMode, is the number of modes used for the numerical aproximation
%

nReal       = 10;
nMode       = 2^(10);                                                      

%%% Assembling the parameter structure
%

params.X.max= L_geo;
params.X.min= l_geo;
params.H    = H;
params.nMode= nMode;
params.func = func;
params.dim  = 1;

%% Central Solver ---------------------------------------------------------

for n = 1:nReal

    %%% computing the N(0,1) Gaussian random field 
    %
    % * FourMeth uses the Fourier method
    % * RandMeth uses the Randomization method
    % * HybMeth uses the Hybrid method
    % * FWMeth uses the Fourier-Wavelet method
    %
    
    % u   = FourMeth(x, params);
    u   = RandMeth(x, params);
    % u  = HybMeth(x, params);
    % u   = FWMeth(x, params);
    
    %%% transforming the random field 
    %
    % * the log-hydraulic conductivity field Y with N(muY, si2Y)
    % * the hydraulic conductivity field K with LogN(myK, si2K)
    %
    
    Y   = sqrt(si2Y).*u + (muY);                                              
    K   = exp(Y);                                                            
 
    %%% postprocessing 
    %
    % * estimPDF determines the histogram of the one-point distribution
    % * estimVariogram1D determines the empirical variogram
    % 
    
    % [myPDF(n,:) myVar] = estimPDF(u, 100);
    [gam(n,:) lag] = estimVariogram(x, u, 'X');
    gammaAv = mean(gam);
        
    %kurt(n) = kurtosis(u);
    %kurtAv = mean(kurt);

    %%% plotting and saving
    %
    % Here the plotting and saving routines, which are called in every
    % iteration, are specified.
    %
    
    drawnow;
    plot(lag, gammaAv, lag, Variogram(lag, [1 L_geo l_geo H], func));
%     plot(lag, gammaAv, lag, Variogram(lag, [1 L_geo l_geo 101/2], func));
    title([ num2str(n) ' Realizations']);
    Variogram(0.5, [1 L_geo l_geo H], func)
    % save([Path '/RandomField_' num2str(n) '.data'], 'u', '-ascii');
    % save(u, num2str(n), '-ascii');

        
end

%% Postprocessing ---------------------------------------------------------
    
%     % averaging over several estimated variograms 
%     gammaAv(i,:) = mean(gam);
%     
%     drawnow;
%     plot(lag, gammaAv, lag, Variogram(lag, [si corrLen H], func));
%     title([ num2str(n) ' Realizations' ]);
    
% lag(:, 1) = [];
% gammaAv(:,1) = [];
% 
% lag = reshape(lag, 1, nDec*length(lag(1,:)));
% gammaAv = reshape(gammaAv, 1, nDec*length(gammaAv(1,:)));
% 
% [lag, index] = sort(lag);
% gammaAv = gammaAv(index);

% gammaAv(find(lag = ))

% averaging over several estimated covariance functions 
% covAv = mean(cov);
% averaging over several estimated random corralation functions 
% myPDFAv = mean(myPDF);
% estimating the one-point PDF
% [myPDF myVar] = estimPDF(u, 200);
% fitting of the estimated variograms
% estim = VariogramFit(lag, gammaAv, flag)

%% Plotting ---------------------------------------------------------------

figure
hold on;
grid on;
 
% plot(x, u)
plot(x, CumDistFunc(x, params))
% plot(lag, gammaAv);
% plot(lag, gammaAv, lag, Variogram(lag, [si corrLen H], flag));
% plot(lag, gammaAv, lag, Variogram(lag, [si pi.*corrLen H], flag), lag, Variogram(lag, estim, flag));

end
