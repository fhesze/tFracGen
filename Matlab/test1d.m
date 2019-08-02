function test1d

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
% Copyright (C) 2013  Falk Heße, Vladyslav Prykhodko.
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

func        = 'tFracExp';

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
L_num       = 100.*L_geo;
l_num       = L_num./(nGrid - 1);

dx          = 0.01; 
x           = [0:dx:10];

%%% Numerical parameters
%
% These are the numerical parameters
%
% * nReal, is the number of realizations
% * nMode, is the number of modes used for the numerical aproximation
%

nReal       = 1;
nMode       = 2^(10); 

k           = 10;

%%% Assembling the parameter structure
%

params.X.max= L_geo;
params.X.min= l_geo;
params.H    = H;
params.nMode= nMode;
params.func = func;
params.dim  = 1;

%% Central Solver ---------------------------------------------------------

S = SpecDens(x, params);
CDF_estim = cumsum(2*S*dx);
CDF = CumDistFunc(x, params);

% params.func = 'tFracExp';
params.H = 0.25;
S3 = SpecDens(x, params);

S2 = k.*g(x, params);

%% Plotting ---------------------------------------------------------------

hold on;
grid on;

% plot(log(x), log(S1), log(x), log(S2), log(x), log(S3));
% plot(x, CDF);
plot(x, CDF, x, CDF_estim);
% plot(lag, gammaAv, lag, Variogram(lag, [si corrLen H], flag));
% plot(lag, gammaAv, lag, Variogram(lag, [si pi.*corrLen H], flag), lag, Variogram(lag, estim, flag));

end

function y = g(x, params)
    
    lMaxX   = params.X.max(1);
    nl      = 10./lMaxX;
    
    y = nl./(pi*(nl.^2 + x.^2));

end

function y = f(x, params)

flag= params.func;
num = params.nMode;
H   = params.H;
nl  = 1;

x(find(x == 0)) = 10.^(-6);

    switch flag
        case 'Gauss'
            y = 1./sqrt(2*pi)*exp(-0.5.*(x.^2/2));
            % y = SpecDens(x, params);
        case 'Exp'
            y = 1./(pi*(1+x.^2));
        case 'tFracGauss'
            y = (H*nl^(2*H)*pi^(H-0.5))./((x.^2).^(1/2+H)).*(gamma(0.5+H)-gamma(0.5+H)*(1 - gammainc(x.^2/(pi*nl^2), 0.5+H)));
        case 'tFracExp'
            y = 2*H/(pi*(1+2*H)*nl)*hypergeomLaplace(1, 0.5+H, 1.5+H, -x.^2/nl^2);
    end        
end
