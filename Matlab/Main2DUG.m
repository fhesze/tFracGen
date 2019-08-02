function Main2DUG

%% This function generates a 2D Gaussian random field with specified
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
% * func, type of the model function (Gauss, Exp, TruncGauss, TruncExp)
%

muK         = 0.416666;                                                   
si2K        = 0.29;                                                        
muY         = log(muK) - 0.5.*log(1 + si2K./muK.^2);                       
si2Y        = log(1 + si2K./muK.^2);                                           
L_geo       = [.1 .1];                                                          
l_geo       = [0 0];                                                                                      
H           = 0.33;                                                      

func        = 'Exp';    

%%% Geometrical parameters
%
% These are the geometrical parameters, which might change according to 
% the implementation.
%
% * L_num, are the maximum numerical lenght scales (domain size)
% * nGrid, are the number of numerical grid points
% * l_num, are the according minimum numerical lenght scales (resolution)
%

nGrid       = [101 101];
L_num       = [1 1];
l_num       = L_num./(nGrid - 1);

x(1,:)      = ((1:nGrid(1)) - 1)*l_num(1);
x(2,:)      = ((1:nGrid(2)) - 1)*l_num(2);

[XI YI]     = meshgrid(x(1,:), x(2,:));
xVec        = reshape(XI, nGrid(1)*nGrid(2), 1);
yVec        = reshape(YI, nGrid(1)*nGrid(2), 1);

gamma       = 0;
% R           = [cos(gamma) -sin(gamma); sin(gamma) cos(gamma)];
% xR          = R*x;
xVecR       = xVec.*cos(-gamma) - yVec.*sin(-gamma);
yVecR       = yVec.*cos(-gamma) + xVec.*sin(-gamma);

alpha       = -pi/3;
T           = [1 tan(alpha); 0 1];
xT          = T*x;
xVecT       = xVec + yVec.*tan(alpha);
yVecT       = yVec;


%%% Numerical parameters
%
% These are the numerical parameters
%
% * nReal, is the number of realizations
% * nMode, is the number of modes used for the numerical aproximation
%

nReal       = 1000;                                                        
nMode       = 2^(10);                                                      

params.X.max= L_geo(1);
params.X.min= l_geo(1);
params.Y.max= L_geo(2);
params.Y.min= l_geo(2);
params.H    = H;
params.nMode= nMode;
params.func = func;
params.dim  = 2;

%%-- Erzeugung von Zufallsfeld u ------------------------------------------

% u = zeros(nGrid(1), nGrid(2));

for n = 1:nReal
    
    %%% computing the N(0,1) Gaussian random field 
    %
    % * irrRandom2D uses the Randomization method
    % * irrHybrid2D uses the Hybrid method
    %
    
    u = irrRandom2D(xVecR, yVecR, params);  
     
    %%% transforming the random field 
    %
    % * the log-hydraulic conductivity field Y with N(muY, si2Y)
    % * the hydraulic conductivity field K with LogN(myK, si2K)
    %
    
    Y = sqrt(si2Y).*u + (muY);                                           
    K = exp(Y);                                                           
         
    %%% postprocessing 
    %
    % * estimPDF determines the histogram of the on-point distribution
    % * estimVariogram1D determines the empirical variogram
    % 
     
    % [myPDF(n,:) myVar] = estimPDF(u, 100);
    [gamX(n,:) lag] = estimVariogram([x(1,:) x(2,:)], reshape(u, nGrid(1), nGrid(2)), 'X');                    
    gammaXAv        = mean(gamX);
    [gamY(n,:) lag] = estimVariogram([x(1,:) x(2,:)], reshape(u, nGrid(1), nGrid(2)), 'Y');                              
    gammaYAv        = mean(gamY);                                         
    
    %%% plotting and saving
    %
    % Here the plotting and saving routines, which are called in every
    % iteration, are specified.
    %

    drawnow; grid on;
    % plot(lag, gammaXAv, lag, Variogram(lag, [1 L_geo(2) 0 H], func));
    plot(lag, gammaXAv, lag, gammaYAv, lag, Variogram(lag, [1 L_geo(1) 0 H], func));
    title([ num2str(n) ' Realizations']);

    % dir_path        = 'Fourier/TruncGauss/Field_Dec0.1_No';
    % save_path       = [home_path dir_path num2str(n) '.data'];
    % save(save_path, 'u', '-ascii');
    
end

%%-- Plotting -------------------------------------------------------------

% plot3(xVec, yVec, u);
% scatter3(xVec, yVec, u);
% contour(x(1,:), x(2,:), u, 20);
contour(x, y, reshape(u, nGrid(1), nGrid(2)), 20);

end

%% Auxillary Functions ----------------------------------------------------

function u = irrRandom2D(x, y, params)

    nMode = params.nMode;
    
    Ksi = randn(nMode, 2);
    Nu = getRandomSet(params);

    for i = 1:length(x)
        u(i) = calculate_u(x(i), y(i), nMode, Ksi, Nu);
    end


end

function u = calculate_u(x, y, n0, Ksi, Nu)

    theta = (Nu(:,1).*x + Nu(:,2).*y);
    u = (1/sqrt(n0))*sum((Ksi(:,1).*cos(theta) + Ksi(:,2).*sin(theta)));
    
end
