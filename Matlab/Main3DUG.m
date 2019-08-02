function Main3DUG

%% This function generates a 3D Gaussian random field with specified
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
L_geo       = [1 1 1];  
l_geo       = [0 0 0];
H           = 0.33;                                                        
func        = 'tFracGauss';                                                                 

%%% Geometrical parameters
%
% These are the geometrical parameters, which might change according to 
% the implementation.
%
% * L_num, are the maximum numerical lenght scales
% * nGrid, are the number of numerical grid points
% * l_num, are the according minimum numerical lenght scales
%

nGrid       = [101 101 101];
L_num       = [10 10 10].*L_geo;
l_num       = L_num./(nGrid - 1);

x(1,:)      = ((1:nGrid(1)) - 1)*l_num(1);
x(2,:)      = ((1:nGrid(2)) - 1)*l_num(2);
x(3,:)      = ((1:nGrid(3)) - 1)*l_num(3);

%%% Numerical parameters
%
% These are the numerical parameters
%
% * nReal, is the number of realizations
% * nMode, is the number of modes used for the numerical aproximation
%

nReal       = 1000;                                                        
nMode       = 2^(10);   

params.X.max= L_geo;
params.X.min= l_geo;
params.H    = H;
params.nMode= nMode;
params.func = func;
params.dim  = 3;

%% Central Solver ---------------------------------------------------------

for n = 1:nReal
    
 
    %%% computing the N(0,1) Gaussian random field 
    %
    % * irrRandom3D uses the Randomization method
    %
    
    u = irrRandom3D(x, params);

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
    [gamX(n,:) lag]  = estimVariogram(x, u, 'X');
    gammaXAv         = mean(gamX);   
    [gamY(n,:) lag]  = estimVariogram(x, u, 'Y');
    gammaYAv         = mean(gamY);   
    [gamZ(n,:) lag]  = estimVariogram(x, u, 'Z');
    gammaZAv         = mean(gamZ);
    % estim           = VariogramFit(lag, gammaAv, func)
        
    %%% plotting and saving
    %
    % Here the plotting and saving routines, which are called in every
    % iteration, are specified.
    %
    
    drawnow;
    plot(lag, gammaXAv, lag, gammaYAv, lag, gammaZAv, lag, Variogram(lag, [1 L_geo(1) 0 H], func));
    title([ num2str(n) ' Realizations']);
    
end

%% Postprocessing ---------------------------------------------------------

% estimating the one-point PDF
% [myPDF myVar] = estimPDF(Y, 200);

%% Plotting ---------------------------------------------------------------

data = smooth3(u,'box',5);
patch(isocaps(data,.5),'FaceColor','interp','EdgeColor','none');
p1 = patch(isosurface(data,.5),'FaceColor','blue','EdgeColor','none');
isonormals(data,p1);
view(3); 
axis vis3d tight
camlight left; 
lighting phong


% plot(lag, gammaAv, lag, Variogram(lag, [si2 L_geoX H], func));
% plot(lag, gammaAv);
% plot(lag, gammaAv, lag, VariogramGauss(lag, [1 L_geo]), lag, VariogramExp(lag, [1 L_geo]));
% plot(myVar, myPDF, myVar, 1./sqrt(2.*pi).*exp(-(myVar).^2./(2.^2)));

end

%% Auxillary Functions ----------------------------------------------------

function u = irrRandom3D(x, params)

    nMode = params.nMode;
    
    Ksi = randn(nMode, 2);
    Nu = getRandomSet(params);

    for i = 1:length(x(1,:))
        for j = 1:length(x(2,:))
            for k = 1:length(x(3,:))
                u(i,j,k) = calculate_u(x(1,i), x(2,j), x(3,k), nMode, Ksi, Nu);
            end
        end
    end

end

function u = calculate_u(x, y, z, n0, Ksi, Nu) 

    theta = Nu(:,1).*x + Nu(:,2).*y + Nu(:,3).*z;
    u = (1/sqrt(n0))*sum((Ksi(:,1).*cos(theta) + Ksi(:,2).*sin(theta)));
    
end
