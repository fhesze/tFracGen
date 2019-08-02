function u = RandMeth(x, params)

%% This function computes a Gaussian random field according to the
%% Randomization method.
%
% Copyright (C) 2013  Falk Heße, Vladyslav Prykhodko.
% The terms of the license agreement can be found in the README.txt.
%

if params.dim == 1
    [u] = Random1D(x(1,:), params);
end

if params.dim == 2
    [u] = Random2D(x(1,:), x(2,:), params);
end

if params.dim == 3
    [u] = Random3D(x(1,:), x(2,:), x(3,:), params);
end


end

%% One-dimensional fields -------------------------------------------------
%
%

function u = Random1D(x, params)

flag    = params.func;
lMax 	= params.X.max;
lMin 	= params.X.min;
H       = params.H;
n0      = params.nMode;

% predefined for efficiency
u 	= zeros(1, length(x));            % solution
theta 	= zeros(n0, length(x));       % all n0 weighted positions

%%-- Erzeugung von Zufallsfeld u ------------------------------------------

% gaussion-distributed random variable for the amplitudes
ksi1 	= repmat(randn(n0, 1),1,length(x));
ksi2 	= repmat(randn(n0, 1),1,length(x));
    
% function-specific random variable for sampling the Fourier domain
nu 	= getRandomSet(params);

% assembling the solution
u 	= 1./sqrt(n0).*sum(ksi1.*cos(nu'*x) + ksi2.*sin(nu'*x));

end

%% Two-dimensional fields -------------------------------------------------
%
%

function u = Random2D(x, y, params)


[XI, YI]    = meshgrid(x, y);
                                                        
u    	= zeros(size(XI, 1), size(XI,2));
n0   	= params.nMode;

Ksi         = randn(n0, 2);                                                % Gaussion random variable for sampling
Nu          = getRandomSet(params);                                        % user-specific random variable

for j = 1:n0
    
    theta   = (Nu(j,1).*XI + Nu(j,2).*YI);
    A       = Ksi(j,1)*cos(theta);
    B       = Ksi(j,2)*sin(theta);
    u       = u + sqrt(1/n0).*(A + B);
        
end
    
end


%% Three-dimensional fields -----------------------------------------------
%
%

function u = Random3D(x, y, z, params)

nGrid(1) = length(x);
nGrid(2) = length(y);
nGrid(3) = length(z);

nMode = params.nMode;

    % get random numbers
    % gaussion-distributed random variable for sampling
    Ksi = randn(nMode, 2);
    
    % uniformly-distributed random variable for the correlation function
    % Nu = getRandomSet3D([L_geo H nMode], func);
    Nu = getRandomSet(params);
    
    for k = 1:nGrid(3)
    	u(:,:,k) = Rand3D(x, y, z(k), params, Ksi, Nu);
    end
    
end

function u = Rand3D(x, y, z, params, Ksi, Nu)

n0 = params.nMode;  

[XI, YI]    = meshgrid(x, y);
u = zeros(size(XI, 1), size(XI,2));

for j = 1:n0
    
    theta = (Nu(j,1).*XI + Nu(j,2).*YI + Nu(j,3).*z);
    A = Ksi(j,1)*cos(theta);
    B = Ksi(j,2)*sin(theta);
    u = u + (1/sqrt(n0))*(A + B);
        
end
    
end
