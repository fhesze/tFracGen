function variogram = Variogram(x, params, flag)

%% This function computes the Variogram function.
%
% Copyright (C) 2013  Falk Heße, Vladyslav Prykhodko.
% The terms of the license agreement can be found in the README.txt.
%

si    = params(1);
lMax  = params(2);
lMin  = params(3);
H     = params(4);

t     = abs(x./lMax);
kappa = H;

switch flag
    case 'Gauss'
        variogram = si.*(1 - exp(-t.^2./2) );
    case 'Exp'
        variogram = si.*(1 - exp(-t) );
    case 'Lor'
        variogram = si.*(1 - 1./(t.^2 + 1) );
    case 'Matern'
        variogram = si.*VariogramMatern(x, [lMax kappa]);
    case 'Power'
        variogram = si.*x.^(2.*H);
    case 'Kolm'
        si = -2.^(5./3).*pi.^(2./3);
        variogram = si.*gamma(-2./3).*x.^(2./3);
    case 'tFracGauss'
        variogram = si.*VariogramTrunc(x, [lMax lMin H], flag);
    case 'tFracExp'
        variogram = si.*VariogramTrunc(x, [lMax lMin H], flag);
end

end

function variogram = VariogramMatern(x, params)

lMax    = params(1);
kappa   = params(2);

t       = sqrt(2.*kappa).*x./lMax;

variogram = 1 - t.^kappa./(gamma(kappa).*2.^(kappa - 1)).*besselk(kappa,t);

variogram(find(x == 0)) = 0;

end

function variogram = VariogramTrunc(x, params, flag)

lMax = params(1);
lMin = params(2);

% parameters of the powervariogram
H = params(3);
H2 = 2.*H;

%%-- zentraler Algorithmus ------------------------------------------------

gammaMax = VariogramTruncSingle(x, [lMax H], flag);
gammaMin = VariogramTruncSingle(x, [lMin H], flag);
variogram = (lMax^H2.*gammaMax - lMin^H2.*gammaMin)./(lMax.^H2 - lMin.^H2);
    
end

function variogram = VariogramTruncSingle(x, params, flag)

lambda  = params(1);
H       = params(2);

t       = abs(x./lambda);

if lambda ~= 0  
    switch flag
        case 'tFracGauss'        
            t = t.^2.*pi./4;
        case 'tFracExp'
            H = 2.*H;
    end   
        % myGamma = gamma(1 - H).*gammainc(t, 1 - H, 'upper');
        variogram = 1 - exp(-t) + t.^(H).*myGammaInc(1 - H, t);
else
    variogram = 0;
end
        
end
