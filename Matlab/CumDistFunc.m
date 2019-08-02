function CDF = CumDistFunc(k, params)

%% This function computes the CDF of the associated PDF.
%
% Copyright (C) 2013  Falk Heße, Vladyslav Prykhodko.
% The terms of the license agreement can be found in the README.txt.
%

if params.dim == 1 
    CDF = myCumDistFunc1D(k, params);
end

if params.dim == 2
    CDF = myCumDistFunc2D(k, params);
end

if params.dim == 3
    CDF = myCumDistFunc3D(k, params);  
end

end

%% One-dimensional fields -------------------------------------------------
%
%

function CDF = myCumDistFunc1D(x, params)

    lMax = params.X.max;
    lMin = params.X.min;
    H    = params.H;
    flag = params.func;
    
    H2   = 2.*H;

    CDFMax = CumDistFunc1D(x, [lMax H], flag);
    CDFMin = 0;
    CDF    = CDFMax;
    
    switch flag
        case 'tFracGauss'
            if (lMin > 0)
                CDFMin = CumDistFunc1D(x, [lMin H], flag);
                CDF = (lMax.^H2.*CDFMax - lMin.^H2.*CDFMin)./(lMax.^H2 - lMin.^H2);
            end
        case 'tFracExp'
            if (lMin > 0)
                CDFMin = CumDistFunc1D(x, [lMin H], flag);
                CDF = (lMax.^H2.*CDFMax - lMin.^H2.*CDFMin)./(lMax.^H2 - lMin.^H2);
            end
    end
    
end

function CDF = CumDistFunc1D(k, params, flag)

lambda  = params(1);
H       = params(2);

H2      = 2.*H;
t       = k.*lambda;

switch flag
    case 'Gauss'
        CDF = erf(t./sqrt(2));
    case 'Exp'
        CDF = 2./pi.*atan(t);
    case 'Lor'
        CDF = sign(t).*(1 - exp(-sign(t).*t));
    case 'Kolm'
        c0 = 0.001;
        for i = 1:length(k)
            if k(i) < c0
                CDF(i) = 0;
            else
                CDF(i) = 1 - (c0.^(2/3).*k(i).^(-2/3));
            end
        end
    case 'tFracGauss'
        var = t.^2./pi;
        A = erf(sqrt(var));
        si = var.^(-H)./sqrt(pi);
        B = si.*(gamma(0.5 + H) - myGammaInc(0.5 + H, var));
        CDF = sign(k).*(A - B);
        CDF = A
    case 'tFracExp'
        A = atan(t);
        B = t./(1 + H2).*hypergeomLaplace(1, 0.5 + H, 1.5 + H, -t.^2);
        CDF = 2./pi.*(A - B);
  
end

CDF(find(k == 0)) = 0;

end

%% Two-dimensional fields -------------------------------------------------
%
%

function CDF = myCumDistFunc2D(x, params)

    lMax = params.X.max;
    lMin = params.X.min;
    H     = params.H;
    flag  = params.func;

    H2    = 2.*H;

    CDFMax = CumDistFunc2D(x, [lMax H], flag);
    CDFMin = 0;
    CDF    = CDFMax;
    
    switch flag
        case 'tFracGauss'
            if (lMin > 0)
                CDFMin = CumDistFunc2D(x, [lMin H], flag);
                CDF = (lMax.^H2.*CDFMax - lMin.^H2.*CDFMin)./(lMax.^H2 - lMin.^H2);
            end
        case 'tFracExp'
            if (lMin > 0)
                CDFMin = CumDistFunc2D(x, [lMin H], flag);
                CDF = (lMax.^H2.*CDFMax - lMin.^H2.*CDFMin)./(lMax.^H2 - lMin.^H2);
            end
    end

end

function CDF = CumDistFunc2D(k, params, flag)

lambda  = params(1);
H       = params(2);

nUp     = 1./params(1);
H2      = 2.*H;
H1      = H + 1;
t       = k.*lambda;

switch flag
    case 'Gauss'
        % CDF = erf(k./(2).*lambda);
        lambda = lambda.*pi./2;
        CDF = sign(k).*(1 - exp(-(lambda*k/2).^2));
    case 'Exp'
        % not yet implemented        
    case 'tFracGauss'
        tau = t.^2./pi;
        A = exp(-tau);
        B = (gamma(H1) - myGammaInc(H1, tau))./tau.^H;
        CDF = sign(k).*(1 - A - B);
    case 'tFracExp'        
        tau = t.^2;
        si = sign(k).*gamma(1 + H)./(gamma(2 + H).*sqrt(1 + tau).*2);
        A = (2 + H2).*(sqrt(1 + tau) - 1);
        B = tau.*hypergeomLaplace(1, 0.5 + H, 2 + H, -tau);
        CDF = si.*(A - B);
end

CDF(find(k == 0)) = 0;

end

%% Three-dimensional fields -----------------------------------------------
%
%

function CDF = myCumDistFunc3D(x, params)

lMax    = params.X.max(1);
lMin    = params.X.min(1);
H       = params.H;
flag    = params.func;
    
H2      = 2.*H;

CDFMax  = CumDistFunc3D(x, [lMax H], flag);
CDFMin  = 0;
CDF     = CDFMax;
    
switch flag
    case 'tFracGauss'
        if (lMin > 0)
            CDFMin = CumDistFunc3D(x, [lMin H], flag);
            CDF = (lMax.^H2.*CDFMax - lMin.^H2.*CDFMin)./(lMax.^H2 - lMin.^H2);
        end
    case 'tFracExp'
        if (lMin > 0)
            CDFMin = CumDistFunc3D(x, [lMin H], flag);
            CDF = (lMax.^H2.*CDFMax - lMin.^H2.*CDFMin)./(lMax.^H2 - lMin.^H2);
        end
end

end

function CDF = CumDistFunc3D(k, params, flag)

lambda  = params(1);
nl      = 1./params(1);
H       = params(2);

H2      = H.*2;
H05     = H - 0.5;
H15     = H + 1.5;
t       = k.*lambda;

switch flag
    case 'Gauss'
        % CDF = -2*exp(-(k.^2)/(pi*nl^2)).*k/(nl*pi) + erf(k/(nl*sqrt(pi)));
        CDF = -2.*exp(-t.^2/pi).*t./pi + erf(t./sqrt(pi));
    case 'Exp'
        % CDF = 2*nl/pi*( atan(k/nl)/nl - (k./(nl^2 + k.^2)) );
        CDF = 2./pi.*( atan(t) - t./(1 + t.^2) );
    case 'tFracGauss'        
        vf1 = -(2*exp(-((k.^2)/(nl^2 *pi))).*k)/(pi*nl);
        vf2 = erf(k/(nl*sqrt(pi)));
        vf3 = -2*(nl^(2*H))*( pi^(-(1/2) + H))*k.^(-2 *H);
        vf4 = gamma(1.5+H)-gamma(1.5+H)*(1-gammainc((k.^2)/(pi*nl^2),1.5+H));
        CDF = vf1 + vf2 + vf3.*vf4;
        
        % A = erf(t./sqrt(pi)) - 2.*t./pi.*exp(-t.^2./pi);
        % si = - 2./t.^H2.*pi^H05;
        % B = gamma(H15) - gamma(H15).*(1 - gamma(H15).*gammainc(t.^2/pi, H15));
        
        % CDF = A - si.*B;
    case 'tFracExp'
        % A = (3 + H2).*(nl.^2).*(-k + nl.*atan(k./nl));
        B = (1 + H2).*(k.^3).*hypergeomEuler(1, 3/2 + H, 5/2 + H, -(k.^2)/(nl^2));
        si = 2/((3 + H2)*(nl^3)*pi);
        
        A = (3 + H2).*nl.^3.*(atan(t) - t);
        % B = (1 + H2).*k.^3.*hypergeomEuler(1, 3/2 + H, 5/2 + H, -t.^2);
        % si = 2/((3 + H2)*(nl^3)*pi);
        CDF = 2*((3+2*H)*(nl^2)*(-k+nl*atan(k/nl))+ (1+2*H)*(k.^3).*hypergeomEuler(1,3/2+H,5/2+H,-(k.^2)/(nl^2)))/((3+2*H)*(nl^3)*pi);
        %  = si*(A + B);
end

CDF(find(k == 0)) = 0;

end
