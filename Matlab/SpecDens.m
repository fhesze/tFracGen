function u = SpecDens(k, params)

%% This function computes a spectral density of a Gaussian random field.
%
% Copyright (C) 2013  Falk Heße, Vladyslav Prykhodko.
% The terms of the license agreement can be found in the README.txt.
%

if params.dim == 1
    [u] = mySpecDens1D(k, params);
end

if params.dim == 2
    [u] = SpecDens2D(k(1,:), k(2,:), params);
end

if params.dim == 3
    [u] = SpecDens3D(k(1,:), k(2,:), k(3,:), params);
end


end

%% One-dimensional fields -------------------------------------------------
%
%

function S = mySpecDens1D(k, params)

    lMax = params.X.max;
    lMin = params.X.min;
    H    = params.H;
    flag = params.func;
    
    H2   = 2.*H;

    SMax = SpecDens1D(k, [lMax H], flag);
    SMin = 0;
    S    = SMax;
    
    switch flag
        case 'tFracGauss'
            if (lMin > 0)
                SMin = SpecDens1D(k, [lMin H], flag);
                S = (lMax.^H2.*SMax - lMin.^H2.*SMin)./(lMax.^H2 - lMin.^H2);
            end
        case 'tFracExp'
            if (lMin > 0)
                SMin = SpecDens1D(k, [lMin H], flag);
                S = (lMax.^H2.*SMax - lMin.^H2.*SMin)./(lMax.^H2 - lMin.^H2);
            end
    end
    
end

function S = SpecDens1D(k, params, flag)

lMaxX  = params(1);
nl      = 1./lMaxX;
H       = params(2);

H05     = 1/2 + H;
H2      = 2.*H;

k(find(k == 0)) = 10.^(-6);
        
    switch flag
        case 'Gauss'
            S = lMaxX./sqrt(2*pi)*exp(-((k.*lMaxX).^2/2));
        case 'Exp'
            S = nl./(pi*(nl.^2 + k.^2));            
        case 'Kolm'
            S = abs(k).^(-5/3);
        case 'tFracGauss'
            si = (H*nl^(2*H)*pi^(H-0.5))./(k.^(1+2*H));
            S = si.*(gamma(0.5+H) - gamma(0.5+H)*(1 - gammainc(k.^2/(pi*nl^2), 0.5+H)));
        case 'tFracExp'
            si = 2*H/(pi*(1+2*H)*nl);
            S = si.*hypergeomLaplace(1, 0.5+H, 1.5+H, -k.^2/nl^2);
    end 

end

%% Two-dimensional fields -------------------------------------------------
%
%
% 
% function S = mySpecDens2D(x, params)
% 
%     lMaxX = params.X.max;
%     lMaxY = params.Y.max;
%     lMinX = params.X.min;
%     lMinY = params.Y.min;
%     H     = params.H;
%     flag  = params.func;
%     
%     H2    = 2.*H;
% 
%     SMax  = SpecDens1D(x, [lMaxX lMaxY H], flag);
%     SMin  = 0;
%     S     = SMax;
%     
%     switch flag
%         case 'tFracGauss'
%             if (lMin > 0)
%                 SMin = SpecDens1D(x, [lMinX lMinY H], flag);
%                 S = (lMaxX.^H2.*SMax - lMinX.^H2.*SMin)./(lMaxX.^H2 - lMinX.^H2);
%             end
%         case 'tFracExp'
%             if (lMin > 0)
%                 SMin = SpecDens1D(x, [lMinX lMinY H], flag);
%                 S = (lMaxX.^H2.*SMax - lMinX.^H2.*SMin)./(lMaxX.^H2 - lMinX.^H2);
%             end
%     end
%     
% end
% 
% function S = SpecDens2D(k, params)
% 
% 
% end

%% Three-dimensional fields -----------------------------------------------
%
%

function S = mySpecDens2D(x, params)

    lMaxX = params.X.max(1);
    lMaxY = params.X.max(2);
    lMaxZ = params.X.max(3);
    lMinX = params.X.min(1);
    lMinY = params.X.min(2);
    lMinZ = params.X.min(3);
    H     = params.H;
    flag  = params.func;
    
    H2    = 2.*H;

    SMax  = SpecDens1D(x, [lMaxX lMaxY lMaxZ H], flag);
    SMin  = 0;
    S     = SMax;
    
    switch flag
        case 'tFracGauss'
            if (lMin > 0)
                SMin = SpecDens1D(x, [lMinX lMinY lMinZ H], flag);
                S = (lMaxX.^H2.*SMax - lMinX.^H2.*SMin)./(lMaxX.^H2 - lMinX.^H2);
            end
        case 'tFracExp'
            if (lMin > 0)
                SMin = SpecDens1D(x, [lMinX lMinY lMinZ H], flag);
                S = (lMaxX.^H2.*SMax - lMinX.^H2.*SMin)./(lMaxX.^H2 - lMinX.^H2);
            end
    end
    
end

function S = SpecDens3D(kx, ky, kz, params, flag)

lMaxX = params(1);
lMaxY = params(2);
lMaxZ = params(3);
nX    = 1./lMaxX;
nY    = 1./lMaxY;
nZ    = 1./lMaxZ;
H     = params(4);

H15   = 3/2 + H;
H2    = 2.*H;

kx(find(kx == 0)) = 10.^(-6);
ky(find(ky == 0)) = 10.^(-6);
kz(find(kz == 0)) = 10.^(-6);

switch flag
    case 'TruncGauss'
        var = 1./nX.^2./pi.*(kx.^2 + ky.^2 + kz.^2);

        si = 8.*H./((nX.*pi).^3.*var.^(3/2 + H));
        S = si.*gamma(3/2 + H).*gammainc(var, 3/2 + H, 'lower');
    case 'Gauss'
        % S = lambda.^2./pi.*exp(-((kXI).^2 + (kYI).^2).*(lambda./2).^2);
        % S = 0.5./sqrt(pi).*exp(-kx.^2./(lambda.*2).^2);
        % S = CovGauss(x, [0.5./sqrt(pi) lambda.*2]);
 
end

end
