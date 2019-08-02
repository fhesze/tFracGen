function u = getRandomSet(params)

%% This function draws a set radom deviates from the specified PDF.
%
% Copyright (C) 2013  Falk Heße, Vladyslav Prykhodko.
% The terms of the license agreement can be found in the README.txt.
%

if params.dim == 1
    [u] = getRandomSet1D(params);
end

if params.dim == 2
    [u] = getRandomSet2D(params);
end

if params.dim == 3
    [u] = getRandomSet3D(params);
end


end

%% One-dimensional fields -------------------------------------------------
%
%

function Nu = getRandomSet1D(params)

% Nu = getInverse(params); 
Nu = getReject(params); 
    
end

%% Two-dimensional fields -------------------------------------------------
%
%

function Nu = getRandomSet2D(params)

lMaxX   = params.X.max(1);
lMaxY   = params.Y.max(1);

Gamma   = rand(params.nMode, 2);
Nu      = zeros(params.nMode, 2);
    
switch params.func
    case 'Gauss'
%         r = randn(params.nMode, 1);
%         Nu(:,1) = r.*cos(2.*pi.*Gamma(:,1))./lMaxX;
%         Nu(:,2) = r.*sin(2.*pi.*Gamma(:,1))./lMaxY;
        Nu(:,1) = randn(params.nMode, 1)./lMaxX;
        Nu(:,2) = randn(params.nMode, 1)./lMaxY;
    case 'Exp'
        r = sqrt(1./Gamma(:,2).^2 - 1);
        Nu(:,1) = r.*cos(2.*pi.*Gamma(:,1))./lMaxX;
        Nu(:,2) = r.*sin(2.*pi.*Gamma(:,1))./lMaxY;
    case 'tFracGauss'
        r = getInverse(params);
        Nu(:,1) = r.*cos(2.*pi.*Gamma(:,1)')./lMaxX;
        Nu(:,2) = r.*sin(2.*pi.*Gamma(:,1)')./lMaxY;
    case 'tFracExp'
        r = getInverse(params);
        Nu(:,1) = r.*cos(2.*pi.*Gamma(:,1)')./lMaxX;
        Nu(:,2) = r.*sin(2.*pi.*Gamma(:,1)')./lMaxY;
end
  
end

%% Three-dimensional fields -----------------------------------------------
%
%

function Nu = getRandomSet3D(params)

lMax    = params.X.max;
H       = params.H;
n0      = params.nMode;  

Gamma   = rand(n0, 3);
Nu      = zeros(n0, 3);
r       = zeros(n0, 1);

switch params.func
    case 'Gauss'
        Nu(:,1) = randn(n0, 1)./lMax(1);
        Nu(:,2) = randn(n0, 1)./lMax(2);
        Nu(:,3) = randn(n0, 1)./lMax(3);
    otherwise
        fii     = 2*pi*Gamma(:,2);
        th      = acos(1 - 2*Gamma(:,3));

        r       = getInverse(params);
        Nu(:,1) = r'.*cos(fii).*sin(th)./lMax(1);
        Nu(:,2) = r'.*sin(fii).*sin(th)./lMax(2); 
        Nu(:,3) = r'.*cos(th)./lMax(3);
end

end
