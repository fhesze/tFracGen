function u = FourMeth(x, params)

%% This function computes a Gaussian random field according to the
%% Fourier method.
%
% Copyright (C) 2013  Falk Heße, Vladyslav Prykhodko.
% The terms of the license agreement can be found in the README.txt.
%

if params.dim == 1
    [u] = Fourier1D(x, params);
end

if params.dim == 2
    [u] = Fourier2D(x, params);
end

% if params.dim == 3
%     [u] = Fourier3D(x(1,:), x(2,:), x(3,:), params);
% end


end

%% One-dimensional fields -------------------------------------------------
%
%

function u = Fourier1D(x ,params)

lMax = params.X.max;
lMin = params.X.min;
corrLen = params.X.max;
H       = params.H;
flag    = params.func;
Nk      = params.nMode;

switch flag
    case 'Gauss'
        Lk = 10.*lMax;
    case 'Exp'
        Lk = 10.*lMax;
    case 'Kolm'
        Lk = 100.*lMax;
    case 'tFracGauss'
        Lk = 100.*lMax;
    case 'tFracExp'
        Lk = 100.*lMax;
end

dk = Lk/(Nk-1);
i = 1:1:Nk;
k(1:1:Nk) = (i-1)*dk;

% gaussion-distributed random variable for sampling
ksi1 = randn(Nk, 1);
ksi2 = randn(Nk, 1);
    
u = zeros(1, length(x));

theta = k'*x;

for j = 1:Nk

    A = ksi1(j)*cos(theta(j,:));
    B = ksi2(j)*sin(theta(j,:));
    
    switch flag
        case 'Gauss'  
            S = SpecDens(k(j), params);
        case 'Exp'
            S = SpecDens(k(j), params);
        case 'Kolm'
            S = SpecDens(k(j), params);
        case 'tFracGauss'    
%             yl=(C*pi^(H-0.5)./(2*(k(j).^2)^(1/2+H)))*(gamma(0.5+H)-gamma(0.5+H)*(1-gammainc(k(j).^2/(pi*nl^2),0.5+H))) ;
%             yu=(C*pi^(H-0.5)./(2*(k(j).^2)^(1/2+H)))*(gamma(0.5+H)-gamma(0.5+H)*(1-gammainc(k(j).^2/(pi*nu^2),0.5+H))) ;
%             S = yl - yu;
            S = SpecDens(k(j), params);
        case 'tFracExp'
%             yl=(C/(pi*(1+2*H)*nl^(1+2*H))).*hypergeomEuler(1,0.5+H,1.5+H,-k(j).^2/nl^2);
%             yu=(C/(pi*(1+2*H)*nu^(1+2*H))).*hypergeomEuler(1,0.5+H,1.5+H,-k(j).^2/nu^2);
%             S = yl - yu;
            %y= 1/2*C*((-k.^2).^(-(1/2) - H)).*(pi^(-(1/2) - 1/2)).*mybetainc(-(k.^2/nl^2), 1/2 + H, 0);
            S = SpecDens(k(j), params);
    end
    u = u + sqrt(2.*dk).*sqrt(S)*(A + B);
end


end

%% Two-dimensional fields -------------------------------------------------
%
%

function u = Fourier2D(x, params)

nl      = 1/params.X.max;
n       = nl;
H       = params.H;
flag    = params.func;
n0      = params.nMode;

[XI,YI] = meshgrid(x(1,:), x(2,:));

%Diskretisierung von Spektralraum
Nx      = n0;
dkx     = 0.1;
kx      = 0:dkx:(dkx*(Nx-1));

Ny      = n0;
dky     = dkx;
ky      = 0:dky:(dky*(Ny-1));

u       = zeros(length(x(1,:)),length(x(2,:)));
Ksi1    = randn(n0, n0);
Ksi2    = randn(n0, n0);

for i = 1:Nx  
    for j = 1:Ny
        
        switch flag              
            case 'Gauss'
                S = 4*exp(-(kx(i)^2 + ky(j)^2)/(n*pi^2))/(pi^3*n);  
            case 'Exp'
                S = 4*gamma(3/2)/(n^2*pi^(3/2)*(1 + (kx(i)^2 + ky(j)^2)/n^2)^(3/2));
            % case 'Lor'
            %       S = SpecDens(kx(i),ky(j), [0.5 corrLen H], 'Lor', 1); 
            case 'Kolm'
                S = (kx(i)^2+ky(j)^2)^(-5/3);
            case 'tFracGauss'
                S = (H*(kx(i)^2 + ky(j)^2).^(-1 - H)).*((1./nl^2).^(-H)).*(pi^(-1 + H)).*(gamma(1+H)-gamma(1+H)*(1-gammainc((kx(i)^2 + ky(j)^2)/(pi*nl^2),1+H)));
                if (kx(i)==0)&&(ky(j)==0)
                     S = H/((1 + H)*nl^2*pi^2);
                end
            case 'tFracExp'
                S = (H*hypergeomEuler(3/2, 1 + H, 2 + H, -((kx(i)^2+ky(j)^2)/nl^2)) )./((2 + 2*H)* nl^2*pi);
        end
        
        theta = 2.*pi.*(kx(i).*XI + ky(j).*YI);
        si = sqrt(2*dkx*dky*S);
        u = u + si.*(Ksi1(i,j).*cos(theta) + Ksi2(i,j).*sin(theta));
        
    end
end


end

