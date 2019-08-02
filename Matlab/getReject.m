function y = getReject(params)

%% This function draws radom deviates according to a PDF by the Rejection
% method.
%
% Copyright (C) 2013  Falk Heße, Vladyslav Prykhodko.
% The terms of the license agreement can be found in the README.txt.
%

if params.dim == 1
    y = getReject1D(params);       
end

% if params.dim == 2
%     y = getReject2D(params);
% end
% 
% if params.dim == 3
%     y = getReject3D(params);
% end

end

%% One-dimensional fields -------------------------------------------------
%
%

function nu = getReject1D(params)

lMaxX   = params.X.max(1);
flag    = params.func;
num     = params.nMode;
H       = params.H;
nl      = 1./lMaxX;

k   = 5;
ngef= 0;
nu  = zeros(1,num);

while ngef < num
    u = rand;
    x0 = tan(pi*(rand - 1/2))./lMaxX;
    if (u*k*g(x0, params)) < SpecDens(x0, params)
        ngef = ngef + 1;
        nu(ngef) = x0;
    end
end

%%%-- auxillary functions -------------------------------------------------

% envelopping spetral-density function
function y = g(x, params)
    
    lMaxX   = params.X.max(1);
    nl      = 1./lMaxX;
    y       = nl./(pi.*(nl.^2 + x.^2));
    
end

%%%-- plotting functions --------------------------------------------------

% xx  = -100:0.01:100;
% plot(xx,f(xx),'g',xx,g(xx),'b',xx,k*g(xx),'r');
% hold on;
% plot(fzuf,f(fzuf),'.')
% hold off

end

%% Two-dimensional fields -------------------------------------------------
%
%

function nu = getReject2D(params)

lMaxX   = params.X.max(1);
flag    = params.func;
num     = params.nMode;
H       = params.H;
nl      = 1./lMaxX;

end


