function u = HybMeth(x, params)

%% This function computes a Gaussian random field according to the
%% Hybrid method.
%
% Copyright (C) 2013  Falk Heße, Vladyslav Prykhodko.
% The terms of the license agreement can be found in the README.txt.
%

if params.dim == 1
    [u] = Hybrid1D(x(1,:), params);
end

if params.dim == 2
    [u] = Hybrid2D(x(1,:), x(2,:), params);
end

% if params.dim == 3
%     [u] = Hybrid3D(x(1,:), x(2,:), x(3,:), params);
% end


end

%% One-dimensional fields -------------------------------------------------
%
%

function u = Hybrid1D(x ,params)

%% Defining Parameters ----------------------------------------------------

lMax 	= params.X.max;
lMin 	= params.X.min;
H       = params.H;
n       = params.nMode;                                                    % number of grid points in k

nP      = 10;                                                              % number of partitions
n0      = ceil(n./nP);                                                     % number of modes per partition

% predefined for efficiency
u       = zeros(1, length(x));                                             % solution

edges 	= logspace(1, 5, nP);
a       = [0 edges];
e       = [edges 100000000000];
    

%% Central Solver ---------------------------------------------------------

for k = 1:nP
    
    nu      = getHybridSet1D(params, [a(k) e(k)]);
    ksi     = randn(n0, 2);
    
    % gaussion-distributed random variable for the amplitudes
    ksi1 	= repmat(randn(n0, 1),1,length(x));
    ksi2 	= repmat(randn(n0, 1),1,length(x));

    si2(k)  = (myCumDistFunc(e(k), params) - myCumDistFunc(a(k), params));
    
    u 	= u + sqrt(si2(k)./n0).*sum(ksi1.*cos(nu'*x) + ksi2.*sin(nu'*x));
    
end

end

function Nu = getHybridSet1D(params, edge)
    
    n       = params.nMode;
    nP      = 10;                                                              % number of partitions
    n0      = ceil(n./nP); 

    gam     = rand(n0);

    Nu = getInverseB(params, edge, gam);

end

function y = getInverseB(params, edges, x0)

    for i = 1:length(x0)
        
        y(i) = getInvScal(params, edges, x0(i));
        
    end

end

function y = getInvScal(params, edges, x0)

    lMax = params.X.max;
    lMin = params.X.min;
    H = params.H;
    type = params.func;
    H2 = 2.*H;
    
    switch type
        case 'Gauss' 
            y = fzero(@Gauss, [0.5]);
        case 'Exp' 
            y = fzero(@Exp, [tan(x0.*pi./2)]);
        case 'tFracGauss'
            y = fzero(@TruncGauss, [edges(1) edges(2)]);
        case 'tFracExp'
            y = fzero(@TruncExp, [edges(1) edges(2)]);   
    end

    function f = Gauss(x)
        
        h = CumDistFunc(x, [lMax H], type);
        
%         h_a = myCumDistFunc(edges(1), [lMax lMin H], type);
%         h_e = myCumDistFunc(edges(2), [lMax lMin H], type);
%         
%         if x < edges(1)
%             h = 0;
%         elseif x > edges(2)
%             h = 1;
%         else
%             h = (myCumDistFunc(x, [lMax lMin H], type) - h_a)./(h_e - h_a);
%         end
      
        f = h - x0;
        
    end

    function f = Exp(x)

        h = CumDistFunc(x, [lMax H], type);
        
        f = h - x0;
        
    end

    function f = TruncGauss(x)
                
        h_a = myCumDistFunc(edges(1), params);
        h_e = myCumDistFunc(edges(2), params);
        
        if x < edges(1)
            h = 0;
        elseif x > edges(2)
            h = 1;
        else
            h = (myCumDistFunc(x, params) - h_a)./(h_e - h_a);
        end
 
        f = h - x0;
        
    end

    function f = TruncExp(x)
        
        h_a = myCumDistFunc(edges(1), params);
        h_e = myCumDistFunc(edges(2), params);
        
        if x < edges(1)
            h = 0;
        elseif x > edges(2)
            h = 1;
        else
            h = (myCumDistFunc(x, params) - h_a)./(h_e - h_a);
        end
        
        f = h - x0;
        
    end

end

function CDF = myCumDistFunc(x, params)

    lMax = params.X.max;
    lMin = 0;
    H   = params.H;
    flag = params.func;
    H2  = 2.*H;

    CDFMax = CumDistFunc(x, [lMax H], flag, 1);
    if lMin == 0
        CDFMin = 0;
    else
        CDFMin = CumDistFunc(x, [lMin H], flag, 1);
    end

    CDF = (lMax.^H2.*CDFMax - lMin.^H2.*CDFMin)./(lMax.^H2 - lMin.^H2);

end


%% Two-dimensional fields -------------------------------------------------
%
%

function u = Hybrid2D(x ,params)

%% Defining Parameters ----------------------------------------------------

lMax 	= params.X.max;
lMin 	= params.X.min;
H       = params.H;
n       = params.nMode;                                                    % number of grid points in k

nP      = 10;                                                              % number of partitions
n0      = ceil(n./nP);                                                     % number of modes per partition

% predefined for efficiency
u       = zeros(1, length(x));                                             % solution

edges 	= logspace(1, 5, nP);
a       = [0 edges];
e       = [edges 100000000000];
    

%% Central Solver ---------------------------------------------------------

for k = 1:nP
    
    uSum    = zeros(1, length(x));
    theta 	= zeros(n0, length(x));
    
    nu      = getHybridSet1D(params, [a(k) e(k)]);
    ksi     = randn(n0, 2);
    
    si2(k)  = (myCumDistFunc(e(k), params) - myCumDistFunc(a(k), params));
    
    for j = 1:n0
        
        theta   = nu(j)'*x;
        A       = ksi(j,1).*cos(theta);
        B       = ksi(j,2).*sin(theta);
        uSum    = uSum + sqrt(si2(k)./n0).*(A + B);
        
    end
    
    u       = u + uSum;
    
end

end
