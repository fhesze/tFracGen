function y = getInverse(params)

%% This function draws radom deviates according to a PDF by inverting the
% associated CDF.
%
% Copyright (C) 2013  Falk Heße, Vladyslav Prykhodko.
% The terms of the license agreement can be found in the README.txt.
%

if params.dim == 1
    for i = 1:params.nMode 
        y(i) = getInverse1D(params);       
    end      
end

if params.dim == 2
    for i = 1:params.nMode 
        y(i) = getInverse2D(params);
    end
end

if params.dim == 3
    for i = 1:params.nMode 
        y(i) = getInverse3D(params);
    end
end

end


%% One-dimensional fields -------------------------------------------------
%
%

function y = getInverse1D(params)
    
    x0 = 2.*rand - 1;
    y = fzero(@my_fun, [tan(x0.*pi./2)]);   

    function f = my_fun(x) 
        f = CumDistFunc(x, params) - x0;  
    end

end

%% Two-dimensional fields -------------------------------------------------
%
%

function y = getInverse2D(params)
    
    x0 = 2.*rand - 1;
    y = fzero(@my_fun, [tan(x0.*pi./2)]);   

    function f = my_fun(x) 
        f = CumDistFunc(x, params) - x0;     
    end

end

%% Three-dimensional fields -----------------------------------------------
%
%

function y = getInverse3D(params)
    
    x0 = rand;

    switch params.func   
        case 'Gauss'
            y = fzero(@my_fun, 1);         
        case 'Exp'
            y = fzero(@my_fun, 1);           
        case 'tFracGauss'
            y = fzero(@my_fun, [10^-10 10^10]);          
        case 'tFracExp'
            y = fzero(@my_fun, 1);
    end

    function f = my_fun(x) 
        f = CumDistFunc(x, params) - x0;
    end

end

