function y=hypergeomEuler(a,b,c,z)

%% This function computes Gaussian hypergeometric function 2F1
%% z can be a number or array of dimensions 2 and 3
%% method based of Euler integral representation of 
%% a hypergeometric function using adaptive Gauss-Kronrod quadrature
%
% Copyright (C) 2013  Falk Heße, Vladyslav Prykhodko.
% The terms of the license agreement can be found in the README.txt.
%

y=zeros(size(z));    
for i=1:size(y,1)
    for j=1:size(y,2)
       for k=1:size(y,3)   
             F=@(t) (t.^(b-1)).*((1-t).^(c-b-1)).*((1-z(i,j,k).*t).^(-a));   
             y(i,j,k)=gamma(c)/(gamma(b)*gamma(c-b))*quadgk(F,0,1);
       end
    end
end
