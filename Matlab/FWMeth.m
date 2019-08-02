function u=FWMeth(x, params)

%% This function computes a Gaussian random field according to the
%% Fourier-Wavelet method.
%
% Copyright (C) 2013  Falk Heße, Vladyslav Prykhodko.
% The terms of the license agreement can be found in the README.txt.

%Parameter of the method:
M=40;
b=10;
m=0:1:(M-1);
j=(-b+1):1:b;
%generate random numbers:
n=zeros(1,length(x)*length(m)*length(j));
for ix=1:1:length(x)
    for im=1:1:length(m)
        for ij=1:1:length(j)
           n((ix-1)*length(m)*length(j)+(im-1)*length(j) + ij ) = floor((2^(m(im)))*x(ix))+j(ij);
        end
     end
end
wn=unique(n);
gamma_mat=randn(M,length(wn));
%generate single realization of random field:
u=zeros(1,length(x));
for im=1:1:length(m)
        [ksi,fm]=fintegral(m(im),params);
        for ij=1:1:length(j)
           ste1 = (2^(m(im)))*x-floor((2^(m(im)))*x)-j(ij);
           ste=floor((2^(m(im)))*x)+j(ij);
           r=findvektor1(wn,ste);
           gm=gamma_mat(im,r);
           ff=interp1(ksi,real(fm),ste1);
           u=u+gm.*ff;
        end
end
end

function [ksi,fm]=fintegral(m,params)
dksi=0.1;
a=1/(2*dksi);
r=8;
N=2^r;
dk=2*a/N;
l=1:1:N;
k=-a+(l-0.5)*dk;
j=1:1:N;
g=zeros(1,length(k));
for i=1:1:length(k)
g(i)=gg(k(i),m,params);
end
ksi= -N*dksi*0.5+(j-1)*dksi;
G=dk*g.*exp(-2*pi*1i*( (N-1)*0.25 - (l-1)*0.5    ));
fm=zeros(1,N);
for j=1:1:N
fm(j)=exp(pi*1i*(j-1)*(1-1/N))*sum(G.*exp( -2*pi*1i*(j-1)*(l-1)/N  ));
end
end

function spek=gg(k,m1,params)
%Spektral density parameters:
lMax = params.X.max;
lMin = params.X.min;
H = params.H;
nl=1/lMax;
nu=1/lMin;
n=nl;
C=H*nl^(2*H);
flag=params.func;
    function y=spektrum(kk,flag)
        switch flag
    case 'Kolm'
         y=abs(kk).^(-5/3);
    case 'Exp'
         y=C/((n^(2*H+1))*pi*(1+kk.^2/n^2));
    case 'Gauss'  
         y=C/((n^(2*H+1))*pi)*exp(-kk.^2/(pi*n^2));
    case 'tFracExpo'
         yl=(C/(pi*(1+2*H)*nl^(1+2*H))).*hypergeomEuler(1,0.5+H,1.5+H,-kk.^2/nl^2);
         yu=(C/(pi*(1+2*H)*nu^(1+2*H))).*hypergeomEuler(1,0.5+H,1.5+H,-kk.^2/nu^2);
         y=yl-yu;
         %y= 1/2*C*((-k.^2).^(-(1/2) - H)).*(pi^(-(1/2) - 1/2)).*mybetainc(-(k.^2/nl^2), 1/2 + H, 0);
    case 'tFracGauss'    
         yl=(C*pi^(H-0.5)./(2*(kk.^2)^(1/2+H)))*(gamma(0.5+H)-gamma(0.5+H)*(1-gammainc(kk.^2/(pi*nl^2),0.5+H))) ;
         yu=(C*pi^(H-0.5)./(2*(kk.^2)^(1/2+H)))*(gamma(0.5+H)-gamma(0.5+H)*(1-gammainc(kk.^2/(pi*nu^2),0.5+H))) ;
         y=yl-yu;
        end
    end
spek = (2^(m1/2))*sqrt(spektrum((2^m1)*k,flag))*ftansformmeyer(k);
end

function y3 = ftansformmeyer(k) %%fi-dach
    y3 = -1i*sign(k)*exp(1i*pi*k)*bfunc(abs(k));
end

function y4 = bfunc(k)
    if ((k<=2/3) && (k>1/3)) 
        y4=sin(pi/2*vfunc(3*k-1)); 
    elseif ((k<=4/3) && (k>2/3)) 
        y4=cos(pi/2*vfunc(3*k/2-1));
    else
        y4=0;
    end  

end

function y5 = vfunc(xx)
    %p = 2;
    y5 = 2*( max(0,xx)^2 + max(0,xx-1)^2  -2*max(0,xx-1/2)^2 );
end

function b=findvektor1(X,a)
    b = zeros(1,length(a));
    for i=1:1:length(a)
        b(i) = find(X==a(i));
    end
end   


