function F = hypergeomLaplace(a,b,c,x)

% t = x.*(b - a) - c;
% y = (2*a)./(sqrt(t.^2 - 4*a*x*(c - b)) - t);
% j21 = a*(1 - y).^2 + (c - a)*y.^2 - b*(x.^2).*(y.^2).*((1 - y).^2)./(1 - x.*y).^2;
% h = ((2*pi)^(1/2)).*(beta(a, c - a)^-1).*(j21.^(-1/2)).*(y.^a).*((1 - y).^(c - a)).*((1 - x.*y).^(-b));

if x<-10^17
    F=-1./x;
else
    F=hypergeomLaplaceRAW(a,b,c,x)./hypergeomLaplaceRAW(a,b,c,0);
end

end

function F1=hypergeomLaplaceRAW(a,b,c,x)
t = x.*(b - a) - c;
y = (2*a)./(sqrt(t.^2 - 4*a*x*(c - b)) - t);
j21 = a*(1 - y).^2 + (c - a)*y.^2 - b*(x.^2).*(y.^2).*((1 - y).^2)./(1 - x.*y).^2;
F1 = ((2*pi)^(1/2)).*(beta(a, c - a)^-1).*(j21.^(-1/2)).*(y.^a).*((1 - y).^(c - a)).*((1 - x.*y).^(-b));
end 