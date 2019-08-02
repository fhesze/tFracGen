function gamma_inc = myGammaInc(a, x)

%% Definition of the incomplete Gamma function according to Mathematica

gamma_inc = gamma(a).*gammainc(x, a, 'upper');

end