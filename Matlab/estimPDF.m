function [myPDF myVar] = estimPDF(data, containers)
      
%% This function estimates the PDF function of a Gaussian random ield.
%
% Copyright (C) 2013  Falk Heße, Vladyslav Prykhodko.
% The terms of the license agreement can be found in the README.txt.
%

% ----- preprocessing the data --------------------------------------------

[xLen yLen] = size(data);

% reshape mesh matrix --> mesh vector
vectZ = reshape(data, xLen*yLen, 1);

%-- random data sampling --------------------------------------------------

% number of stations
stations = 10000;

% random sample of measure stations
selectedStat = floor(rand(stations, 1)*length(vectZ)) + 1;

% coordinates of stations
% dataSel = data(selectedStat);
dataSel = data;

% sorting out outliers
dataSel = sort(dataSel);
dataSel((ceil(1.00.*length(dataSel))):end) = [];

% defining the edges of the containers
edges = linspace(0, round(containers./10));
% edges = [0 logspace(-2, 7, containers)];
    
% ----- calculate emperical PDF -------------------------------------------

% myPDF = hist(dataSel, containers);
% myVar = linspace(min(dataSel), max(dataSel), containers);
myPDF = histc(dataSel, edges);
myVar = edges;
myPDF = myPDF/trapz(myVar, myPDF);
	
end