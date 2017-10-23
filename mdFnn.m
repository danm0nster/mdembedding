function [fnnPerc, embTimes] = mdFnn(data, tau, maxEmb, doPlot, noSamples, Rtol, Atol)
% The mdfnn function computes the percentage of false nearest neigbors
% for multidimensional time series as a function of embedding dimension
% based on Kennel, M. B., Brown, R., & Abarbanel, H. D. (1992). Determining
% embedding dimension for phase-space reconstruction using a geometrical
% construction. Physical review A, 45, 3403.
%
% Inputs:
%  data = an m*n matrix, were m is the number of data points and n is the
%  number of dimensions of the time series.
%
%  tau = time delay for embedding.
%
%  maxEmb = maximum number of embedding dimensions that are considered.
%  Default is 10.
%
%  doPlot = plot false nearest neighbors if 1, supress plot if 0.
%  Default is 1.
%
%  noSamples = number of randomly drawn coordinates from phase-space used
%  to estimate the percentage of false-nearest neighors.
%  Deflaut is 500.
%
%  Rtol = First distance criterion for separating false neighbors.
%  Default is 10.
%
%  Atol = Second distance criterion for separating false neighbors.
%  Default is 2.
%
% Outputs:
%  fnnPerc = percentage of false neighbors for each embedding.
%  
%  embTimes = number of times the multidimensional time series was embedded
%  with delay tau. Note, that embTimes = 1 means no embedding, and
%  embTimes = 2 means embedding the multidimensional time series once.
%  Hence, the resulting falase neibor percentages result from embeddeding
%  data (embTimes - 1) times
%
% Version: 1.0,, 22 October 2017
% by Sebastian Wallot, Max Planck Insitute for empirical aesthetics
% & Dan Monster, Department of Economics, Aarhus University
%
% Reference:
%  Wallot, S., & Monster, D. (under review). Calculation of average mutual
%  information (AMI) and false-nearest neighbors (FNN) for the estimation
%  of embedding parameters of multidimensional time-series in Matlab. ???

% check inpus
if exist('data','var')
else
error('No input data specified.');
end

if exist('tau','var')
else
    error('No time delay tau specified.');
end

if exist('maxEmb','var')
else
    maxEmb = 10;
end

if exist('doPlot','var')
else
    doPlot = 1;
end

if exist('noSamples','var')
else
    noSamples = 500;
end

if exist('Rtol','var')
else
    Rtol = 10;
end

if exist('Atol','var')
else
    Atol = 2;
end

dims = size(data,2); % get dimensionality of time series
fnnPerc = 100; % first FNN
Ra = sum(var(data)); % estimate of attractor size
if (length(data)-tau*(maxEmb-1)) < noSamples % check whether enough data points exist for random sampling
    noSamples = length(data)-tau*(maxEmb-1);
    samps = 1:1:noSamples;
else
 samps = 1:1:(length(data)-tau*(maxEmb-1)); % generate random sample
 samps = sortrows([rand(length(samps),1) samps']);
 samps = samps(1:noSamples,2);
end
embData = [];
for i = 1:maxEmb % embedd data
    embData = [embData data(1+(i-1)*tau:end-(maxEmb-i)*tau, 1:dims)];
    dists = pdist2(embData.^2,embData.^2);
    r2d1 = [];
    yRd1 = [];
    for j = 1:noSamples % get nearest neighbors and distances
       [temp, coord] = sort(dists(:,samps(j)));
       r2d1(j) = temp(2);
       yRd1(j) = coord(2);
    end
    if i == 1
        r2d = r2d1;
        yRd = r2d1;
    else
        fnnTemp = [];
        for j = 1:length(r2d1)
            temp = dists(:,samps(j));
           if sqrt((temp(yRd1(j)) - r2d(j))/r2d(j)) > Rtol || abs(temp(yRd1(j)) - r2d(j))/Ra > Atol % check whether neighbors are false
               fnnTemp(j) = 1;
           else
               fnnTemp(j) = 0;
           end
        end
        fnnPerc(i) = 100*sum(fnnTemp)/length(fnnTemp); % compute percentage of FNN
        r2d = r2d1;
        yRd = r2d1;
    end
end
embTimes = 1:1:maxEmb;
if doPlot == 1 % plot results
plot(fnnPerc,'--o')
title('False nearest neighbour function')
xlabel('Number of embeddings')
ylabel('Percentage of false nearest neighbors')
text(floor(maxEmb/2),90,['Parameters: maxEmb = ' num2str(maxEmb) ', tau = ' num2str(tau) '.'])
s=size(data);
text(floor(maxEmb/2),80,['Data: # datapoints =' num2str(length(data)) ', # dimsensions = ' num2str(s(2))  '.'])
end
end