function [fnnPerc, embTimes] = mdFnn(data, tau, varargin)
%MDFNN The mdfnn function computes the percentage of false nearest
%   neighbors for multidimensional time series as a function of embedding
%   dimension based on Kennel, M. B., Brown, R., & Abarbanel, H. D. (1992).
%   Determining embedding dimension for phase-space reconstruction using a
%   geometrical construction. Physical review A, 45, 3403.
%
%   Inputs:
%   Required arguments:
%    data = an m*n matrix, were m is the number of data points and n is the
%    number of dimensions of the time series.
%
%    tau = time delay for embedding.
%
%   Optional arguments:
%    maxEmb = maximum number of embedding dimensions that are considered.
%     Default is 10.
%
%    doPlot = plot false nearest neighbors if 1, supress plot if 0.
%     Default is 1.
%
%    noSamples = number of randomly drawn coordinates from phase-space used
%     to estimate the percentage of false-nearest neighors.
%     eflaut is 500.
%
%    Rtol = First distance criterion for separating false neighbors.
%     Default is 10.
%
%    Atol = Second distance criterion for separating false neighbors.
%     Default is 2.
%
%   Outputs:
%    fnnPerc = percentage of false neighbors for each embedding.
%  
%    embTimes = number of times the multidimensional time series was
%     embedded with delay tau. Note, that embTimes = 1 means no embedding,
%     and embTimes = 2 means embedding the multidimensional time series
%     once. Hence, the resulting falase neibor percentages result from
%     embeddeding data (embTimes - 1) times
%
%   Version: 1.0, 22 June 2018
%   Authors:
%     Sebastian Wallot, Max Planck Insitute for Empirical Aesthetics
%     Dan Moenster, Aarhus University
%
%   Reference:
%     Wallot, S., \& M{\o}nster, D. (2018). Calculation of average mutual
%     information (AMI) and false-nearest neighbors (FNN) for the
%     estimation of embedding parameters of multidimensional time-series in
%     Matlab. Front. Psychol. - Quantitative Psychology and Measurement
%     (under review)

%
% Parse and validate the input
%
parser = inputParser;

% Required parameter tau, must be positive integer
checkTau = @(x) validateattributes(x, {'numeric'}, {'positive', 'numel', 1});

% Optional parameter: maxEmb
defaultMaxEmb = 10;
checkMaxEmb = @(x) isnumeric(x) && isscalar(x) && (x > 1) && (mod(x, 1) == 0);

% Optional parameter: doPlot, must be either 0 or 1 or true or false
defaultDoPlot = 1;
checkDoPlot = @(x) ((x == 0) || (x == 1)) || ((x ==  true) || (x == false));

% Optional parameter: numSamples
defaultNumSamples = 500;
checkNumSamples = @(x) isnumeric(x) && isscalar(x) && (x > 1) && (mod(x, 1) == 0);

% Optional parameter: Rtol
checkRtol = @(x) validateattributes(x, {'numeric'}, {'positive', 'numel', 1});
defaultRtol = 10;

% Optional parameter: Atol
checkAtol = @(x) validateattributes(x, {'numeric'}, {'positive', 'numel', 1});
defaultAtol = 2;

% Do the parsing
addRequired(parser, 'data', @checkData);
addRequired(parser, 'tau', checkTau);
addOptional(parser, 'maxEmb', defaultMaxEmb, checkMaxEmb);
addOptional(parser, 'doPlot', defaultDoPlot, checkDoPlot);
addOptional(parser, 'numSamples', defaultNumSamples, checkNumSamples);
addOptional(parser, 'Rtol', defaultRtol, checkRtol);
addOptional(parser, 'Atol', defaultAtol, checkAtol);
parse(parser, data, tau, varargin{:});

% Get the optional arguments if provided. Otherwise the specified defaults
% are used.
maxEmb = parser.Results.maxEmb;
doPlot = parser.Results.doPlot;
numSamples = parser.Results.numSamples;
Rtol = parser.Results.Rtol;
Atol = parser.Results.Atol;

%
% Function to check input data, used by the input parser
%
function check = checkData(x)
   check = false;
   if (~isnumeric(x))
       error('Input is not numeric');
   elseif (numel(x) <= 1)
       error('Input must be a vector or matrix');
   else
       check = true;
   end
end

%
% Now do the actual work to find FNN
%
dims = size(data,2); % get dimensionality of time series
fnnPerc = 100; % first FNN
Ra = sum(var(data)); % estimate of attractor size
if (length(data)-tau*(maxEmb-1)) < numSamples % check whether enough data points exist for random sampling
    numSamples = length(data)-tau*(maxEmb-1);
    samps = 1:1:numSamples;
else
    samps = 1:1:(length(data)-tau*(maxEmb-1)); % generate random sample
    samps = sortrows([rand(length(samps),1) samps']);
    samps = samps(1:numSamples,2);
end
embData = [];
for i = 1:maxEmb % embed data
    embData = [embData data(1+(i-1)*tau:end-(maxEmb-i)*tau, 1:dims)];
    dists = pdist2(embData.^2,embData.^2);
    r2d1 = [];
    yRd1 = [];
    for j = 1:numSamples % get nearest neighbors and distances
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
            % check whether neighbors are false
            if sqrt((temp(yRd1(j)) - r2d(j))/r2d(j)) > Rtol || abs(temp(yRd1(j)) - r2d(j))/Ra > Atol
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
if (doPlot)  % plot results
    figure()
    plot(fnnPerc,'--o')
    title('False nearest neighbour function')
    xlabel('Number of embeddings')
    ylabel('Percentage of false nearest neighbors')
    text(floor(maxEmb/2),90,['Parameters: maxEmb = ' num2str(maxEmb) ', tau = ' num2str(tau) '.'])
    s=size(data);
    text(floor(maxEmb/2),80,['Data: # datapoints =' num2str(length(data)) ', # dimensions = ' num2str(s(2))  '.'])
end
end