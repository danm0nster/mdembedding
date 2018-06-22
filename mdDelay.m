function delay = mdDelay(data, varargin)
%MDDELAY Estimates time delay for embedding of multivariate times series.
%   The function plots the mutual information for multivariate times series
%   data, so the user can estimate the optimal value of the time delay for
%   embedding the data. The function also returns an estimate of the
%   optimal time delay, using simple methods, such as the mean of the lag
%   for which the auto mutual information for each of the variables
%   (columns) is less than a threshold, such as 1/e.
%
%   This function currently implements the uniform multivariate embedding
%   method.
%
%   Inputs:
%   Required arguments:
%     data - a matrix with multiple timeseries, one in each column.
%
%   Optional arguments:
%     maxLag: The maximum time lag for which AMI is computed. Default = 10.
%     numBins: The number of bins used to construct the histograms for
%       computing AMI. Default = 10.
%     threshold: The threshold value to select the delay when AMI drops
%       below threshold. Default = exp(-1)
%     plottype: Determines how the AMI is plotted. Possible values are
%     'mean', 'all', 'both', 'none'. Default = 'mean'
%
%   Outputs:
%     delay: The estimated optimal delay given the input and critera.
%
%   Version: 1.0, 22 June 2018
%   Authors:
%     Sebastian Wallot, Max Planck Insitute for empirical aesthetics
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

% Optional parameter: plottype
defaultPlotType = 'mean';
validPlotTypes = {'mean', 'all', 'both', 'none'};
checkPlotType = @(x) any(validatestring(x, validPlotTypes));

% Optional parameter: numBins
defaultNumBins = 10;
checkNumBins = @(x) validateattributes(x, {'numeric'}, {'positive', 'numel', 1});

% Optional parameter: maxLag
defaultMaxLag = 10;
checkMaxLag = @(x) validateattributes(x, {'numeric'}, {'positive', 'numel', 1});

% Optional parameter: threshold
defaultThreshold = exp(-1);
% checkThreshold = @(x) validateattributes(x, {'numeric'}, {'positive'});

addRequired(parser, 'data', @checkdata);
addOptional(parser, 'plottype', defaultPlotType, checkPlotType);
addOptional(parser, 'numBins', defaultNumBins, checkNumBins);
addOptional(parser, 'maxLag', defaultMaxLag, checkMaxLag);
addOptional(parser, 'threshold', defaultThreshold, @checkThreshold);
parse(parser, data, varargin{:});

% Get the optional arguments if provided. Otherwise the specified defaults
% are used.
numBins = parser.Results.numBins;
maxLag = parser.Results.maxLag;
threshold = parser.Results.threshold;

[~, ncol] = size(data);

%
% Calculation of the mutual information as a function of time lag
%

% Allocate a matrix, where each column will be the auto mutual information
% as a function of time lag [0; maxlag] for a variable in the input data.
auto_mi = zeros(maxLag + 1, ncol);

% Allocate a vector to hold the estimated optimal time lag for each
% dimension.
lags = zeros(1, ncol);

for c=1:ncol
    auto_mi(:,c) = autoMI(data(:, c), numBins, maxLag);
    lags(c) = findFirstBelowThreshold(auto_mi(:,c), threshold);
end

%
% Call the relevant plotting function
%
disp(['Plot type: ', parser.Results.plottype])
plotType = char(parser.Results.plottype);
switch plotType
    case 'mean'
        plotMeanMI(auto_mi, threshold);
    case 'all'
        plotAllMI(auto_mi, threshold);
    case 'both'
        plotMeanMI(auto_mi, threshold);
        plotAllMI(auto_mi, threshold);
    case 'none'
end
%
% Return the estimated optimal time lag
%
delay = mean(lags);
end

function check = checkdata(x)
   check = false;
   if (~isnumeric(x))
       error('Input is not numeric');
   elseif (numel(x) <= 1)
       error('Input must be a vector or matrix');
   else
       check = true;
   end
end

function check = checkThreshold(x)
   check = false;
   if (~isnumeric(x))
       error('The parameter threshold must be numeric');
   elseif (numel(x) ~= 1)
       error('The parameter threshold must be a scalar');
   elseif (x < 0)
       error('threshold must be a positive number');
   else
       check = true;
   end
end

function lag = findFirstBelowThreshold(ami, threshold)
    % First find the first element below the threshold. Then test whether
    % an element below the threshold was found, and recover if this is not
    % the case.
    idx = find(ami < threshold, 1, 'first');
    if isempty(idx)
        disp('No value below threshold found. Will use minimum instead');
        % If there is more than one elemtent that has the minimum value
        % the min() function returns the first one.
        [~, idx] = min(ami);
    end
        % A value of the index idx = 1 corresponds to lag = 0, so 1 is
        % subtracted from the index to get the lag.
        lag = idx - 1;       
end

% TODO: Add function findFirstLocalMinimum(ami)

function plotMeanMI(ami, threshold)
    [nlag, ncol] = size(ami);
    maxlag = nlag - 1;
    % Compute a vector with the mean of each row.
    y = mean(ami, 2);
    % Vector with standard deviation of each row.
    stddev = std(ami, 0, 2);  
    % Construct a vector with lags for x-axis.
    x = (0:maxlag)';
    figure();
    hold off
    % Plot shaded area indicating standard deviation if it is non-zero.
    if ~max(stddev) == 0
        yu = y + stddev;
        yl = y - stddev;
        % To make this work the vectors need to be transposed to
        % become row vectors.
        fill([x' fliplr(x')], [yu' fliplr(yl')], [.9 .9 .9], 'linestyle', 'none')
        hold on
    end
    plot(x, y, 'b')
    refline(0, threshold)
    limits = ylim;
    ylim([0, limits(2)]);
    xlabel('Time lag')
    ylabel('Mutual Information')
end

function plotAllMI(ami, threshold)
    [nlag, ncol] = size(ami);
    maxlag = nlag - 1;
    figure();
    hold off
    for c = 1:ncol
        plot(0:maxlag, ami(:, c), 'b');
        hold on
    end
    refline(0, threshold)
    limits = ylim;
    ylim([0, limits(2)]);
    xlabel('Time lag')
    ylabel('Mutual Information')
end
