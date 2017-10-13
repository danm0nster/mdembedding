function delay = mdDelay(data, varargin)
%MDDELAY Estimates time delay for embedding of multivariate times series.
%   The function plots the mutual information for multivariate times series
%   data, so the user can estimate the optimal value of the time delay for
%   embedding the data. The function also returns an estimate of the
%   optimal time delay, using simple methods, such as the mean of the lag
%   for which the auto mutual information for each of the variables
%   (columns) is less than 1/e. A single value. This is the uniform
%   multivariate embedding method.
%
%   Other methods may be added in a later version.
%
%   1 Brief description 2 Full-syntax call 3 Input &
%   output explanations
%     3.1 Alternative 1: List of inputs and outputs
%   Required arguments:
%     data - a matrix with multiple timeseries, one in each column.
%
%     3.2 Alternative 2: Expanded-syntax style
%   4 Notes
%   5 Examples
%   6 Methods and properties
%   7 See also statement
%   8 Repeat of full-syntax call
%   9 Showdemo
%
%   Authors: Sebastian Wallot and Dan M{\o}nster

%
% Parse and validate the input
%
par = inputParser;

% Optional p
defaultPlot = 'mean';
validPlots = {'mean', 'all', 'none'};
checkPlot = @(x) any(validatestring(x,validPlots));

addOptional(par,'plot',defaultPlot,checkPlot)
addRequired(par,'data', @checkdata);
   parse(par,data);

% TODO: These constants should be made parameters of the function.
nbins = 10;
maxlag = 10;
threshold = 1/exp(1);
[~, ncol] = size(data);

%
% Calculation of the mutual information as a function of time lag
%

% Allocate a matrix, where each column will be the auto mutual information
% as a function of time lag for a variable in the input data.
auto_mi = zeros(maxlag, ncol);

% Allocate a vector to hold the estimated optimal time lag for each
% dimension.
lags = zeros(1,ncol);

for c=1:ncol
    info = mi(data(:,c), nbins, maxlag, 'silent');
    % For each lag there is a 2x2 matrix with identical elements, since it
    % is the AUTO mutual information. We onlt need one of these, e.g.,
    % (1,1), and the result is squeezed to get rid of the extra dimensions.
    auto_mi(:,c) = squeeze(info(1,1,:));
    lags(c) = findFirstBelowThreshold(auto_mi, threshold);
end
plotMeanMI(auto_mi);
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

function lag = findFirstBelowThreshold(ami, threshold)
    % First find the first element below the threshold. Then test whether
    % an element below the threshold was found, and recover if this is not
    % the case.
    idx = find(ami < threshold, 1, 'first');
    if isempty(idx)
        disp("No value below threshold found. Will use minium instead");
        % If there is more than one elemtent that has the minimum value
        % the min() function returns the first one.
        [~, idx] = min(ami);
    end
        % A value of the index idx = 1 corresponds to lag = 0, so 1 is
        % subtracted from the index to get the lag.
        lag = idx - 1;       
end

function plotMeanMI(mi)
    [nlag, ncol] = size(mi);
    for c = 1:ncol
        plot(0:nlag, mi(:,c), 'b');
    end
end

% TODO: Write local function to calculate multual information
