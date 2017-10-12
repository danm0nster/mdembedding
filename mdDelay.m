function delay = mdDelay(data)
%MDDELAY Estimates time delay for embedding of multivariate times series.
%   The function calculates the value of the time delay for which the auto
%   mutual information for each of the variables (columns) is less than
%   1/e, and returns the mean value of these as the estimate of the optimal
%   time delay. This is the uniform multivariate embedding method.
%
%   Other methods will be added in a later version.
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
addRequired(par,'data', @checkdata);
   parse(par,data);

% TODO: These constants should be made parameters of the function.
nbins = 10;
maxlag = 10;
threshold = 1/exp(1);
[~, ncol] = size(data);
% Allocate a vector to hold the optimal time lag for each dimension.
lags = zeros(ncol);
for c=1:ncol
    info = mi(data(:,c), nbins, maxlag, 'silent');
    % For each lag there is a 2x2 matrix with identical elements, since it
    % is the AUTO mutual information. We onlt need one of these, e.g.,
    % (1,1), and the result is squeezed to get rid of the extra dimensions.
    auto_mi = squeeze(info(1,1,:));
    lags(c) = find(auto_mi < threshold, 1, 'first');
end
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
