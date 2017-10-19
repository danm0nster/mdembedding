function ami = autoMI(x, nbins, maxlag)
    ami = zeros(maxlag+1,1);
    % Compute a histogram of x
    p_i = histcounts(x, nbins);
    % Normalize to make it a probability
    p_i = p_i / sum(p_i);
    % Save the original input vector
    xorig = x;
    % Loop over time lags from 0 to maxlag
    for tau=0:maxlag
        % Make a time delayed version of x
        xtau = circshift(xorig,-tau);
        % Remove the last tau elements of this vector and of x
        xtau = xtau(1:end-tau);
        x = xorig(1:end-tau);
        % Compute the two-dimensional histogram of x and xtau
        p_ij = histcounts2(x, xtau, [nbins nbins],'Normalization','probability');
        % Normalize the rows to get the joint probability the x is
        % in bin i and xtau is in bin j.
        % row_sum = sum(p_ij, 2);
        % If any rows sum to zero, replace the sum by 1
        % row_sum(row_sum == 0) = 1;
        % Divide each element with the row sum. 
        % p_ij = diag(1 ./ row_sum) * p_ij;
        % Compute the auto mutual information by summing over all
        % combinations of i and j
        for i=1:nbins
            for j=1:nbins
                if (p_i(i) ~= 0 && p_i(j) ~= 0 && p_ij(i,j) ~= 0)
                    ami(tau+1) = ami(tau+1) + p_ij(i,j) * log(p_ij(i,j) / (p_i(i) * p_i(j)));
                end
            end
        end
    end
end
