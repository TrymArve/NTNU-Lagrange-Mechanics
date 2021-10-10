function[t_conv, i_conv] = TimeOfConvergence(t,x,s,tol,relative)
% This function finds at what time and index a verctor-series x converges
% according to your definition("s","tol")

% t - timestamps for x values
% x - vector(columns) of values that are converging
%     (rows show the progression of the values)
% s - number of consecutive samplepoints required to be within tolerance(must be: s > 1, and integer)
% tol - maximum deviation from mean allowed(absolute or relative)
%       (tolerance can be one number that should be used for all values,
%        or a vector(column) containing tolerances for each element of x)
% absolute - absolute tolerance instead of relative (any input argument)

% t_conv - the time at which x has converged (according to specifications)
% i_conv - the index at which x has converged (according to specifications)


if nargin == 4
    Relative = 0;
elseif nargin == 5
    Relative = 1;
else
    error('Check number of inputs!')
end

t_conv = 'N';
i_conv = NaN;
L = length(x(1,:));
first = 1;
    for i = L-s+1:-1:1
        
        X = x(:,i:end);
        Mean = (sum(X')/(L-i+1))';
        Deviation = X - Mean;
        max_deviation =  max(abs(Deviation'))';
            if Relative
                max_deviation = max_deviation./Mean;
            end
        violates = (max_deviation > tol);
            if ~isempty(find(violates==1,1)) && ~first
                t_conv = t(i);
                i_conv = i;
                break;
            end
        first = 0;

    end

    if t_conv == 'N'
        t_conv = 0;
        disp(' ')
        text = '---WARNING: The vector does not converge within the given time frame';
        disp(text)
        disp(' ')
    end
    
end