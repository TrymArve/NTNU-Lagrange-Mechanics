function[t_conv, i_conv] = TimeOfConvergence(t,x,tol,N,relative)

%OPS! old version is now called TimeOfConvergenceOLD

% This function finds at what time and index a verctor-series x converges
% according to your definition("s","tol")

% t - timestamps for x values(row or column)
% x - vector(columns) of values that are converging
%     (rows show the progression of the values)
% N - Number of point that must be within tolerance to be considered
%      converged
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

t_conv = NaN;
i_conv = NaN;
[n,L] = size(x);
[m,~]= size(tol);
if n ~= m; error('NOTE: Tolerance vector must have the same length as the state vector');end
if L < N;
PMS = zeros(n,1);
first = 1;
    for i = L:-1:1
        
        X = x(:,i);
        PMS = PMS + X;
        Mean = PMS/(L-i+1);
        Deviation = X - Mean;
        if Relative
            Deviation = Deviation./Mean;
        end
        violates = (Deviation > tol);
        if ~isempty(find(violates==1,1))
            if first
                break;
            end
            t_conv = t(i);
            i_conv = i;
            break;
        end
        first = 0;

    end

    if isnan(i_conv)
%         disp(' ')
%         text = '---WARNING: The vector does not converge within the given time frame';
%         disp(text)
%         disp(' ')
          i_conv = '---WARNING: The vector does not converge within the given time frame';
    end
    
end