function [alpha, r] = optimize_for_alpha(lower, upper, x, p)
% Inputs: search direction (p), lower bound (l), upper bound (u), x, such
% that x satisfies l <= x <= u

% Returns step size (alpha) and index r, such that alpha < infinity
% otherwise r = 0

% Divide p into positive and negative segments:

negative_ix = find(p < 0);
positive_ix = find(p > 0);

% Define the ratios for upper and lower bounds:
upper_bound_ratio = (upper(positive_ix) - x(positive_ix))./p(positive_ix);
lower_bound_ratio = (lower(negative_ix) - x(negative_ix))./p(negative_ix);

% Define the constraints (in form of minima) for upper and lower bound ratios:
[min_upper_bound_ratio, min_upper_ix] = min(upper_bound_ratio);
[min_lower_bound_ratio, min_lower_ix] = min(lower_bound_ratio);

upper_r = positive_ix(min_upper_ix)
lower_r = negative_ix(min_lower_ix)


% Find which of the two bounds (min. upper bound ration or min. lower bound 
% ratio are an actual constraint. Use the limiting ratio as the alpha to
% maintain feasibility:
if min_upper_bound_ratio > min_lower_bound_ratio
    alpha =  min_lower_bound_ratio;
    r = lower_r;
else
    alpha = min_upper_bound_ratio;
    r = upper_r;
end

% Declare that if alpha == infinity, r = 0
if isinf(alpha)
    r = 0
      
end
end