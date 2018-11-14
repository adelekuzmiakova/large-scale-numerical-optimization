function [alpha, r] = optimize_for_alpha(p, lower, upper, x)
pos_p_indices = find(p > 0);
neg_p_indices = find(p < 0);

upper_bound_ratio_vector = (upper(pos_p_indices) - x(pos_p_indices))./p(pos_p_indices)
lower_bound_ratio_vector = (lower(neg_p_indices) - x(neg_p_indices))./p(neg_p_indices)

[min_upper_bound_ratio, index_up] = min(upper_bound_ratio_vector)
[min_lower_bound_ratio, index_low] = min(lower_bound_ratio_vector)

r_u = pos_p_indices(index_up)
r_l = neg_p_indices(index_low)

if min_upper_bound_ratio <= min_lower_bound_ratio
    alpha = min_upper_bound_ratio
    r = r_u
else
    alpha =  min_lower_bound_ratio
    r = r_l 
    
end

if isinf(alpha)
    r = 0
end
end