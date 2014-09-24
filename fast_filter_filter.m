function s_forward = fast_filter_filter(forward, forwardJN_std, nstd);
% smooths out the filter for displaying or calculation purposes.
% Faster than filter_filter, but less fancy.
global zero_negatives
% Scales the filter everywhere by a sigmoid in forward/forwardJN_std, with
% inflection point at nstd, and a dynamic range from nstd - .5 to nstd + .5
if nstd > 0
    epsilon = 10^-8; % To prevent division by 0.

    factor = (1 + tanh(2.*(abs(forward)-abs(nstd*forwardJN_std))./(epsilon + abs(forwardJN_std))))/2;

    s_forward = factor .* forward;
else
    s_forward = forward;
end

if strcmp(zero_negatives,'yes')
    s_forward(find(s_forward < 0)) = 0;
end