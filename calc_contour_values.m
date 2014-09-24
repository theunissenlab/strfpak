function [fracvalues] = calc_contour_values(stim_stat, fracpower);
% Calculates the contour values that capture a percent of the total power
% given in the array fracpower.  It returns fracvalues which has the same
% dimension as fracpower
debug_plot = 0;

% Reshape modulation spectrum to make it into a vector. 
X = reshape(stim_stat, 1, prod(size(stim_stat)));

% Sort into ascending order
XS = fliplr(sort(X));

% Find cumulative sum in normalized units
cumxs = cumsum(XS)./sum(XS);

nf = length(fracpower);
fracvalues = zeros(1, nf);
ind = cell(1, nf);

for i=1:nf
    diffxs = abs(cumxs-fracpower(i));
    mindiff = min(diffxs);    
    ind{i} = find(diffxs==mindiff);
    fracvalues(i) = XS(ind{i}(1));
end

% Debugging stuff
if (debug_plot)
    current_fig = gcf;
    figure;
    plot(cumxs);
    hold;
    for i=1:nf
        plot ([ind{i}(1) ind{i}(1)], [0 cumxs(ind{i}(1))], 'k');
        plot ([0 ind{i}(1)], [cumxs(ind{i}(1)) cumxs(ind{i}(1))], 'k');
    end
    figure(current_fig);
end
    
