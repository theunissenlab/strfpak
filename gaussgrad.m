function [dy] = gaussgrad(beta,x)

%   if length(beta) != 3
%      error('beta parameter in gaussfit of wrong dimensions');
%   end


nl = length(x);
dy = zeros(3,nl);

for i=1:nl
    dx = x(i) - beta(2);
    expval = exp(-0.5*dx*dx/(beta(3)*beta(3)));
    dy(1,i) = expval;
    dy(2,i) = beta(1)*dx*expval/(beta(3)*beta(3));
    dy(3,i) = beta(1)*dx*dx*expval/(beta(3)*beta(3)*beta(3));
end

