function [y] = gaussfit(beta,x)
% y = gaussfit(beta, x)
%    
% Created by FET.
%
    if length(beta) ~= 3
        error('beta parameter in gaussfit of wrong dimensions');
    end

    nl = length(x);
    y = zeros(1,nl);
    
    for i=1:nl
      if (beta(3) <= 0)
          y(1,i) = 0.0;
      else
          expval = ((x(i) - beta(2))*(x(i) - beta(2)))/(beta(3)*beta(3));
          y(1,i) = beta(1)*exp(-0.5*expval);
      end
    end
% ===========================================================
% End of gaussfit
% ===========================================================

