% function w=hanning2(n)
%
% create a circularly symmetric 2d hanning window
% ripped off matlab hanning 1d code
%
function w=hanning2(varargin)

% Select the sampling option
if nargin == 1,
   sflag = 'symmetric';
else
   sflag = lower(varargin{2});
end
n=varargin{1};

% Allow partial strings for sampling options
allsflags = {'symmetric','periodic'};
sflagindex = strmatch(sflag, allsflags);
if length(sflagindex)~=1         % catch 0 or 2 matches
   error('Sampling flag must be either ''symmetric'' or ''periodic''.');
end
sflag = allsflags{sflagindex};
% Evaluate the window
switch sflag,
case 'periodic'
   % Includes the first zero sample
   w = [zeros(1,n); zeros(n-1,1) sym_hanning(n-1)];
case 'symmetric'
   % Does not include the first and last zero sample
   w = sym_hanning(n);
end

%---------------------------------------------------------------------
function w = sym_hanning(n)
%SYM_HANNING   Symmetric Hanning window. 
%   SYM_HANNING Returns an exactly symmetric N point window by evaluating
%   the first half and then flipping the same samples over the other half.

[xx,yy]=meshgrid((1:n)-(n+1)/2,(1:n)-(n+1)/2);

d=sqrt(xx.^2+yy.^2);

w = .5*(1 - cos(2.*pi*(d-(n+1)/2)/(n+1))); 
w(d>n/2)=0;







