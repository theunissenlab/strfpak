% function pal=redblue(expon);
%
% passed to colormap() function to set blue-->white-->red palette
%
% expon is the exponent of coloration. 
% (ie, 1= linear, <1 more white, >1 more color)
% 
function pal=redblue(expon);

if ~exist('expon','var'),
   expon=0.75;
end

skbottom=[linspace(0,1,32).^expon' linspace(0,1,32).^expon' ones(32,1)];
sktop=[ones(32,1) linspace(1,0,32).^expon' linspace(1,0,32).^expon'];
pal=[skbottom;sktop];
