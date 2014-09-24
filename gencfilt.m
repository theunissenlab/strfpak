% function [cfilt,cfiltconj,cmask,cmaskconj]=gencfilt(Xmax[=16],Ymax[=Xmax])
%
% select half of bins from a Xmax X Ymax matrix to remove redundancy
% from a 2d Fourier power representation
%
% example of bins in cfilt where points in a 16x16 fourier space get mapped 
%
% 0   0   0   0   0   0   0   0  29  38  48  59  71  84   98  113
% 0   0   8   9  11  14  18  23  30  39  49  60  72  85   99  114
% 0   0   0  10  12  15  19  24  31  40  50  61  73  86  100  115
% 0   0   0   0  13  16  20  25  32  41  51  62  74  87  101  116
% 0   0   0   0   0  17  21  26  33  42  52  63  75  88  102  117
% 0   0   0   0   0   0  22  27  34  43  53  64  76  89  103  118
% 0   0   0   0   0   0   0  28  35  44  54  65  77  90  104  119
% 0   0   0   0   0   0   0   0  36  45  55  66  78  91  105  120
% 1   0   0   0   0   0   0   0  37  46  56  67  79  92  106  121
% 2   0   0   0   0   0   0   0   0  47  57  68  80  93  107  122
% 3   0   0   0   0   0   0   0   0   0  58  69  81  94  108  123
% 4   0   0   0   0   0   0   0   0   0   0  70  82  95  109  124
% 5   0   0   0   0   0   0   0   0   0   0   0  83  96  110  125
% 6   0   0   0   0   0   0   0   0   0   0   0   0  97  111  126
% 7   0   0   0   0   0   0   0   0   0   0   0   0   0  112  127
% 0   0   0   0   0   0   0   0   0   0   0   0   0   0    0  128
%
% created SVD 2/2003
%
function [cfilt,cfiltconj,cmask,cmaskconj]=gencfilt(Xmax,Ymax)

if ~exist('Xmax','var'),
   Xmax=16
end
if ~exist('Ymax','var'),
   Ymax=Xmax;
end

pixcount=Xmax*Ymax;
xc=round((Xmax+1)/2);
yc=round((Ymax+1)/2);

if Xmax/2==round(Xmax/2),
   cmask=triu(ones(Xmax,Ymax))-diag([ones(Xmax./2,1);zeros(Xmax./2,1)]);
   cmask(1,1:(Xmax/2))=0;
   cmask((Ymax/2+1):Ymax-1,1)=1;
else
   cmask=triu(ones(Xmax,Ymax)) - diag([ones(xc-1,1);zeros(xc,1)]);
   cmask(1,Xmax)=0;  % throw away last thing to make it even
end

% find indices of non-zero points
cfilt=find(reshape(cmask,pixcount,1));

% flip around center points to get conjugate index matrix
[ii,jj]=ind2sub([Xmax Ymax],cfilt);
ii(find(ii>1))=2*xc-ii(find(ii>1));
jj(find(jj>1))=2*yc-jj(find(jj>1));

cfiltconj=sub2ind([Xmax Ymax],ii,jj);

cmaskconj=zeros(Xmax,Ymax);
cmaskconj(cfiltconj)=1;
