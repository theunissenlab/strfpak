% function [H,lambda]=normalizetime(SR,tSA,nfactors,lambda)
% 
% created SVD 4/12/02 - hacked from normalize.m
% modified SVD 4/15/02 - scale lambdas according to eigenvalues
% 
% Inputs (basically the output from movxc.m or summed output from
%        multiple runs of it):
%      N - number of spatial dimensions
%      U - total number of time lags (eg, length([maxlag(1):maxlag(2)]))
%      M - number of response sets (or response dimensions)
% SR - raw spike-triggered average (correlation between stim and
%      resp) (N x U x M)
% sSA - spatial autocorrelation matrix (N x N x M or N x N and use
%       the same one for each response set--speeds things up)
% tSA - temporal autocorrelation matrix (U x U x M or leave empty
%       [], to skip temporal decorrelation)
% nfactors - number of regularization parameters to test (ranging
%            (10^3 through 10^-3)
% 
% Outputs:
% H - decorrelated kernel estimate(s), N x U x nfactors x M matrix
% lambda - regulatization parameters corresponding to each estimate
%          of H
% 
function [H,lambda]=normalizetime(SR,tSA,nfactors,lambda);

global BATQUEUEID

TMINSFS=0.0001;

disp('normalizetime.m: Doing single temporal normalization');

s=size(SR);
N=s(1);
U=s(2);
M=prod(s(3:end));
TSMOOTH=0;
if TSMOOTH>0,
   fprintf('TSMOOTH=%d\n',TSMOOTH');
end

if ~exist('lambda'),
   if ~exist('nfactors'),
      nfactors=10;
   end
   lambda=10.^(linspace(5,-1,nfactors));
else
   nfactors=length(lambda);
end
H=zeros([N,U,nfactors,s(3:end)]);

% as in XC, each response set is essentially independent here
for respidx=1:M,
   % decorr time if tSA included
   % normalize tSA to avoid double stimulus power over-compensation
   
   ttSA=zeros(U);
   for uu=1:U,
      ttSA(:,uu)=tSA((U+1-uu):(U*2-uu),respidx);
   end
   
   % compute SVD on ttSA and then pseudo-inverse
   [tU tS tV] = svd(ttSA);
   
   if 0
      teigs=diag(tS)./sum(diag(tS));
      tSAinv=tV*diag(1./diag(tS).*(teigs>TMINSFS))*tU';
      
      tSR=SR(:,:,respidx)*tSAinv';
      if TSMOOTH==1,
         tSR=conv2(tSR,[0.2 0.6 0.2],'same');
      end
      H(:,:,2,respidx)=tSR;
   end
   
   tD=diag(tS);
   tDnz=tD(find(tD>0));
   tDz=tD(find(tD<=0));
   
   n=length(tDnz);
   if nfactors>2,
      %lambda=10.^(log10(lambda./lambda(1))+log10(sum(tDnz.^2)));
      lambda=tDnz(round(n*linspace(0.1,0.9,nfactors))).^2;
   else
      lambda=[tDnz(round(n*0.4)).^2 tDnz(round(n*0.5)).^2  ];
   end
   
   
   for sfsidx=1:nfactors,
      % scale by regulatization paramater here.
      tDi=1./(tDnz+lambda(sfsidx)./tDnz);
      tSAinv=tV*diag([tDi;tDz])*tU';
      H(:,:,sfsidx,respidx)=SR(:,:,respidx)*tSAinv';
      if TSMOOTH==1,
         H(:,:,sfsidx,respidx)=conv2(H(:,:,sfsidx,respidx),...
                                     [0.2 0.6 0.2],'same');
      end
      %imagesc(tSAinv);
      %pause;
   end
end

