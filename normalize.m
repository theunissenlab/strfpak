% function [H,neigs]=normalize(SR,sSA,tSA,tolval)
%
% created SVD 3/20/02 - hacked from kerndodecorr.m
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
% tolval - vector of tolerance values (ie, % stimulus variance to
%          exclude from decorrelated STRF)
%
% Outputs:
% H - decorrelated kernel estimate(s), N x U x length(mineigs) x M matrix
% neigs - number of eigenvectors included in corresponding
%         decorrelated kernel estimate
%
function [H,neigs]=normalize(SR,sSA,tSA,tolval);

disp('normalize.m:');
global BATQUEUEID

TMINSFS=0.0001;

s=size(SR);
N=s(1);
U=s(2);
M=prod(s(3:end));

if size(sSA,1)~=N,
   disp('Spatial dimension must match in SR and sSA!');
   return
end

if exist('tSA') & ~isempty(tSA),
   decorrtime=1;
   if size(tSA,1)~=U*2-1,
      disp('Temporal dimension must match in SR and tSA!');
      return
   end
else
   decorrtime=0;
end

testcount=length(tolval);
H=zeros([N,U,testcount,s(3:end)]);

% as in XC, each response set is essentially independent here
for respidx=1:M,
   
   if M>1,
      if respidx==1,
         fprintf('respidx=(%d) 1',M);
      elseif respidx==M
         fprintf(' %d\n',respidx);
      else
         fprintf(' %d',respidx);
      end
   end
   
   % decorr time if tSA included
   if decorrtime,
      % normalize tSA to avoid double stimulus power over-compensation
      %m=mean(diag(sSA(:,:,respidx)));
      %tSA(:,respidx)=tSA(:,respidx)./m;
      m=max(tSA(:,respidx));
      sSA(:,:,respidx)=sSA(:,:,respidx)./m;
      
      ttSA=zeros(U);
      for uu=1:U,
         ttSA(:,uu)=tSA((U+1-uu):(U*2-uu),respidx);
      end
      
      % compute SVD on tSA and then pseudo-inverse
      [tU tS tV] = svd(ttSA);
      
      teigs=diag(tS)./sum(diag(tS));
      tSAinv=tV*diag(1./diag(tS).*(teigs>TMINSFS))*tU';
      
      tSR=SR(:,:,respidx)*tSAinv';
   else
      tSR=SR(:,:,respidx);
   end
   
   if respidx==1 | size(sSA,3)>1 | size(sSA,4)>1,
      [sU sS sV] = svd(sSA(:,:,respidx));
   end
   
   % spatial decorr
   % first precompute some intermediate matrices to speed things up
   
   % left side of decorr
   dsS=diag(sS);
   dsSi=zeros(size(dsS));
   dsSi(dsS>0)=1./dsS(dsS>0); % no power at these spatial
                              % dimension, set to 0
   dsSi(isinf(dsSi))=0;
   sVsSinv=sV*diag(dsSi);
   
   % right side:
   % temporalily transform kernel into eigenvector
   sUtSR=sU'*tSR;
   
   % row vector fraction of total eigenvalues
   seigs=dsS.^2./sum(dsS.^2);
   
   neigs=zeros(size(tolval));
   
   for tidx=1:testcount,
      neigs(tidx)=length(find(cumsum(seigs)<1-tolval(tidx)));
      
      if neigs(tidx)==0,
         H(:,:,tidx,respidx)=tSR;
      else
         % assume mineigs are ordered from highest to lowest...
         sUtSR0=sUtSR;
         
         sUtSR0(neigs(tidx)+1:end,:)=0;
         
         %fprintf('tidx=%d mineigs=%.4f count=%d\n',...
         %        tidx,mineigs(tidx),length(find(seigs>mineigs(tidx))));
         H(:,:,tidx,respidx)=sVsSinv * sUtSR0;
      end
   end
   
   % update queue if active
   dbsetqueue(BATQUEUEID,1);
end




