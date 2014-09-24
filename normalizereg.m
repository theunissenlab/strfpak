% function [H,lambda]=normalizereg(SR,sSA,tSA,nfactors,lambdaomag,topSR,smoothtime)
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
function [H,lambda]=normalizereg(SR,sSA,tSA,nfactors,lambdaomag,topSR,smoothtime);

disp('normalizereg.m:');
global BATQUEUEID

TMINSFS=0.0001;  % hard coded minimum eigenvector for temporal
                 % decorr. hacky but it seems to work.

s=size(SR);
N=s(1);
U=s(2);
M=prod(s(3:end));
if ~exist('topSR','var'),
   topSR=1;
end

if size(sSA,1)~=N,
   disp('Spatial dimension must match in SR and sSA!');
   return
end

if exist('tSA','var') & ~isempty(tSA),
   decorrtime=1;
   if size(tSA,1)~=U*2-1,
      disp('Temporal dimension must match in SR and tSA!');
      return
   end
   if ~exist('smoothtime','var'),
      smoothtime=1;
   end
else
   decorrtime=0;
   if ~exist('smoothtime','var'),
      smoothtime=0;
   end
end
if ~exist('lambdaomag','var') | isempty(lambdaomag),
   lambdaomag=4;
end

if smoothtime>0,
   fprintf('(smoothtime=%d)',smoothtime');
end

H=zeros([N,U,nfactors,s(3:end)]);

% as in XC, each response set is essentially independent here
fprintf('respidx(/%d)=',M);
for respidx=1:M,
   
   fprintf('%d ',respidx);
   
   % decorr time if tSA included
   if decorrtime,
      
      % normalize sSA to avoid double stimulus power
      % over-compensation. this is derivable from the idea that
      % tsSA(x1,t1,x2,t2)=sSA(x1,x2) * tSA(t1,t2)
      m=max(tSA(:,respidx));
      
      % create at toeplitz tSA matrix to get full temporal
      % decorr. seems kind of hokey but avoids edge effects for
      % short kernels
      ttSA=zeros(U);
      for uu=1:U,
         ttSA(:,uu)=tSA((U+1-uu):(U*2-uu),1);
         %ttSA(:,uu)=mean(tSA((U+1-uu):(U*2-uu),:),2);
      end
      
      % compute SVD on tSA and then pseudo-inverse
      [tU tS tV] = svd(ttSA);
      teigs=diag(tS)./sum(diag(tS));
      tSAinv=tV*diag(1./diag(tS).*(teigs>TMINSFS))*tU';
      
      tSR=SR(:,:,respidx)*tSAinv';
      if smoothtime==1,
         tSR=conv2(tSR,[0.2 0.6 0.2],'same');
      end
   else
      tSR=SR(:,:,respidx);
      m=1;
   end
   
   % spatial decorr
   % first precompute some intermediate matrices to speed things up
   % do some preliminary stuff to speed up decorrelation algebra
   
   % do SVD and parse out zero and non-zero eigenvalues
   if respidx==1 | size(sSA,3)>1 | size(sSA,4)>1,
      %[sU sS sV] = svd(sSA(:,:,respidx)./m);
      [sU sS sV] = svd(sSA(:,:,respidx));
      sD=diag(sS);
      sDnz=sD(find(sD>0));
      sDz=sD(find(sD<=0));
   end
   
   % convert sta into svd domain.
   tSR=sU'*tSR;
   
   if ~exist('nfactors','var'),
      nfactors=1;
      lambda=0;
   else
      % scale reg parm to match stats of AC power
      % this case: first term scales so that first ac is equiv
      % to third???
      
      if topSR,
         lambda=sD(1).*sD(3).*10.^(-[-1 linspace(0,lambdaomag,nfactors-1)]);
      else
         % this case: first term scales so that first ac is equiv
         % to third???
         lambda=sD(1).*sD(3).*10.^(-linspace(0,lambdaomag,nfactors));
         
         % alt: make first reg parm a bit further out
         %lambda=sqrt(prod(sD(1:4)).*10.^(-linspace(0,4,nfactors));
      end
      
      %fprintf(' (lambda(2)=%.2f) sfsidx=',lambda(2));
   end
   
   for sfsidx=1:nfactors,
      
      % hackish thing: first guess is STA w/o spatial correlations removed
      if lambda(sfsidx)>=lambda(1) & topSR,
         H(:,:,sfsidx,respidx)=sV*tSR;
         %fprintf('(topSR) ');
      else
         % scale by regularization paramater here.
         sDi=1./(sDnz+lambda(sfsidx)./sDnz);
         
         % option: force max sDi >= sDi(1)
         sDi(find(sDi<sDi(1)))=sDi(1);
         
         SAspaceinv=sV*diag([sDi;sDz]);
         H(:,:,sfsidx,respidx)=SAspaceinv*tSR;
         if smoothtime==2,
            H(:,:,sfsidx,respidx)=conv2(H(:,:,sfsidx,respidx),...
                                        [0.2 0.6 0.2],'same');
         end
      end
      
      if mod(sfsidx,10)==0,
         %fprintf('%d . ',sfsidx);
         if exist('BATQUEUEID') & BATQUEUEID>0,
            dbsetqueue(BATQUEUEID,1);
         end
      end
   end
   
   %fprintf('\n');
end
fprintf('\n');




