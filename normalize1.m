% function [H,neigs]=normalize1(SR,sSA1,tSA,mineigs,neigs)
%
% created SVD 6/5/02 - hacked from normalize.m
%
% Inputs (basically the output from movxc.m or summed output from
%        multiple runs of it):
%      N - number of spatial dimensions
%      U - total number of time lags (eg, length([maxlag(1):maxlag(2)]))
%      M - number of response sets (or response dimensions)
% SR - raw spike-triggered average (correlation between stim and
%      resp) (N x U x M)
% sSA1 - spatial autocorrelation matrix (N x M or N x 1 and use
%       the same one for each response set--speeds things up)
% tSA - temporal autocorrelation matrix (U x U x M or leave empty
%       [], to skip temporal decorrelation)
%
% Outputs:
% H - decorrelated kernel estimate(s), N x U x length(mineigs) x M matrix
% neigs - number of eigenvectors included in corresponding
%         decorrelated kernel estimate
%
function [H,neigs,lambda]=normalize1(SR,sSA1,tSA,mineigs,neigs);

disp('normalize1.m:');
global BATQUEUEID

TMINSFS=0.0001;

s=size(SR);
N=s(1);
U=s(2);
M=prod(s(3:end));
if ~exist('neigs'),
   testcount=length(mineigs);
   useneigs=0;
   neigs=zeros(testcount,1);
else
   testcount=length(neigs);
   useneigs=1;
end

lambda=10.^linspace(-3.5,3.5,length(neigs));
%lambda=zeros(testcount,1);

if size(sSA1,1)~=N,
   disp('Spatial dimension must match in SR and sSA1!');
   return
end

if exist('tSA') & ~isempty(tSA),
   decorrtime=1;
   if size(tSA,1)~=U*2-1,
      disp('Temporal dimension must match in SR and tSA!');
      return
   end
   TSMOOTH=2;
else
   decorrtime=0;
   TSMOOTH=0;
end
if TSMOOTH>0,
   fprintf('TSMOOTH=%d\n',TSMOOTH');
end

H=zeros([N,U,testcount,s(3:end)]);

% as in XC, each response set is essentially independent here
for respidx=1:M,
   
   fprintf('respidx=%d/%d\n',respidx,M);
   
   % decorr time if tSA included
   if decorrtime,
      % normalize tSA to avoid double stimulus power over-compensation
      m=max(tSA(:,respidx));
      sSA1(:,respidx)=sSA1(:,respidx)./m;
      
      ttSA=zeros(U);
      for uu=1:U,
         ttSA(:,uu)=tSA((U+1-uu):(U*2-uu),respidx);
      end
      
      % compute SVD on tSA and then pseudo-inverse
      [tU tS tV] = svd(ttSA);
      
      teigs=diag(tS)./sum(diag(tS));
      tSAinv=tV*diag(1./diag(tS).*(teigs>TMINSFS))*tU';
      
      tSR=SR(:,:,respidx)*tSAinv';
      if TSMOOTH==1,
         tSR=conv2(tSR,[0.2 0.6 0.2],'same');
      end
   else
      tSR=SR(:,:,respidx);
   end
   
   if respidx==1 | size(sSA1,3)>1 | size(sSA1,4)>1,
      sP = sSA1(:,respidx);
   end
   
   %keyboard
   
   [sPsort,sPsortidx]=sort(sP);
   sPsortidx=flipud(sPsortidx);
   
   dPi=zeros(size(sP));
   dPi(find(sP>0))=1./sP(find(sP>0));
   dPi=repmat(dPi,[1 U]);
   neigs(find(neigs>N))=N;
   
   for tidx=1:testcount,
      if neigs(tidx)==0,
         H(:,:,tidx,respidx)=tSR;
         lambda(tidx)=0;
      else
         
         % alternative regularization strategry: neigs(tidx) is scaled to
         % 1/2 its value
         val1=sP(sPsortidx(1));
         val2=sP(sPsortidx(neigs(tidx)));
         z=0.5./(-1./val1^2 + 0.5./val2.^2);
         %keyboard
         %lambda(tidx)=z;
         z=lambda(tidx);
         
         dPi=zeros(size(sP));
         dPi(find(sP>0))=1./(sP(find(sP>0)) + z./sP(find(sP>0)));
         dPi=repmat(dPi,[1 U]);
         
         H(:,:,tidx,respidx)=tSR .* dPi;
         
         % no regularization. just do sharp cutoff
         %H(sPsortidx(1:neigs(tidx)),:,tidx,respidx)=...
         %    tSR(sPsortidx(1:neigs(tidx)),:) .* ...
         %    dPi(sPsortidx(1:neigs(tidx)),:);
      end
      if TSMOOTH==2,
         H(:,:,tidx,respidx)=conv2(H(:,:,tidx,respidx),...
                                   [0.2 0.6 0.2 ],'same');
      end
   end
   
   % update queue if active
   if exist('BATQUEUEID') & BATQUEUEID>0,
      dbsetqueue(BATQUEUEID,1);
   end
end




