% function [SR,mS,mR,tSA,sSA2]=normalize0(SR,n,mS,mR,tSA,sSA2,meansub)
%
% created SVD 3/26/02
% modified SVD 7/12/02 - removed sSA1, added mS,mR
%
% Inputs (basically the output from movxc.m or summed output from
%        multiple runs of it):
%      N - number of spatial dimensions
%      U - total number of time lags (eg, length([maxlag(1):maxlag(2)]))
%      M - number of response sets (or response dimensions)
%      R - number of resamples
% SR - raw spike-triggered average (correlation between stim and
%      resp) (N x U x M x R)
% n   - sample counts (M x R)
% mS  - stimulus mean (N x M x R)
% mR  - response mean (M x R)
% tSA - temporal autocorrelation matrix (U x U x M x R or leave empty
%       [], to skip temporal decorrelation)
% sSA2 - spatial autocorrelation matrix (N x N x M x R or N x N and use
%       the same one for each response set--speeds things up)
%
function [SR,mS,mR,tSA,sSA2]=normalize0(SR,n,mS,mR,tSA,sSA2,meansub)

N=size(SR,1);
M=size(n,1);
R=size(n,2);
A=size(n,3);

if ~exist('meansub'),
   meansub=1;
end

for ii=1:M,
   for kk=1:A,
      
      if meansub,
         tSall=sum(mS(:,ii,:,kk),3)./sum(n(ii,:,kk));
         tRall=ones(1,size(SR,2))*sum(mR(ii,:,kk))./sum(n(ii,:,kk));
      else
         tSall=0;
         tRall=0;
      end
      
      for jj=1:R,
         % normalize the means
         mS(:,ii,jj,kk)=mS(:,ii,jj,kk)./n(ii,jj,kk);
         mR(ii,jj,kk)=mR(ii,jj,kk)./n(ii,jj,kk);
         if meansub,
            tS=mS(:,ii,jj,kk);
            tR=ones(1,size(SR,2))*mR(ii,jj,kk);
         else
            tS=0;
            tR=0;
         end
         
         %SR(:,:,ii,jj,kk)=SR(:,:,ii,jj,kk)./n(ii,jj,kk) - tSall*tRall;
         SR(:,:,ii,jj,kk)=SR(:,:,ii,jj,kk)./n(ii,jj,kk) - tS*tR;
         
         if ii<=size(tSA,2) & jj<=size(tSA,3),
            tSA(:,ii,jj,kk)=tSA(:,ii,jj,kk)./n(ii,jj,kk) - tS'*tS./N;
            if exist('sSA2'),
               sSA2(:,:,ii,jj,kk)=sSA2(:,:,ii,jj,kk)./n(ii,jj,kk) - tS * tS';
            else
               sSA2=[];
            end
         end
      end
   end
end

      
