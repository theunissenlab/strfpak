% [mH,eH,mSR,eSR,mS,mR]=xcspacetimesep(stim,resp,params);
% 
% function to compute unbiased kernel from preprocessed stimulus and
% response matrices (ie, stim T x X and resp T x M supplied or
% stim X x T and resp M x T)
% 
% returns:
% mH, eH: mean and stderr on kernel estimates (X x tbincount x M x sfscount)
% mSR, eSR: mean and stderr on raw (biased) XC (X x tbincount x M x sfscount)
% mS, mR: mean stim/response across resamples
% 
% params (defaults):
% params.resampcount=getparm(params,'resampcount',20);
% params.decorrspace=getparm(params,'decorrspace',2);
% params.decorrtime=getparm(params,'decorrtime',1);
% params.smoothtime=getparm(params,'smoothtime',2);
% params.maxlag=getparm(params,'maxlag',15);
% params.tolval=getparm(params,'tolval',...
%                       [0.1 0.05 0.01 0.005 0.001 0.0005 0.0001]);
% 
% required programs:
% normalize0.m, normalize.m, resampsegs.m, getparm.m
% 
% created SVD 6/12/03 - ripped off of xccore.m
% 
function [mH,eH,mSR,eSR,mS,mR]=xcspacetimesep(stim,resp,params);

%
% set params to default values if they don't exist
%
if ~exist('params','var'),
   disp('using default params.');
   params=[];
end
params.resampcount=getparm(params,'resampcount',1);
params.decorrspace=getparm(params,'decorrspace',2);
params.decorrtime=getparm(params,'decorrtime',1);
params.smoothtime=getparm(params,'smoothtime',2);
if isfield(params,'maxlag') & length(params.maxlag)>=2,
   % do nothing, this is a good format for running cellxc
else
   params.maxlag=getparm(params,'maxlag',300);
   params.maxlag=[-params.maxlag params.maxlag];
end
params.tolval=getparm(params,'tolval',...
                             [0.1 0.05 0.01 0.005 0.001 0.0005 0.0001 0.00005]);


global BATQUEUEID

%
% make sure stim is T x N and resp is T x respcount
%
if size(stim,2)>size(stim,1),
   stim=stim';
end
if size(resp,2)>size(resp,1),
   resp=resp';
end

%
% define resampling regimes. if resampcount is 1, just make one regime.
%
if params.resampcount>1,
   [rstartidx,rendidx]=resampsegs(resp,params.resampcount);
else
   rstartidx=1;
   rendidx=size(resp,1);
end

disp('initializing matrices');
respcount=size(resp,2);
corrmtxcount=respcount;
singlesSA=0;
spacecount=size(stim,2);

% first time, set kernel & ac matrices to zero
SR=zeros(spacecount,diff(params.maxlag)+1,respcount,params.resampcount);
n=zeros(respcount,params.resampcount);
mS=zeros(spacecount,respcount,params.resampcount);
mR=zeros(respcount,params.resampcount);
tSA=zeros(diff(params.maxlag)*2+1,corrmtxcount,params.resampcount);
sSA2=zeros(spacecount,spacecount,corrmtxcount,params.resampcount);

% global stimulus temporal autocorr. works better than resamped?
tSA1=zeros(diff(params.maxlag)*2+1,corrmtxcount,params.resampcount);
nSA1=0;
firstseg=0;

fprintf('resamp=');
for resampidx=1:params.resampcount,
   fprintf('%d ',resampidx);
   for respidx=1:respcount,
      tr=resp(:,respidx);
      tr([1:rstartidx(resampidx,respidx)-1 ...
          rendidx(resampidx,respidx)+1:end],:)=nan;
      %tr(rstartidx(resampidx,respidx):rendidx(resampidx,respidx),:)=nan;
      
      rgoodidx=find(~isnan(tr));
      tmR=sum(resp(rgoodidx,respidx))';
      tmS=sum(stim(rgoodidx,:))';
      tn=length(rgoodidx);
      
      tSR=zeros(spacecount,diff(params.maxlag)+1);
      for tt=params.maxlag(1):params.maxlag(2),
         %fprintf('.');
         trg=rgoodidx;
         trg=trg(find(trg-tt>0 & trg-tt<size(stim,1)));
         tSR(:,tt-params.maxlag(1)+1)=(resp(trg,respidx)'*stim(trg-tt,:))';
      end
      
      if respidx<=corrmtxcount & length(rgoodidx)>0,
         if params.decorrspace==2,
            tsSA2=stim(rgoodidx,:)'*stim(rgoodidx,:);
            %fprintf('*');
         end
         ttSA=zeros(diff(params.maxlag)*2+1,1);
         for xx=1:spacecount,
            tstim=stim(rgoodidx,xx);
            ttSA=ttSA+xcorr(tstim,diff(params.maxlag),'biased') ./ ...
                 spacecount.*length(rgoodidx);
         end
         %fprintf('*\n');
      end
      
      for rr=1:params.resampcount,
         % add means to all resamp channels. this is to
         % make it so that everything is centered around
         % the same DC!
         if params.resampcount==1 | rr~=resampidx,
            % otherwise add outputs to running total
            SR(:,:,respidx,rr)=SR(:,:,respidx,rr)+tSR;
            mS(:,respidx,rr)=mS(:,respidx,rr)+tmS;
            mR(respidx,rr)=mR(respidx,rr)+tmR;
            n(respidx,rr)=n(respidx,rr)+tn;
            tSA(:,respidx,rr)=tSA(:,respidx,rr)+ttSA;
            
            if params.decorrspace==2 & respidx<=corrmtxcount,
               sSA2(:,:,respidx,rr)=sSA2(:,:,respidx,rr)+tsSA2;
            end
         end
      end
   end
   
   % update queue if active
   dbsetqueue(BATQUEUEID,1);
end
fprintf('\n');

for xx=1:spacecount,
   
   tstim=stim(:,xx)-mean(stim(:,xx));
   tSA1=tSA1+repmat(xcorr(tstim,diff(params.maxlag),'biased')./ ...
                    spacecount.*size(stim,1),...
                    [1 corrmtxcount params.resampcount]);
end
nSA1=nSA1+size(stim,1);

% zeroth order normalization. ie, divide sums by the number of
% samples to get appropriate means
[SR,mS,mR,tSA,sSA2]=normalize0(SR,n,mS,mR,tSA,sSA2);
tSA1=tSA1./nSA1;

fprintf('params.decorrspace=%d\n',params.decorrspace);
if params.decorrspace,
   [H,neigs]=normalize(SR,sSA2,tSA,params.tolval);
else
   H=reshape(SR,spacecount,diff(params.maxlag)+1,1,respcount,params.resampcount);
end

% alternative: regularized decorrelation. doesn't work with tolvals
% defined as percent variance to exclude
%[H,lambda]=normalizereg(SR,sSA2,tSA1,params.sfscount,params.sfsstep,...
%                        1,params.smoothtime); % topSR=1


mH=mean(H,5);
eH=std(H,1,5) .* sqrt((params.resampcount-1)/params.resampcount) .* ...
   sqrt(params.resampcount);

mSR=mean(SR,4);
eSR=std(SR,1,4).* sqrt((params.resampcount-1)/params.resampcount) .* ...
    sqrt(params.resampcount);

% take mean stim/response across resamples
mS=mean(mS,3);
mR=mean(mR,2);
ntot=n;

