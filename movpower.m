% function fmoviep=movpower(movie,startidx,stopidx,bwindow,btakepower,.. 
%                           bzdc,btempnl); 
%
% transform a stimulus movie into the fourier power domain.
% 
% startidx,stopidx - start and stop frames (default 0,0, keep all)
% bwindow - apply hanning window before fft (default 0, off)
% btakepower - take fourier amplitude to this power (default 1)
% bzdc - subtract this value from movie before window, then add back
%        before fourier transform... unless btempnl~=0, in which case
%        negative component gets stuck in an weak sf channel (def 0)
% btemnpnl - experimental code for replacing high sf channels with
%            nl-transformed temporal information (default 0)
%
% CREATED SVD 4/25/02  Hacked from movphasesep.m
%
function fmoviep=movpower(movie,startidx,stopidx,bwindow,...
                          btakepower,bzdc,btempnl);

Xmax=size(movie,1);
Ymax=size(movie,2);
movlen=size(movie,3);

if not(exist('startidx','var')) | startidx < 1 | startidx > movlen,
   startidx=1;
end
if not(exist('stopidx','var')) | stopidx < 1 | stopidx > movlen,
   stopidx=movlen;
end
if not(exist('bwindow','var')),
   bwindow=0;
end
if not(exist('btakepower','var')),
   btakepower=1;
end
if not(exist('bzdc','var')),
   bzdc=0;
end
if not(exist('btempnl','var')),
   btempnl=0;
end
movlen=stopidx-startidx+1;

pixcount=Xmax*Ymax;
xc=round((Xmax+1)/2);
yc=round((Ymax+1)/2);

movie=movie(:,:,startidx:stopidx);

if bzdc & btempnl==1, 
   fprintf('movpower.m: bzdc=%.1f\n',bzdc);
   movie=movie-bzdc;
   bdc=0;
elseif bzdc,
   %fprintf('movpower.m: bzdc=%.1f\n',bzdc);
   movie=movie-bzdc;
   bdc=bzdc;
else
   bdc=0;
end
if bwindow,
   %filt=(hanning(Xmax,'periodic') * hanning(Ymax,'periodic')');
   filt=hanning2(Xmax,'periodic');
   movie=movie .* repmat(filt,[1 1 stopidx-startidx+1]) + bdc;
end


fmovie=fft(movie,[],1);
fmovie=cat(1,fmovie(xc:Xmax,:,:),fmovie(1:xc-1,:,:));
fmovie=fft(fmovie,[],2);
fmovie=cat(2,fmovie(:,yc:Ymax,:),fmovie(:,1:yc-1,:));

% if pixels are even in number, copy over extreme left pixels to
% imaginary part of real-only top row
if 0 & Xmax/2==round(Xmax/2),
   fmovie(1,2:Xmax,:)=fmovie(1,2:Xmax,:)+...
       i*imag(reshape(fmovie(2:Xmax,1,:),1,Xmax-1,movlen));
end

if 0,
   % do nothing. change to 1 to skip this use of btempnl
elseif btempnl==1,
   fprintf('movpower.m: btempnl=%.1f\n',btempnl);
   dc=fmovie(xc,yc,:);
   fmovie(xc,yc,:)=fmovie(xc,yc,:).*(fmovie(xc,yc,:)>0);
   fmovie(1,Ymax,:)=dc.*(dc<0);
elseif btempnl==2,
   fprintf('movpower.m: btempnl=%.1f\n',btempnl);
   dc=fmovie(xc,yc,:);
   fmovie(xc,yc,:)=0;
   
   % normalize by square root of mean power in each image. a crude
   % attempt at contrast gain control.
   vv=sqrt(mean(mean(abs(fmovie))));
   fmovie(xc,yc,:)=dc.*(dc>0);
   fmovie(1,Ymax,:)=dc.*(dc<0);
   fmovie=fmovie./repmat(vv+(vv==0),[Xmax Ymax]);
end

% reshape to space X time
fmovie=reshape(fmovie,pixcount,movlen);

% only save half of the bins since the other half are redundant
cfilt=gencfilt(Xmax,Ymax);

% take abs to get phase-free fft amplitude. optionally scale by
% taking a power of each coefficient.
if btakepower>0 & btakepower~=1,
   % channel value = abs(fft(movie)).^btakepower.  default is just 1.
   
   %fprintf('movpower.m: btakepower=%.1f\n',btakepower);
   fmoviep=abs(fmovie(cfilt,:)).^btakepower;
elseif btakepower<0,
   % nonlinear normalize by over contrast of image
   fprintf('movpower.m: btakepower=%.1f\n',btakepower);
   
   dcidx=find(cfilt==(xc-1).*Xmax+yc);
   fmoviep=abs(fmovie(cfilt,:)).^0.5;
   
   fms=std(fmoviep([1:dcidx-1 dcidx+1:end],:));
   
   if size(fmoviep,2)>10,
      fmsstd=std(fms);
   else
      % subst mean for std measure if there aren't that many frames
      fmsstd=mean(fms);
   end
   
   fms0=(tanh(fms./fmsstd.*abs(btakepower)).*fmsstd)./fms./(abs(btakepower));
   chancount=size(fmoviep,1);
   fmoviep=fmoviep.*repmat(fms0,[chancount,1]);
   
elseif btakepower<0,
   % take power but saturate if more than btakepower stds above
   % mean ... too kludgy?
   fprintf('movpower.m: btakepower=%.1f\n',btakepower);
   
   fmoviep=abs(fmovie(cfilt,:));
   if size(fmoviep,2)>10,
      fms=std(fmoviep,0,2);
   else
      % subst mean for std measure if there aren't that many frames
      fms=mean(fmoviep,2);
   end
   chancount=length(fms);
   for ii=1:chancount,
      %plot(fmoviep(ii,:),...
      %     tanh(fmoviep(ii,:)./fms(ii).*abs(btakepower)).*fms(ii),'.');
      fmoviep(ii,:)=tanh(fmoviep(ii,:)./fms(ii).*abs(btakepower)).*fms(ii);
   end
   
else
   fmoviep=abs(fmovie(cfilt,:));
end

% dcidx is the frame where the dc power is contained
dcidx=find(cfilt==(xc-1).*Xmax+yc);

if 0 & bzdc, % old SVD 4/19/04
   % normalize by dc of each frame
   dc=squeeze(fmoviep(dcidx,:));  % save DC for each frame
   dc(find(dc==0))=1;  % avoid division by zero, if it ever happens.
   
   fmoviep=fmoviep./repmat(dc, [length(cfilt) 1]);
   fmoviep(dcidx,:)=0;
end

if btempnl,
   
elseif 1 & btempnl,
   disp('gaussianizing!');
   chancount=size(fmoviep,1);
   stimcount=size(fmoviep,2);
   m=mean(fmoviep,2);
   v=std(fmoviep,0,2);
   for ii=1:chancount,
      outdata=sort(randn(1,stimcount,1)) .* v(ii)+m(ii);
      [aa,tdata]=sort(fmoviep(ii,:));
      fmoviep(ii,tdata)=outdata;
   end
   
elseif 1 & btempnl,
   
   % use first entries in pfft space--these are high freqs
   fprintf('<TNL>');
   
   for ii=1:abs(btempnl),
      
      % leave it squared for the time being
      if ii==1,
         %td=squeeze(mean(mean(movie(:,:,(ii+1):end) - ...
         %                     movie(:,:,1:end-ii),1),2).^2)';
         td=abs(fmoviep(dcidx,(ii+1):end)-fmoviep(dcidx,1:end-ii));
         td=[mean(td) td];
      elseif ii==2,
         td=abs(fmoviep(dcidx,3:end)-...
                0.5.*fmoviep(dcidx,2:end-1)-...
                0.5.*fmoviep(dcidx,1:end-2));
         td=[mean(td) mean(td) td];
      elseif ii==3,
         td=abs(0.5.*fmoviep(dcidx,3:end)+...
                0.5.*fmoviep(dcidx,2:end-1)-...
                fmoviep(dcidx,1:end-2));
         td=[mean(td) mean(td) td];
      end
      
      if btempnl>0,
         fmoviep(ii,:)=td./abs(btempnl.*6);
      else
         fmoviep(ii,:)=shuffle(td)./abs(btempnl.*6);
      end
   end
   
elseif 0 & btempnl,
   
   fprintf('<TNL2>');

   threshold=0.5;
   [firstbins,lastbins]=movgetepochs(movie,threshold);
   if firstbins(1)==1,
      firstbins=firstbins(2:end);
   end
   fmoviep(1:2,:)=0;
   fmoviep(1,firstbins)=abs(fmoviep(dcidx,firstbins)-fmoviep(dcidx,firstbins-1));
   fmoviep(2,lastbins)=abs(fmoviep(dcidx,lastbins)-fmoviep(dcidx,lastbins-1));
   
end




