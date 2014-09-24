function [info, infoup, infodown, fpxy,cxy,cxyup, cxydown, cxy_notnormalized,cxy_notnormalizedup,cxy_notnormalizedlo, cxypsth,cxypsthup, cxypsthdown, infopsth, infopsthup, infopsthdown]=SNRinfo_no_cutoff(ntrials, flag, varargin)
%[info, infoup, infodown, fpxy,cxy,cxyup, cxydown, cxy_notnormalized,cxy_notnormalizedup,cxy_notnormalizedlo]=SNRinfo(ntrials, flag, varargin)
%infoup and down are upper and lower bounds.cxyup and down are upper and lower bounds for coherence.
%cxy_notnormalized means just the coherence without transforming to coherence per single trial.
%this is called by plot_est to get the coherence and SNR info between two halves of the psth (flag = 2)
%as well as the normalized SNR info between a mean rate and a single spike trial (flag = 1).
%if flag = 3, this is if the rate to cohere with is not the actual mean.
%in this flag = 3 case, the two halves of the psth must also be columns of x, along with the mean psth and the rate to be compared
%             STRFPAK: STRF Estimation Software
% Copyright �2003. The Regents of the University of California (Regents).
% All Rights Reserved.
% Created by Theunissen Lab and Gallant Lab, Department of Psychology, Un
% -iversity of California, Berkeley.
%
% Permission to use, copy, and modify this software and its documentation
% for educational, research, and not-for-profit purposes, without fee and
% without a signed licensing agreement, is hereby granted, provided that
% the above copyright notice, this paragraph and the following two paragr
% -aphs appear in all copies and modifications. Contact The Office of Tec
% -hnology Licensing, UC Berkeley, 2150 Shattuck Avenue, Suite 510,
% Berkeley, CA 94720-1620, (510) 643-7201, for commercial licensing
% opportunities.
%
%IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
%SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
%ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
%REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
%LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
%PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY,
%PROVIDED HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PRO
%-VIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

% Created and Modified by Anne and FET, 2003.
%
%

if ~(and(flag>=1, flag<=3))
    disp('flag must be one two or three');
    return
end
[x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers] = mtparam(varargin);
global running_in_script_mode
if 0%strcmp(running_in_script_mode,'yes') %This is the calculation that gets re-used, so you might want to cache this one anyways.
    [y, fpxy, cxyo, cxyo_u, cxyo_l stP]=mtchd_JN(x(:,1:2),nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers);

else
    [y, fpxy, cxyo, cxyo_u, cxyo_l stP]=do_locally_cached_calc(get_local_cache_dir,'mtchd_JN',x(:,1:2),nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers);
end

cxy_notnormalized=cxyo(:,1,2);
cxy_notnormalizedup=cxyo_u(:,1,2);
cxy_notnormalizedlo=cxyo_l(:,1,2);

if flag==3
    if size(x,2)~=4
        display('must have four columns for x, spikepre, meanpsth, psth1half, psth2half');
        return
    else
        if 0%strcmp(running_in_script_mode,'yes')
            [y, fpxy, cxyo, cxyo_u, cxyo_l stP]=mtchd_JN(x(:,3:4),nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers);

        else
            [y, fpxy, cxyo, cxyo_u, cxyo_l stP]=do_locally_cached_calc(get_local_cache_dir,'mtchd_JN',x(:,3:4),nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers);
        end
        cxy_notnormalizedpsth=cxyo(:,1,2);
        cxy_notnormalizedpsthup=cxyo_u(:,1,2);
        cxy_notnormalizedpsthlo=cxyo_l(:,1,2);

        %the following sets to zero all the values where lower bound <= 0.
        indexzero=find(cxy_notnormalizedpsthlo<=0);
        %         if ~isempty(indexzero)
        %             cxy_notnormalizedpsthup(indexzero)=0;
        %             cxy_notnormalizedpsth(indexzero)=0;
        %             cxy_notnormalizedpsthlo(indexzero)=0;
        %         end

        cxy_notnormalizedpsth=cxy_notnormalizedpsth.^2;
        cxy_notnormalizedpsthup=cxy_notnormalizedpsthup.^2;
        cxy_notnormalizedpsthlo=cxy_notnormalizedpsthlo.^2;
    end
end


%the following sets to zero all the values where lower bound <= 0.
% indexzero=find(cxy_notnormalizedlo<=0);
% if ~isempty(indexzero)
%     cxy_notnormalizedup(indexzero)=0;
%     cxy_notnormalized(indexzero)=0;
%     cxy_notnormalizedlo(indexzero)=0;
% end

cxy_notnormalized=cxy_notnormalized.^2;
cxy_notnormalizedup=cxy_notnormalizedup.^2;
cxy_notnormalizedlo=cxy_notnormalizedlo.^2;


cxydown=cxy_notnormalizedlo;
cxyup=cxy_notnormalizedup;
cxy=cxy_notnormalized;

%changes coherences to that of one trial (normalization) except for where cxy=0;
if flag ==2
    index=find(cxydown~=0);
    kdown=(-ntrials+ntrials*sqrt(1./cxydown(index)))/2;
    cxynewdown=1./(kdown+1);
    cxydown(index)=cxynewdown;
    clear cxynewdown;

    index=find(cxy~=0);
    k=(-ntrials+ntrials*sqrt(1./cxy(index)))/2;
    cxynew=1./(k+1);
    cxy(index)=cxynew;
    clear cxynew

    index=find(cxyup~=0);
    kup=(-ntrials+ntrials*sqrt(1./cxyup(index)))/2;
    cxynewup=1./(kup+1);
    cxyup(index)=cxynewup;
    clear cxynewup

elseif flag==3
    cxypsthdown=cxy_notnormalizedpsthlo;
    cxypsthup=cxy_notnormalizedpsthup;
    cxypsth=cxy_notnormalizedpsth;
    index=find(cxypsthdown~=0);
    kdown=(-ntrials+ntrials*sqrt(1./cxypsthdown(index)))/2;
    cxypsthnewdown=1./(kdown+1);
    cxypsthdown(index)=cxypsthnewdown;
    clear cxypsthnewdown;

    index=find(cxypsth~=0);
    k=(-ntrials+ntrials*sqrt(1./cxypsth(index)))/2;
    cxypsthnew=1./(k+1);
    cxypsth(index)=cxypsthnew;
    clear cxypsthnew

    index=find(cxypsthup~=0);
    kup=(-ntrials+ntrials*sqrt(1./cxypsthup(index)))/2;
    cxypsthnewup=1./(kup+1);
    cxypsthup(index)=cxypsthnewup;
    clear cxypsthnewup

    %now cxypsth's are the normalized coherence of one spike train with the actual meanrate.
    cxydown(index)=cxy_notnormalizedlo(index).*(1+sqrt(1./cxy_notnormalizedpsthlo(index)))./(-ntrials+ntrials*sqrt(1./cxy_notnormalizedpsthlo(index))+2);

    cxy(index)=cxy_notnormalized(index).*(1+sqrt(1./cxy_notnormalizedpsth(index)))./(-ntrials+ntrials*sqrt(1./cxy_notnormalizedpsth(index))+2);
    cxyup(index)=cxy_notnormalizedup(index).*(1+sqrt(1./cxy_notnormalizedpsthup(index)))./(-ntrials+ntrials*sqrt(1./cxy_notnormalizedpsthup(index))+2);


elseif flag == 1

    cxydown=cxydown./(cxydown + ntrials*(1-cxydown));
    cxy=cxy./(cxy + ntrials*(1-cxy));
    cxyup=cxyup./(cxyup + ntrials*(1-cxyup));

end


infodown=-fpxy(2)*sum(log2(1-cxydown));
infoup=-fpxy(2)*sum(log2(1-cxyup));
info=-fpxy(2)*sum(log2(1-cxy));

if flag==3
    infopsthdown=-fpxy(2)*sum(log2(1-cxypsthdown));
    infopsthup=-fpxy(2)*sum(log2(1-cxypsthup));
    infopsth=-fpxy(2)*sum(log2(1-cxypsth));

else
    infopsthdown=0;
    infopsthup=0;
    infopsth=0;
    cxypsth=0;
    cxypsthup=0;
    cxypsthdown=0;
end
