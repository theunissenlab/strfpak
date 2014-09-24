function [CS, errFlg] = cal_AutoCorr_for_cache(DS_data, stim_avg,twindow, nband, JN_flag)
%
%  [CS, t] = cal_AutoCorr(DS, stim_avg, twindow, nband, JN_flag)
%     --  Calucate auto correlation matrix of stimulus
%     Input:
%         DS_data: the data struct that contains four fields: 
%               stims  - stimulus file name
%               nlength    - length of time domain
%               ntrials    - num of trials
%              e.g. DS_data(1) = struct('stimfiles', stimulus,
%                'nlength', 1723, 'ntrials', 20);
%         stim_avg: avg stimulus that used to smooth the noise 
%                   If stim_avg is empty, we will call cal_AVG to get it.
%         twindow: the variable to set the time interval to calculate
%                  autocorrelation. e.g. twindow=[-TimeLag TimeLag]
%         nband: the size of spatio domain of the stimulus file
%         JN_flag: the flag that specify whether we calculate JackKnifed CS
%                  The default value of JN_flag = 0(dont calulate) 
%    Output:
%          CS: the autocorrelation matrix, its size is 
%              ((nband*(nband-1))/2+nband) X  (2*twindow +1)
%        
%             STRFPAK: STRF Estimation Software
% Copyright ©2003. The Regents of the University of California (Regents).
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

% Created by JXZ, 2002.
% Dec. 20, 2002 - new algorithm for calculating autocorrelation
%

% ========================================================
% check whether we have valid required input
% ========================================================
running_flag = 0;
errFlg = 0;
if isempty(DS_data)
    errordlg('ERROR: Please enter non-empty data filename',...
        'Input Data Error', 'modal')
    errFlg = 1;
    return
   
end

if ~exist('JN_flag')
    JN_flag = 0;
end

% ========================================================
% check whether stim_avg has been calculated or not
% ========================================================
% if isempty(stim_avg)
%     % calculate stim_avg and avg_psth and psth 
%     [stim_avg, avgr, psth] = cal_AVG(DS, nband);
% end

% ========================================================
% initialize the output and set its size
% ========================================================
% get the total data files
filecount = length(DS_data); 

% temporal axis range
tot_corr = diff(twindow) + 1;

% spatial axis range
spa_corr = (nband * (nband - 1))/2 + nband;

% initialize CS and CSJN variable
CS = zeros(spa_corr, tot_corr);
%CSJN = zeros(filecount, spa_corr, tot_corr);
CS_ns = zeros(1, tot_corr);
%CS_ns_JN = zeros(filecount, spa_corr, tot_corr);

% ========================================================
% do calculation. The algorithm is based on FET's dcp_stim.c
% ========================================================

%Visually give calculate auto_corr status
for fidx = 1:filecount     % loop through all data files
    
    % load stimulus file
    stim_env = DS_data(fidx).stims;
    
    nlen = size(stim_env,2);
    stimval = zeros(nband, nlen);
    
%     % Check if input data are chosen properly.
%     thisLength = size(stim_env, 2);
%     
%     if thisLength < nlen
%         answ= questdlg(['Data Error: Please check your input data by clicking ',...
%                 '"Get Files" Button in the main window: The first data file need ',...
%                 'to be stimuli and the second data file need to be its corresponding',...
%                 ' response file. If you made a mistake, please type "clear all" ',...
%                 ' or hit "reset" button first and then choose input data again.',...
%                 ' Otherwise, we will truncate the longer length. Do you want to continue?'],...
%             'Input Data Error','Yes','No','No');
%         switch answ
%             case 'No'
%                 errFlg = 1;
%                 return
%         end
%     end
    
%     nlen = min(nlen, thisLength);

    % subtract mean of stim from each stim 
    for tgood = 1:nlen
        stimval(:, tgood) = stim_env(1:nband, tgood) - stim_avg(1:nband);
    end
    
    % do autocorrelation calculation
    if running_flag == 1
        tempWait = waitbar(0,...
            sprintf('Calculating Auto-correlation for dataset %d',fidx));
    end
    CS_diff = do_cached_calc('small_autocorr',stimval,nband,DS_data(fidx).ntrials,size(CS),twindow(2));
    CS = CS + CS_diff;
    if running_flag == 1
        close(tempWait)
    end 
    
    % Count the total trials for later normalization 
    lengthVec = ones(1, nlen);
    CS_ns(1, :) = CS_ns(1, :) + DS_data(fidx).ntrials *...
        xcorr(lengthVec, lengthVec, twindow(2));
    
    % clear workspace
    clear stim_env
    clear stimval
    clear lengthVec
    
end		        % END of fidx


disp('Done auto-correlation  calculation');

% ========================================================
% To normalize CS by CS_ns: 
%   if CS_ns ~= 0
%      CS = CS /CS_ns
%   end
% ========================================================
% elminate zero in CS_ns 
nozero_ns = isinf( 1 ./CS_ns) + CS_ns;

% normalize CS matrix 
for i=1:spa_corr
    CS(i, :) = CS(i, :) ./ nozero_ns;
end

% ========================================================
% save stimulus auto-corrlation matrix into file 
% ========================================================

%  Now this section has been moved to calculation_GUI.m

% ========================================================
% END OF CAL_AUTOCORR
% ========================================================

