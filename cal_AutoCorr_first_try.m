function [CS, errFlg] = cal_AutoCorr(running_flag, DS, stim_avg,twindow, nband, JN_flag)
%
%  [CS, t] = cal_AutoCorr(DS, stim_avg, twindow, nband, JN_flag)
%     --  Calucate auto correlation matrix of stimulus
%     Input:
%         DS: the data struct that contains four fields: 
%               stimfiles  - stimulus file name
%               respfiles  - response file name
%               nlength    - length of time domain
%               ntrials    - num of trials
%              e.g. DS{1} = struct('stimfiles', 'stim1.dat', 'respfiles',
%                    'resp1.dat', 'nlength', 1723, 'ntrials', 20);
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
errFlg = 0;
if isempty(DS)
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
if isempty(stim_avg)
    % calculate stim_avg and avg_psth and psth 
    [stim_avg, avgr, psth] = cal_AVG(DS, nband);
end

% ========================================================
% initialize the output and set its size
% ========================================================
% get the total data files
filecount = length(DS); 

% temporal axis range
tot_corr = diff(twindow) + 1;

% spatial axis range
spa_corr = (nband * (nband - 1))/2 + nband;

% initialize CS and CSJN variable
CS_size = [spa_corr, tot_corr];
%CSJN = zeros(filecount, spa_corr, tot_corr);
CS_ns_size = [1, tot_corr];

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
    stim_env = Check_And_Load(DS{fidx}.stimfiles);
    [CS_diff,CS_ns_diff] = one_AutoCorr(stim_env,nband,nlen,stim_avg,DS{fidx}.ntrials);
    CS = CS + CS_diff;
    CS_ns = CS_ns + CS_ns_diff
    clear stim_env    
end		        % END of fidx

clear CS_diff CS_ns_diff
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
currentPath = pwd;
global outputPath
if ~isempty(outputPath)
    cd (outputPath);
else
    disp('Saving output to Output Dir.');
    stat = mkdir('Output');
    cd('Output');
    outputPath = pwd;
end

save('Stim_autocorr.mat', 'CS'); 
cd(currentPath);
% ========================================================
% END OF CAL_AUTOCORR
% ========================================================
