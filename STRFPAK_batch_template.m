% 
%   STRFPAK's batch processing Template
%
%        - You need modify this file to meet your specific situation.
%          For example, you need assign data_path and change the parameters
%          you want to use.
%  
%  Created by Junli Zhang, 1/23/2004 
%  Updated by JXZ, 7/31/04
%      1. add computer OS checking
%      2. modify STRFPAK_core by adding preprocessOption, filteroption
%  JXZ: 9/23/2004
%      Add more comments.
%  JXZ: 9/07/2005
%      Add more calculation parameters for all options.

% ---------------------------------------------------
%  Clear everything in matlab memory and close all figures
% ---------------------------------------------------
clear all
close all

% ---------------------------------------------------
%  Add the current directory to matlab path
% ---------------------------------------------------
addpath(pwd);
current_CD = pwd;

% ---------------------------------------------------
%  NOTE: PLEASE SPECIFY THE PATHES OF YOUR DATA.
%        Since lower version of Matlab can not identify path 
%        separator, we check computer's OS and let you
%        beware when you assign data_path. 
% ---------------------------------------------------
% ---------------------------------------------------
% NOTE: ThIS IS ONLY FOR SPECIFIC LAYOUT. 
%      /DATA/ANIMALS/CELLS/STIMTYPES
% ---------------------------------------------------

if strcmp(computer, 'PCWIN')
    
    animal_list = {...
            'c:\jxz\strfpak\STRFPAK-3.0\DemoData\RawAuditoryData\Bird1',...
        };
    % ---------------------------------------------------
    % NOTE: PLEASE SPECIFY WHERE YOU WANT TO SAVE YOUR RESULT 
    % ---------------------------------------------------
    output_Path = 'c:\jxz\strfpak\strfpak-3.0\Output';
   
else
    animal_list = {...
            '/auto/fdata/fet/calibration_data/blabla0713',...
            %             '/auto/fdata/fet/calibration_data/blabla0903',...
        %             '/auto/fdata/fet/calibration_data/blabla0313',...
        %             '/auto/fdata/fet/calibration_data/blabla1515',...
        %             '/auto/fdata/fet/calibration_data/blabla1819',...
    };
    % ---------------------------------------------------
    % NOTE: PLEASE SPECIFY WHERE YOU WANT TO SAVE YOUR RESULT 
    % ---------------------------------------------------
    output_Path = '/auto/fdata/junli/STRFPAK_Result';
    
end


cell_count = 0;
for animal_indx=1:length(animal_list)
    
     cell_list = dir(animal_list{animal_indx});

     for cell_indx =1:length(cell_list)
         if strcmp(cell_list(cell_indx).name, '.')
             continue;
         elseif strcmp(cell_list(cell_indx).name, '..')
             continue;
         else
               
             cell_count = cell_count +1;
             data_files{cell_count}= fullfile(animal_list{animal_indx}, cell_list(cell_indx).name);
         end
    end
end


% ---------------------------------------------------
%  NOTE: IF YOUR COMPUTING FACILITY SUPPORTS PARALLEL COMPUTING,
%        JUST SET PARALLEL_OK = 1, 0 OTHERWISE.
% ---------------------------------------------------
PARALLEL_OK = 0;


% ---------------------------------------------------
% NOTE: Please specify which preprocessoption you want to use 
%  - preprocess: the preprocess option you choose:
%                1 : spectrogram (STFT)
%                2 : scalogram (WAVELET)
% ---------------------------------------------------
nfiles = length(data_files);
for preprocess = 1:1
    
    for ifiles=1:nfiles
        
        if PARALLEL_OK == 0
            
            if preprocess == 1  % Spectrogram 
                
                % ---------------------------------------------------
                % NOTE: You can modify the following proprocessing parameters:
                %  - fwidthHz : the width of the filter in Hz. 
                %               It defines the window length of the filter
                %  - respsamprate: sampling rate of your spike data in Hz
                %  - ampsamprate: the sampling rate of the amplitude evelope
                %                  in Hz.    
                %  - filteroption: the choice of scaling amplitude
                %                 1:1 - choose linear
                %                 2:2 - choose log
                %                 1:2 - both scale 
                %  - psth_smoothconst: window size of smoothing psth(in ms)
                %  - initialFreq: lower bound of frequency limit
                %  - endFreq:  upper bound of frequency you want to study
                %  - sDim:  dimension size of stimulus except time domain
                % ---------------------------------------------------
                % NOTE: You can modify the following calculation parameters:
                % - TimeLag: time lag used for calculating the STA.
                % - Tol_val: a regularization parameter used to estimate
                % STRF. 
                % - Std_val: another regularization parameter in the STRF
                % estimation that enforces sparseness.
                % - timevary_PSTH: flag for whether removing time-varying mean firing
                % rate from inputs.
                % - smooth_rt: window size for smoothing time-varying mean
                % firing rate.
                % - smoothVect: a vector of smoothing window range for
                % validation stage.
                % ---------------------------------------------------
                fwidthHz_list = [125];
                for ii = 1:length(fwidthHz_list)
                    
                    cd(data_files{ifiles});
                    outputPath_temp = output_Path;
                    preprocessoption = 1;
                    clear global respsamprate fwidthHz ampsamprate initialFreq endFreq sDim filteroption
                    clear global  TimeLag Tol_val Std_val timevary_PSTH smooth_rt smoothVect
                    global initialFreq endFreq respsamprate fwidthHz ampsamprate sDim filteroption
                    global  TimeLag Tol_val Std_val timevary_PSTH sDim smooth_rt smoothVect
                    
                    % preprocessing parameters
                    sDim = '1-D';
                    respsamprate = 1000;
                    fwidthHz = fwidthHz_list(ii);
                    ampsamprate = 1000;
                    filteroption = 2:2;
                    psth_smoothconst = 21;
                    initialFreq = 250;
                    endFreq = 8000;
                    
                    % calculation parameters from users
                    TimeLag = 300;
                    Tol_val = [...
                            0.100000 0.050000 0.010000 0.005000 0.001000 0.000500 0.000100 0.000050 ];
                    Std_val = [0 0.25 0.5 0.75 1 1.5 2 4 6];
                    timevary_PSTH = 1;
                    smooth_rt = 41;
                    smoothVect = [5 8 29];
                    STRFPAK_core;
                end
                
                
            elseif preprocess == 2   % Wavelet Transform
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % You can modify the following parameters:
                %  - wavelength:  It defines the window length of the wavelet 
                %  - NBAND:  the number of frequency samp points 
                %  - respsamprate: sampling rate of your spike data in Hz
                %  - ampsamprate: the sampling rate of the amplitude evelope
                %                  in Hz.    
                %  - filteroption: the choice of scaling amplitude
                %                 1:1 - choose linear
                %                 2:2 - choose log
                %                 1:2 - both scale 
                %  - psth_smoothconst: window size of smoothing psth(in ms)
                 % ---------------------------------------------------
                % NOTE: You can modify the following calculation parameters:
                % - TimeLag: time lag used for calculating the STA.
                % - Tol_val: a regularization parameter used to estimate
                % STRF. 
                % - Std_val: another regularization parameter in the STRF
                % estimation that enforces sparseness.
                % - timevary_PSTH: flag for whether removing time-varying mean firing
                % rate from inputs.
                % - smooth_rt: window size for smoothing time-varying mean
                % firing rate.
                % - smoothVect: a vector of smoothing window range for
                % validation stage.
                % ---------------------------------------------------
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Nfilt_list = [21 41 81];
                for ii =1:length(Nfilt_list)
                    cd(data_files{ifiles});
                    outputPath_temp = output_Path;
                    preprocessoption = 2;
                    clear global respsamprate fwidthHz ampsamprate initialFreq endFreq sDim filteroption
                    global respsamprate fwidthHz ampsamprate wavelength Npoints initialFreq endFreq sDim filteroption
                    clear global smooth_rt TimeLag Tol_val Std_val timevary_PSTH smoothVect
                    global timevary_PSTH smooth_rt TimeLag Tol_val Std_val smoothVect
                    
                    % preprocessing parameters
                    wavelength = 20;
                    respsamprate = 1000;
                    fwidthHz = 125;
                    Npoints = Nfilt_list(ii);
                    ampsamprate = 1000;
                    filteroption = 2:2;
                    psth_smoothconst = 21;
                    initialFreq = 250;
                    endFreq = 8000;
                    sDim = '1-D';
                    
                    % calculation parameters from users
                    TimeLag = 300;
                    Tol_val = [...
                            0.100000 0.050000 0.010000 0.005000 0.001000 0.000500 0.000100 0.000050 ];
                    Std_val = [0 0.25 0.5 0.75 1 1.5 2 4 6];
                    timevary_PSTH = 1;
                    smooth_rt = 41;
                    smoothVect = [5 8 29];
                    
                    STRFPAK_core;
                    
                end   % END of Nfilt_list
            end     % END of preprocess    
            
        elseif PARALLEL_OK == 1
            if preprocess == 1  % Spectrogram 
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % You can modify the following parameters:
                %  - fwidthHz : the width of the filter in Hz. 
                %               It defines the window length of the filter
                %  - respsamprate: sampling rate of your spike data in Hz
                %  - ampsamprate: the sampling rate of the amplitude evelope
                %                  in Hz.    
                %  - filteroption: the choice of scaling amplitude
                %                 1:1 - choose linear
                %                 2:2 - choose log
                %                 1:2 - both scale 
                %  - psth_smoothconst: window size of smoothing psth(in ms)
                 % ---------------------------------------------------
                % NOTE: You can modify the following calculation parameters:
                % - TimeLag: time lag used for calculating the STA.
                % - Tol_val: a regularization parameter used to estimate
                % STRF. 
                % - Std_val: another regularization parameter in the STRF
                % estimation that enforces sparseness.
                % - timevary_PSTH: flag for whether removing time-varying mean firing
                % rate from inputs.
                % - smooth_rt: window size for smoothing time-varying mean
                % firing rate.
                % - smoothVect: a vector of smoothing window range for
                % validation stage.
                % ---------------------------------------------------
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fwidthHz_list = [250];
                for ii = 1:length(fwidthHz_list)
                    
                    dbaddqueuemaster([' cd ',data_files{ifiles},';' ,...
                            ' outputPath_temp = ''', output_Path, ''';', ...
                            ' preprocessoption = 1;', ...
                            ' clear global respsamprate fwidthHz ampsamprate initialFreq endFreq sDim filteroption; ', ...
                            ' global respsamprate fwidthHz ampsamprate initialFreq endFreq sDim filteroption;',...
                            ' respsamprate = 1000; ', ...
                            ' fwidthHz = ', num2str(fwidthHz_list(ii)),'; ',...
                            ' ampsamprate = 1000; ',...
                            ' filteroption = 2:2;  ', ...
                            ' initialFreq = 250;   ',...
                            ' endFreq = 8000;    ',...
                            ' sDim = \''1-D\'';    ',...
                            ' clear global smooth_rt TimeLag Tol_val Std_val timevary_PSTH smoothVect;  ',...
                            ' global timevary_PSTH smooth_rt TimeLag Tol_val Std_val smoothVect;       ',...
                            ' TimeLag = 300;  ',...
                            ' Tol_val = [...   ',...
                            '     0.100000 0.050000 0.010000 0.005000 0.001000 0.000500 0.000100 0.000050 ];  ',...
                            ' Std_val = [0 0.25 0.5 0.75 1 1.5 2 4 6];  ',...
                            ' timevary_PSTH = 1;   ',...
                            ' smooth_rt = 41;   ',...
                            ' smoothVect = [5 8 29];   ',...
                            ' STRFPAK_core'],data_files{ifiles});
                end % END of fwidthHz_list
            elseif preprocess == 2   % Wavelet Transform
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % You can modify the following parameters:
                %  - wavelength:  It defines the window length of the wavelet
                %  - NBAND:  the number of frequency samp points
                %  - respsamprate: sampling rate of your spike data in Hz
                %  - ampsamprate: the sampling rate of the amplitude evelope
                %                  in Hz.
                %  - filteroption: the choice of scaling amplitude
                %                 1:1 - choose linear
                %                 2:2 - choose log
                %                 1:2 - both scale
                 % ---------------------------------------------------
                % NOTE: You can modify the following calculation parameters:
                % - TimeLag: time lag used for calculating the STA.
                % - Tol_val: a regularization parameter used to estimate
                % STRF. 
                % - Std_val: another regularization parameter in the STRF
                % estimation that enforces sparseness.
                % - timevary_PSTH: flag for whether removing time-varying mean firing
                % rate from inputs.
                % - smooth_rt: window size for smoothing time-varying mean
                % firing rate.
                % - smoothVect: a vector of smoothing window range for
                % validation stage.
                % ---------------------------------------------------
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Call Queue Master to submit your processing request
                Nfilt_list = [21 41 81];
                for ii =1:length(Nfilt_list)
                    dbaddqueuemaster([' cd ',data_files{ifiles},';' ...
                            ' outputPath_temp = ''', output_Path, ''';' ...
                            ' preprocessoption = 2;', ...
                            ' clear global respsamprate wavelength NBAND ampsamprate fwidthHz initialFreq endFreq sDim filteroption; ', ...
                            ' global respsamprate fwidthHz ampsamprate wavelength Npoints initialFreq endFreq sDim filteroption;',...
                            ' fwidthHz = 125; ',...
                            ' Npoints = ',num2str(Nfilt_list(ii)),'; ',...
                            ' respsamprate = 1000; ', ...
                            ' ampsamprate = 1000; ',...
                            ' filteroption = 1:2;  ', ...
                            ' initialFreq = 250;   ',...
                            ' endFreq = 8000;    ',....
                            ' sDim = \''1-D\'';    ',...
                            ' clear global smooth_rt TimeLag Tol_val Std_val timevary_PSTH smoothVect;  ',...
                            ' global timevary_PSTH smooth_rt TimeLag Tol_val Std_val smoothVect;       ',...
                            ' TimeLag = 300;  ',...
                            ' Tol_val = [...   ',...
                            '     0.100000 0.050000 0.010000 0.005000 0.001000 0.000500 0.000100 0.000050 ];  ',...
                            ' Std_val = [0 0.25 0.5 0.75 1 1.5 2 4 6];  ',...
                            ' timevary_PSTH = 1;   ',...
                            ' smooth_rt = 41;   ',...
                            ' smoothVect = [5 8 29];   ',...
                            ' STRFPAK_core'],data_files{ifiles});
                end % end of Nfilt_list
            end % end of preprocess
            
        end   % end of parallel_OK
        
    end % end of total number of datafile  
end % END of preprocess

cd(current_CD);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END of STRFPAK_batch_template
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
