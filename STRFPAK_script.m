%  Hello and welcome to STRFPAK_script.  This file is a collection of
%  commands which mirror the instructions performed by the gui version of
%  STRFPAK for the cell in this output directory.
%  This code is dynamically generated as you click through STRFPAK.  
%
%  As of version 4.1, this code looks for two .mat files:
%  STRFPAK_script_parameters.mat and STRFPAK_script_dataset.mat.  The first
%  contains all your options and settings; the second contains the file
%  names of your input datasets and everything specific to that one cell.  
%  We hope it will be easy to "hack" STRFPAK_script_dataset.mat so that 
%  you can do computations for other cells and datasets without having to 
%  run the gui version of STRFPAK each time.  Good luck!
%  Here's an example of the kind of code which makes this script portable:
%  These lines of code assume you have loaded "STRFPAK_script_dataset.mat"
%  and that your original template data set was in
%  /Applications/MATLAB7/work/Theunissen/half_hash/STRFPAK_4.1/DemoData/,
%  and the current data is in pwd.  The current data happens to have the
%  same types of stim and resp names (which makes this easier on the
%  coder), but you get the idea (I hope).
%  for jj = 1:length(rawDS); rawDS{jj}.respfiles = ...
%  strrep(rawDS{jj}.respfiles, '/Applications/MATLAB7/work/Theunissen/half_hash/STRFPAK_4.1/DemoData/',[pwd '/']); end
%  for jj = 1:length(rawDS); rawDS{jj}.stimfiles = ...
%  strrep(rawDS{jj}.stimfiles, '/Applications/MATLAB7/work/Theunissen/half_hash/STRFPAK_4.1/DemoData/',[pwd '/']); end
%  Now make a global variable "outputPath" in pwd, save outputPath and
%  rawDS in the file "STRFPAK_script_dataset.mat" and copy
%  "STRFPAK_script_parameters.mat" and STRFPAK_script.m to the current
%  directory, and run STRFPAK_script.  I know, it's hard the first time,
%  but you can even automate this to suit your database, and soon it will
%  be no problem at all.  Enjoy!
clear global
load STRFPAK_script_parameters.mat %contains all your tol values, smoothing windows, etc.
load STRFPAK_script_dataset.mat %everything specific to the cell in question: filenames of datasets and the output directory.
running_in_script_mode = 'Yes'; % For those pesky statements where there's soemthing in the code you just have to modify.  In the GUI version, this variable doesn't exist.

%%%  NEW in version 4.4:  If STRFPAK typically guesses your input filenames
%%%  corectly, you can try commenting out "load STRFPAK_script_dataset.mat"
%%%  (You still need the line "load STRFPAK_script_parameters.mat")
%%%  and un-commenting the following:

% global rawDS outputPath
% rawDS = guess_rawDS;
% outputPath = fullfile(pwd,'Output');

%%%  End of automatic dataset detection.  Now, the super-lazy way to do
%%%  STRFs is:
%%%      % Run STRFPAK on a demo data set to make sure everything works, 
%%%      % Add the path of the Output directory of this sample cell to your MATLAB path 
%%%      % Make the changes in STRFPAK_script code comments (so that datasets are automatically detected) as above
%%%      % Move to a new data directory
%%%      % type STRFPAK_script


rec_make_dir(outputPath);

%%%^^^ begin preprocess


global respsamprate ampsamprate
if isempty(respsamprate) | isempty(ampsamprate)
    errordlg('Please enter your response data sampling rate.',...
        'Data Input Error', 'modal')
    return
end
global rawDS
if isempty(rawDS)
    errordlg('Please select wave files first.',...
        'Data input error', 'modal')
    return
end

% 2. Ask where you want to put your intermediate results
global outputPath
if length(outputPath) == 0
    initialize_outputPath;
end
if ~exist('running_in_script_mode','var')
    set(handles.outdirshow, 'String', outputPath);
end
% 3. Call ComplexSpetrum to compute specgram for input
%     and assign them to global variable DS
%
numfiles = length(rawDS);
global DS NBAND num_trials
global stimsamprate initialFreq endFreq sDim user_specified_function
sDim='1-D';

%outSpectrum = cell(numfiles, 1);
%stimsamprate = cell(numfiles,1);
currentPath = pwd;

% dB in Noise for the log compression - values below will be set to
% zero.

% do calculation
tempWait = waitbar(0, ['Calculating function ' user_specified_function ', please wait...']);
hashes_of_stims = {};
global the_checksum
the_checksum = '';

for ii=1:numfiles

    waitbar(ii/numfiles, tempWait);

    % 3.1. Take care of Stimulus file
    [input, fs] = wavread(rawDS{ii}.stimfiles);
    [path,name,ext,ver] = fileparts(rawDS{ii}.stimfiles);
    % Spectrogram calculation
    %
    outSpectrum = do_cached_calc(user_specified_function,input);
    if ~isempty(the_checksum)
        hashes_of_stims{ii} = the_checksum; % calculated in do_cached_calc; the_checksum is a global variable
    end

    stimsamprate = fs;

    save(fullfile(outputPath,[name,'_Stim_',num2str(ii),'.mat']), 'outSpectrum');
    if isempty(the_checksum)
        hashes_of_stims{ii} = checksum_from_file(fullfile(outputPath,[name,'_Stim_',num2str(ii),'.mat']));
    end
    % Assign values to global variable DS
    DS{ii}.stimfiles = fullfile(outputPath,[name,'_Stim_',num2str(ii),'.mat']);
    DS{ii}.nlen = round(length(input)*1000/fs) +1;

    %save(fullfile(outputPath,['specgram_Stimlabel_',num2str(ii),'.mat']), 'fo');

    % 3.2. Take care of Response file
    %    If you have multiple trial data, calculate psth first
    %    Then resample it using new amp_samp_rate

    %rawResp = load(rawDS{ii}.respfiles);
    % Modified by Junli, 2003 to read new spike arrivial time file
    %
    [rawResp, trials] = read_spikeTime_2cell(rawDS{ii}.respfiles, round(length(input)*1000/fs) +1);
    [path,name,ext,ver] = fileparts(rawDS{ii}.respfiles);
    %rawResp = read_spikeTime(rawDS{ii}.respfiles, size(outSpectrum, 2));


    % save to the file for each data pair
    save(fullfile(outputPath,[name,'_Spike_time_',num2str(ii),'.mat']), 'rawResp');

    % Assign values to global variable DS
    DS{ii}.respfiles = fullfile(outputPath,[name,'_Spike_time_',num2str(ii),'.mat']);
    DS{ii}.ntrials = trials;

end
%save(fullfile(outputPath,'preprocessed_stim_hashes.mat'),'hashes_of_stims');
put_stim_checksums(hashes_of_stims,DS);

close(tempWait)

NBAND = size(outSpectrum, 1);

global originalDS
originalDS = DS;
%%%^^^ end preprocess


%%%^^^ begin select_validation
%%%^^^ end select_validation

%%%^^^ begin calculate
%%%^^^ end calculate

%%%^^^ begin validation
%%%^^^ end validation

%%%^^^ begin goodness_of_fit
%%%^^^ end goodness_of_fit

%%%^^^ begin prediction
%%%^^^ end prediction


