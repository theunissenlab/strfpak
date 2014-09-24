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
global running_in_script_mode
running_in_script_mode = 'yes'; % For those pesky statements where there's soemthing in the code you just have to modify.  In the GUI version, this variable doesn't exist.

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


