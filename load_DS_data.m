function [DS_data,the_checksum] = load_DS_data(DS,stim_avg,twindow, nband)
%               Pre-loads data to see if its autocorrelation is already
%               cached.  Also computes a checksum of the data that's
%               insensitive to the position of the data; i.e. if you have
%               the same number of trials for two experiments but the
%               stimuli are presented in a different order, we still know
%               the CS is going to be the same.
%               stims  - stimulus file name
%               nlen    - length of time domain
%               ntrials    - num of trials
DS_data = struct([]);
to_hash = struct([]);
for jj = 1:length(DS)

    stim_env = Check_And_Load(DS{jj}.stimfiles);
    to_hash(jj).stims = checksum(stim_env);
    to_hash(jj).nlen = DS{jj}.nlen;
    to_hash(jj).ntrials = DS{jj}.ntrials;
    DS_data(jj).stims = stim_env;
    DS_data(jj).nlen = DS{jj}.nlen;
    DS_data(jj).ntrials = DS{jj}.ntrials;
end

[junk,sorted] = sort({to_hash(:).stims});%  Sort the list according to the stimulus' checksum
the_checksum = checksum([to_hash(sorted)],stim_avg,twindow,nband);%  Generates a checksum which is invariant to shuffling the order of the stims played, so long as ntrials, etc, are the same.
