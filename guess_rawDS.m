function guessDS = guess_rawDS(lookdirstim,lookdirresp)
if ~exist('lookdirstim','var')
    lookdirstim = pwd;
end
if ~exist('lookdirresp','var')
    lookdirresp = pwd;
end
guessDS = {};
diroutstim = dir(lookdirstim);
diroutresp = dir(lookdirresp);
possible_stims = zeros(1,length(diroutstim));
possible_resps = zeros(1,length(diroutresp));
for jj = 1:length(diroutstim)
    if ~diroutstim(jj).isdir
        fname = diroutstim(jj).name;
        if strcmp(fname((end-3):end),'.wav')  | strcmp(fname(1:4),'stim') |strcmp(fname(1:4),'Stim') %Here we see if this file looks like a stim file
            possible_stims(jj) = 1;
        end
    end
end
for jj = 1:length(diroutresp)
    if ~diroutresp(jj).isdir
        fname = diroutresp(jj).name;
        if ~strcmp(fname((end-3):end),'.wav')&~possible_stims(jj)  %Here we see if this file looks like a resp file
            if ~strcmp(fname(end),'~')
            possible_resps(jj) = 1;
            end
        end
    end
end

possible_stim_names = {diroutstim(find(possible_stims)).name};
possible_resp_names = {diroutresp(find(possible_resps)).name};

possible_stim_names = sort(possible_stim_names);
possible_resp_names = sort(possible_resp_names);
resp_not_used = ones(1,length(possible_resp_names));
stim_not_used = ones(1,length(possible_stim_names));

sort_by_number = 1;
guess_value = [];
for jj = 1:length(possible_stim_names)
    putstim = possible_stim_names{jj};
    for kk = find(resp_not_used)
        if stim_not_used(jj)
            putresp = possible_resp_names{kk};
            if strcmp(putstim(1:(end-4)),putresp(1:(end-4))) |  strcmp(putstim(1:(end-4)),putresp) %  If you have xtrfd.wav and xtrfd.dat
                guessDS{end+1}.stimfiles = fullfile(lookdirstim,putstim);
                guessDS{end}.respfiles = fullfile(lookdirresp,putresp);
                resp_not_used(kk) = 0;
                stim_not_used(jj) = 0;
                sort_by_number = 0;
            else
                stimdigits = putstim(find((putstim < 57.5) .* (putstim > 47.5)));  %  This will extract all of the digits in the stim filename, so 'stim12.wav' will become '12'
                respdigits = putresp(find((putresp < 57.5) .* (putresp > 47.5)));
                if strcmp(stimdigits,respdigits)
                    guessDS{end+1}.stimfiles = fullfile(lookdirstim,putstim);
                    guessDS{end}.respfiles = fullfile(lookdirresp,putresp);
                    resp_not_used(kk) = 0;
                    stim_not_used(jj) = 0;
                    guess_value(end+1) = str2double(stimdigits);
                end
            end
        end
    end
end

if sort_by_number
    [junk,order] = sort(guess_value);
    guessDS = {guessDS{order}};
end
if length(guessDS) > 0
    disp(['STRFPAK made a guess that your stim files look like "' putstim '" and your data files look like "' putresp '".']);
else
    disp('STRFPAK could not recognize any likely stim/response pairs.');
end
