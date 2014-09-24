function preload_stims = make_ready_stims(predDS,stim_avg)
preload_stims = {};
for jj = 1:length(predDS)
    filename = predDS{jj}.stimfiles;
    [fndir, fname, fext] = fileparts(filename);
    filename = [fname fext];
    global outputPath
    outfilename = fullfile(outputPath,[filename(1:(end-4)) '_mean_removed.mat']);
    if ~exist(outfilename,'file')
        stim_env = Check_And_Load(predDS{jj}.stimfiles);
        nband = size(stim_env,1);
        nlen = size(stim_env,2);
        stimval = zeros(nband, nlen);
        for tgood = 1:nlen
            stimval(:, tgood) = stim_env(1:nband, tgood) - stim_avg(1:nband);
        end
        save(outfilename,'stimval');
    else
        if nargout > 0
            load(outfilename);
        end
    end
    if nargout > 0
        preload_stims{jj} = stimval;
    end
end
