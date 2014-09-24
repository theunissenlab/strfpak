function out = put_stim_checksums(hashes_of_stims,DS);
global outputPath
rec_make_dir(fullfile(outputPath,'stim_hashes'));
for ii = 1:length(DS)
    dsname = DS{ii}.stimfiles;
    [fdir,fname,fext] = fileparts(dsname);
    dsname = [fname fext];
    filename = fullfile(outputPath,'stim_hashes',dsname);
    if ~exist(filename,'file')
        this_hash = hashes_of_stims{ii};
        save(filename,'this_hash');
    end
end
