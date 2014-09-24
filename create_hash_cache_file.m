function hashes_of_stims = create_hash_cache_file(outputPath,DS);
%  Creates the file with the hashes of the stimuli
hashes_of_stims = {};
rec_make_dir(fullfile(outputPath,'stim_hashes'))
for ii = 1:length(DS)
    filename = fullfile(outputPath,'stim_hashes',DS{ii}.stimfiles);
    if ~exist(filename,'file')
        this_hash = checksum_from_file(DS{ii}.stimfiles);
        save(filename,'this_hash');
    else
        load(filename);
    end
    hashes_of_stims{ii} = this_hash;
end
