function hashes_of_pred_stims = create_hash_cache_file_PDS(outputPath,PDS);
%  Creates the file with the hashes of the stimuli
hashes_of_pred_stims = {};
global DS
[cached_dir,maxsize] = dir_of_caches;
hash_cache_filename = fullfile(outputPath,'preprocessed_stim_hashes.mat');
if ~exist(hash_cache_filename,'file')
    hashes_of_stims = create_hash_cache_file(outputPath,DS);
else
    load(hash_cache_filename); % puts "hashes_of_stims" into the workspace
end

for ii = 1:length(PDS)
    if (length(PDS) >= ii ) & (length(DS) >= ii)
        if strcmp(PDS{ii}.stimfiles,DS{ii}.stimfiles)
            hashes_of_pred_stims{ii} = hashes_of_stims{ii};
        else
            hashes_of_pred_stims{ii} = checksum_from_file(PDS{ii}.stimfiles);
        end
    else
        hashes_of_pred_stims{ii} = checksum_from_file(PDS{ii}.stimfiles);
    end
end
save(fullfile(outputPath,'preprocessed_pred_stim_hashes.mat'),'hashes_of_pred_stims');
