function [out,order] = load_and_sort_stim_caches(hash_cache_filename);
load(hash_cache_filename);  %  Puts "hashes_of_stims" into workspace
global for_validation
if isempty(for_validation)
[out,order] = sort(hashes_of_stims);
else
    hashes_of_stims = {hashes_of_stims{find(for_validation == 0)}};
    [out,order] = sort(hashes_of_stims);
end
