function total_space_used = update_total_space_used(cached_dir);
loaded = load(fullfile(cached_dir,'cache_useage_stats.mat'));
cache_useage_stats = loaded.cache_useage_stats;
total_space_used = 0;
for jj = 1:length(cache_useage_stats)
    if exist(fullfile(cached_dir,[cache_useage_stats(jj).checksum '.mat']),'file')
        total_space_used = total_space_used + cache_useage_stats(jj).space_needed;
    end
end

save(fullfile(cached_dir,'cache_useage_stats.mat'),'cache_useage_stats','total_space_used');
