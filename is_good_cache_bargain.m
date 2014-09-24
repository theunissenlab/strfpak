function is_good_bargain = is_good_cache_bargain(cached_dir,maxsize,the_checksum,time_saved,space_needed)
%  Sees if the cache at hand should replace a less worthy cache file, and does if needed.

if space_needed > maxsize
    is_good_bargain = 0;
else
    active = add_record_to_cache(cached_dir,maxsize,the_checksum,time_saved,space_needed);

    loaded = load(fullfile(cached_dir,'cache_useage_stats.mat'));
    cache_useage_stats = loaded.cache_useage_stats;
    total_space_used = loaded.total_space_used;


    [junk,sorted] = sort([cache_useage_stats(:).bargain]);

    active_in_sorted = find(sorted == active);

    %  Check to see if the sum of the space_neededs of the files which are a worse
    %  bargain than the active one are less than the space_needed of the file of the
    %  active one.
    space_needed_of_worse = sum([cache_useage_stats(sorted(1:(active_in_sorted-1))).space_needed]);
    if space_needed_of_worse < space_needed
        is_good_bargain = 0;
    else
        is_good_bargain = 1;
        space_saved = 0;
        pos_to_delete = 1;
        while space_saved < space_needed
            filename = fullfile(cached_dir,[cache_useage_stats(sorted(pos_to_delete)).checksum '.mat']);
            cache_useage_stats = cache_useage_stats([1:(sorted(pos_to_delete)-1) (sorted(pos_to_delete)+1):end]); 
            if exist(filename,'file')
                delete(filename);
                disp(['deleting file ' filename ]);
                space_saved = space_saved + cache_useage_stats(sorted(pos_to_delete)).space_needed;
            end
            pos_to_delete = pos_to_delete + 1;
        end
        total_space_used = total_space_used - space_saved;
    end
    save(fullfile(cached_dir,'cache_useage_stats.mat'),'cache_useage_stats','total_space_used');

end
