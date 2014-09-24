%function out = do_cache_cleanup(cached_dir,desired_cache_size);
%  Cleans up the cached dir until its contents are less than desired_cache_size.
%  The formula for the "fitness" of the cache dir is:
%  time_saved/space_needed / (time_offset + min(0,age_of_file in days))
[cached_dir_out,maxsize] = dir_of_caches;
    cached_dir = cached_dir_out;
    desired_cache_size = maxsize * .8;
time_offset = .1;  %  If the file wasn't used in this many days, it starts to get increasingly linear penalty
t_start = now;
working_filename = fullfile(cached_dir,'metadata','deletion_in_progress.mat')
if exist(working_filename,'file')
    disp(['Skipping cache cleanup since another process is doing it.' char(10) ...
        'If this message persists and the cache doesn''t get cleaned, delete the file ' char(10) ...
        working_filename '.']);
    return
else
        working = 'yes';
    save(working_filename,'working');
end
total_space_used = get_total_space_used(cached_dir)

if total_space_used < desired_cache_size
    out = 0;
else
    min_to_trim = total_space_used - desired_cache_size;
    disp(['The cache is getting a little full.  I''m going to trim it to ' num2str((10^(-6))*round((10^6)*desired_cache_size)) ' GB or less.']);
    dirout = dir(fullfile(cached_dir,'metadata','*.mat'));
    expected_vars = {'last_used','space_needed','time_saved'};  %  If these aren't present, the filename isn't bona fide, and we won't touch it.
    goodness_vect = [];
    space_used_vect = [];
    checksum_cells = {};
    for jj = 1:length(dirout)
        the_checksum = dirout(jj).name;
        the_checksum = the_checksum(1:(end-4)); %  chops off'.mat'
        metadata_filename = fullfile(cached_dir,'metadata',dirout(jj).name);
        try
            loaded = load(metadata_filename);
            is_ok = 1;
            for kk = 1:length(expected_vars)
                is_ok = is_ok *isfield(loaded,expected_vars{kk});
            end

            if is_ok
                file_age = now - loaded.last_used;
                file_age = min(0,file_age);
                the_goodness= loaded.time_saved / loaded.space_needed / (time_offset + file_age);
                goodness_vect(end+1) =  the_goodness;
                space_used_vect(end+1) = loaded.space_needed;
                checksum_cells{end+1} = the_checksum;
            end
        catch
            disp(['The function "do_cache_cleanup" encountered a corrupt data file in the cache, so it will purge it.' char(10) ...
                'The file is ' metadata_filename ]);
            try
                delete(metadata_filename);
            catch
                disp(lasterr)

            end
            to_delete_file1 = fullfile(cached_dir,[the_checksum '.mat']);
            try
                delete(to_delete_file1);
            catch
                disp(lasterr)

            end
        end
    end

    [junk,good_ind] =sort(goodness_vect);
    size_deleted = 0;
    active_file = 1;
    switch filesep
        case '/'
                while (active_file < (length(good_ind) + .5)) & (size_deleted < min_to_trim)
        the_checksum = checksum_cells{good_ind(active_file)};
        to_delete_file1 = fullfile(cached_dir,[the_checksum '.mat']);
        to_delete_file2 = fullfile(cached_dir,'metadata',[the_checksum '.mat']);
        f_deleted_name = get_f_ref(to_delete_file2);
        try
            system(['rm -f ' to_delete_file1]);
        catch
            disp(lasterr)
        end
        try
            system(['rm -f ' to_delete_file2]);
        catch
            disp(lasterr)
        end
        size_deleted = size_deleted + space_used_vect(good_ind(active_file));
        disp(['Deleting record of ' f_deleted_name ', which had a goodness of ' num2str(goodness_vect(good_ind(active_file))) '.']);
        active_file = active_file + 1;
    end

        case '\'
    while (active_file < (length(good_ind) + .5)) & (size_deleted < min_to_trim)
        the_checksum = checksum_cells{good_ind(active_file)};
        to_delete_file1 = fullfile(cached_dir,[the_checksum '.mat']);
        to_delete_file2 = fullfile(cached_dir,'metadata',[the_checksum '.mat']);
        f_deleted_name = get_f_ref(to_delete_file2);
        try
            delete(to_delete_file1);
        catch
            disp(lasterr)
        end
        try
            delete(to_delete_file2);
        catch
            disp(lasterr)
        end
        size_deleted = size_deleted + space_used_vect(good_ind(active_file));
        disp(['Deleting record of ' f_deleted_name ', which had a goodness of ' num2str(goodness_vect(good_ind(active_file))) '.']);
        active_file = active_file + 1;
    end
    end
    delete(working_filename);
end
