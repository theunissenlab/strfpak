function [strf_checksum,the_dir] = get_strf_checksum(strfFiles);
posslash = findstr(strfFiles,filesep);
the_dir = strfFiles(1:(posslash(end)));
the_name = strfFiles((posslash(end)+1):end);
strf_hash_filename = [the_dir 'hash_of_' the_name];
if exist(strf_hash_filename,'file')

    load(strf_hash_filename);

else
    %disp(['Problem with "get_strf_checksum": I thought the strf checksum was already calculated.' char(10) ...
    %    'Doing it the slow way ...']);  % The strf checksum can take a long time to calculate.  I hope that a faster checksum
    %   will have been already made in the calStrf_script, but if not, we do it here.
    strf_checksum = checksum_from_file(strfFiles);
    save(strf_hash_filename,'strf_checksum');
end
