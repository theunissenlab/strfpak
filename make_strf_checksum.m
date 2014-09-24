function [strf_checksum,the_dir] = make_strf_checksum(strfFiles);
posslash = findstr(strfFiles,filesep);
the_dir = strfFiles(1:(posslash(end)));
the_name = strfFiles((posslash(end)+1):end);
strf_hash_filename = [the_dir 'hash_of_' the_name];
if exist(strf_hash_filename,'file')
    load(strf_hash_filename);
else
    strf_checksum = checksum_from_file(strfFiles);
    save(strf_hash_filename,'strf_checksum');
end
