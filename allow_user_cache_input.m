function using_cache = allow_user_cache_input

[cached_dir,maxsize] = dir_of_caches;

if strcmp(cached_dir,'null')
    using_cache = 0;
    disp(['First type "help Caching_Readme", then type "edit dir_of_caches.m" to get started with caching.'])
else
    using_cache=1;
end
