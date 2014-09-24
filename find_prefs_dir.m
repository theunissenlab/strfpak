function outdir = find_prefs_dir
%  Returns the directory where your preferences are to be stored.  Will use
%  the same directory as dir_of_caches.
outdir = which('dir_of_caches.m');
outdir = fileparts(outdir);
outdir = [outdir filesep];