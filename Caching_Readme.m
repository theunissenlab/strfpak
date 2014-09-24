%%%%%%   Welcome to the magical world of caching!   %%%%%%
%
%
%%    Here's an introduction to what you'll need to know.
%
%  Outline: 1)  Why cache?
%           2)  Caching automatically
%           3)  What caching will do
%           4)  How to set up or stop caching
%           5)  Details
%
%
%%%%%  1)  Why cache?
%
%  You might run the same stimulus sets over and over again in your
%  work.  Wouldn't it be great if you didn't have to re-compute them every
%  time you run a cell through STRFPAK?  
%  With caching, you can save the output of long computations in a central
%  directory, so every time you need them you just have to look up the
%  answers, not re-compute them.
%
%%%%%  2)  Caching automatically
%
%  One of the big worries with caching is that you don't want to have to
%  worry about making sure you're looking up the appropriate answer.
%  Relax: this system computes a checksum for every calculation it does, so
%  you know you'll never look up a wrong result.
%
%%%%%%%%%%%%%  2.1)  What's a checksum?
%
%  A checksum is a 32 character string of random-looking hex numbers which all
%  depend sensitively on every aspect of its input, which can be arbitrary.  For example: 
%
%  checksum('All your base are belong to us.') = B73B4F4BD07013881F92317D3B4C842F
%
%  but
%
%  checksum('All your Base are belong to us.') = 17FBEF6B705093C8FF32F13D3BCC040F
%
%  Notice that even though the input strings are almost identical, every
%  second character between the strings is completely different.  Note
%  also that the input can be arbitrarily long, but the checksum will
%  always be 32 characters long.  Also note that the "every second
%  character" limitation is unique to strings; e.g. try 
%  A = checksum('All your Base are belong to us.',1)
%  B = checksum('All your Base are belong to us.',2)
%  A == B
%  You will see that only two characters in the hashes are the same.
%
%%%%%  3)  What caching will do
%
%  If caching is enabled, any time STRFPAK finds a calculation it
%  thinks it might repeat, it will be saved to the directory specified by
%  "dir_of_caches" under the filename which is the checksum of the function
%  and all the data used as inputs to the function.  If this checksum
%  exists, there's no need to re-compute the calculation.  The file is
%  loaded without repeating the calculation.
%
%%%%%  4)  How to set up or stop caching
%
%  Step 1: Open the file "dir_of_caches.m"  (type "edit dir_of_caches.m"
%  from the directory where STRFPAK is).  This is the preferences file for
%  caching; it tells STRFPAK everything it needs to know about where to put
%  cached results and how much space to use up in the storage of cached
%  results.
%  
%  Step 2: Put the full path (that's "C: etc" for Windows and "/ etc" for
%  Linux/Mac OS X) of the directory where you want to store intermediate
%  results into the line "cache_dir = '".  Note: single quotes (') are needed
%  around the path name.  For best results, choose a directory which every
%  computer working with your data set can see and write to.
%
%         2.1: If you want to use the same code as other STRFPAK users but
%         you want to have a separate directory to save intermediate
%         results, copy "dir_of_caches.m" into a directory you can write
%         to, modify it as per above, and add it to your path AFTER you add
%         STRFPAK 4.x to your path.  Then, MATLAB will look to your version
%         of dir_of_caches before it looks at the version in the shared
%         STRFPAK 4.x folder.  (To see which version of dir_of_caches.m is
%         being used, type "which dir_of_caches.  You'll also want to start
%         STRFPAK from a directory other than the STRFPAK 4.x directory.)
%
%  Step 3: Set the target maximum size (in GB) of the directory where you want to 
%  store intermediate results.  The more generous you are the faster
%  STRFPAK will be.  (I find 10 GB to be plenty, but this really depends on
%  the data you work with.)
%
%  Step 4: If ever you want to stop caching, edit dir_of_caches.m again so
%  that the line "cache_dir = '" becomes "cache_dir = 'null';". 
%
%%%%%  5)  Details (less important than #4, above - start there)
%
%  Keeping your trash/recycle bin/etc. free of temporary files
%
%  If the cache ever gets too full (i.e., exceeds a user-specified disk space limit),
%  STRFPAK will delete a few of the least useful cache files with the
%  "delete" command.  It also uses the "delete" command every time it evaluates a 
%  checksum.  You do not want these files clogging up your recycle
%  bin, so either change your MATlab preferences so that deleted files are
%  removed (not sent to a to-be-deleted directory), or, if you're a MATLAB
%  wizard, change the delete commands in: "concat_for_checksum.m" and
%  "do_cache_cleanup.m" to commands using "!rm -f", or whatever is
%  appropriate for your platform.
%
%  Reducing network traffic to a minimum (for computer clusters)
%  
%  Any time a checksum is computed, a temporary file exists for a split
%  second.  That temporary file should be written to a local hard drive to
%  keep network traffic down.  We always recommend copying data to a local
%  hard drive for speed reasons, but if you want to do calculations
%  directly on a remote fileserver anyway, we strongly advise you do not
%  use the remote system for these temporary files.  To ensure these files
%  get written locally, set the MATLAB environment 'TMP' to a local hard
%  drive.
