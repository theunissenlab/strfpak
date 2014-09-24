function out =md5checksum(varargin)
 testnum = 0;
done = 0;
while done == 0
    tempfile = fullfile(getenv('TMP'),'/tmp',['Temp_hashing_name_' num2str(testnum) '.tmp']);
    if ~exist(tempfile,'file')
        %save(tempfile,'varargin');
        A = ver('MATLAB');
        A = A.Version;
        A = str2num(A(1));
        if A > 6
            to_eval = ['save ' tempfile ' varargin -V6'];
        else
            to_eval = ['save ' tempfile ' varargin'];
        end
        eval(to_eval);
       [h,r] = system(['ls -l ' tempfile]);
       spaces = strfind(r,' ');
       spaces = spaces(find(spaces(2:end) > (1+spaces(1:(end-1)))));
       flen = str2num(r((spaces(4)+1):(spaces(5)-1)));
       [junk,out] = system(['tail -c ' num2str(flen - 75) ' ' tempfile ' | /sw/bin/md5sum']);
       out = out(1:32);
       done = 1;
    else
        testnum = testnum + 1;
    end
end
delete(tempfile);
