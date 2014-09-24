function out = get_f_ref(filename)
%  Returns the name of the function whose result was cached in the metadata
%  file "filename"

out = '*unknown*';

if exist(filename,'file')
    load(filename);
    if exist('function_name','var')
        out = function_name;
    end
end
