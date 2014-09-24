function sorted = sortrawDS(rawDS)
rawDS_names = beautify_filenames(rawDS);
for jj = 1:length(rawDS_names)
    this_name = rawDS_names{jj};
    spaces = findstr(this_name,' ');
    rawDS_names{jj} = [this_name(1:((spaces(1)-1))) ' !'];
end
rawDS_names
[junk,order] = sort(rawDS_names);
order
sorted = {rawDS{order}};
