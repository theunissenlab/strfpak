moddir = '/auto/fhome/pgill/STRFPAK_bleeding_edge';

dirout = dir(fullfile(moddir,'*.m'));
findtext = 'newSlider = get(hObject, ''Value'');';
replacetext = 'newSlider = round(get(hObject, ''Value''));';
for jj = 1:length(dirout)
   fname = fullfile(moddir,dirout(jj).name);
    fid = fopen(fname,'r');
    textdump = char(fread(fid,'char')');
    fclose(fid);
    fid = fopen(fname,'w');
    textdump = strrep(textdump,findtext,replacetext);
    fwrite(fid,textdump,'char');
    fclose(fid);
end
