function script_text = initialize_STRFPAK_script(script_filename);
if ~exist('script_filename','var')
    global outputPath
    initialize_outputPath;
    script_filename = fullfile(outputPath,'STRFPAK_script.m');
end

fid = fopen(which('SP_script_header.m'),'r');
script_text = char(fread(fid,'char')');
fclose(fid);
fid = fopen(script_filename,'w');
fwrite(fid,script_text,'char');
fclose(fid);
