function out = add_to_STRFPAK_script(source_code_filename,subroutine_tag,script_filename);
%  Copies code out of STRFPAK into a script file.
%  Specifically, loads the source code file, looks for %%%$$$ begin subroutine_tag,
%  then splices %%%^^^ subroutine_tag and the source code contents into
%  STRFPAK_script.m.

script_tag = '%%%^^^';  %These are the section delimiters you'll find in STRFPAK_script.m.
source_tag = '%%%$$$';  %These are the delimiters in the source code.
if ~exist('script_filename','var')
    global outputPath
    initialize_outputPath;
    script_filename = fullfile(outputPath,'STRFPAK_script.m');
    if ~exist(script_filename,'file')
        script_text = initialize_STRFPAK_script(script_filename);
    else
        fid = fopen(script_filename,'r');
        script_text = char(fread(fid,'char')');
        fclose(fid);
    end
else
    if ~exist(script_filename,'file')
        script_text = initialize_STRFPAK_script(script_filename);
    else
        fid = fopen(script_filename,'r');
        script_text = char(fread(fid,'char')');
        fclose(fid);
    end
end

insert_position = length(script_text) + 1;  % If the section doesn't yet exist, add it to the end.

%  First, let's see if the section already exists in our intrepid script
%  file.

starting_script_tag = [script_tag ' begin ' subroutine_tag];
ending_script_tag = [script_tag ' end ' subroutine_tag];

script_start = findstr(starting_script_tag,script_text);
if length(script_start) > 0
    script_end = findstr(ending_script_tag,script_text);
    if length(script_end) == 0
        disp(['There has been a serious error in the generation of STRFPAK_script.m:' char(10) ...
            'I found the tag "' starting_script_tag '" but not the corresponding end tag "' ending_script_tag '".' char(10) ...
            'STRFPAK_script probably won''t work now.']);
        insert_position = script_start;
        end_position = length(script_text) + 1;

    else
        if (length(script_start) > 1 ) | (length(script_end) > 1)
            disp(['I found duplicate instances of the tag "' subroutine_tag '" in STRFPAK_script.' char(10) ...
                'I''m modifying the first one, and there''s a chance STRFPAK_script will still work, but I''m not promising anything...']);
            script_start = script_start(1);
            script_end = script_end(1);
        end
        insert_position = script_start + length(starting_script_tag);
        end_position = script_end ;
    end
else
    end_position = length(script_text) + 1;
end

%  Now it's time to read the source code file.

fid = fopen(which(source_code_filename),'r');
source_text = char(fread(fid,'char')');
fclose(fid);

%  Now we see if the source code tag exists where it should.


starting_source_tag = [source_tag ' begin ' subroutine_tag];
ending_source_tag = [source_tag ' end ' subroutine_tag];


source_start = findstr(starting_source_tag,source_text);
if length(source_start) > 0
    source_end = findstr(ending_source_tag,source_text);
    if length(source_end) == 0
        disp(['There has been a serious error in the generation of STRFPAK_script.m:' char(10) ...
            'I found the tag "' starting_source_tag '"  in "' source_code_filename '"  but not the corresponding end tag.' char(10) ...
            'STRFPAK_script probably won''t work now.']);
        copy_start_position = source_start;
        copy_end_position = length(source_text);

    else
        if (length(source_start) > 1 ) | (length(source_end) > 1)
            disp(['I found duplicate instances of the tag "' subroutine_tag '" in "' source_code_filename '".' char(10) ...
                'I''m modifying the first one, and there''s a chance STRFPAK_script will still work, but I''m not promising anything...']);
            source_start = source_start(1);
            source_end = source_end(1);
        end
        copy_start_position = source_start + length(starting_source_tag);
        copy_end_position = source_end - 1;% + length(ending_source_tag);
    end
else
    disp(['Error: couldn''t find the tag "' subroutine_tag '" in file "' source_code_filename '".'  char(10) ...
        'Skipping...']);
    copy_start_position = length(source_text);
    copy_end_position = length(source_text);
end

copied = source_text(copy_start_position:copy_end_position);
%  Now all we have to do is write to STRFPAK_script.m

script_text = [script_text(1:(insert_position)) copied script_text((end_position):end)];

fid = fopen(script_filename,'w');
fwrite(fid,script_text,'char');
fclose(fid);
