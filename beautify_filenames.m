function DS_names = beautify_filenames(DS);
global rawData
DS_stimnames = {};
DS_respnames = {};
DS_names = {};

for jj = 1:length(DS)
    this_raw_pair = DS{jj};
    this_stimfile = this_raw_pair.stimfiles;
    this_respfile = this_raw_pair.respfiles;
    [stimdir,stimname] = fileparts(this_stimfile);
    if rawData
        underscores = findstr(stimname,'_');
        if isempty(underscores)
            underscores = [(length(stimname)+1) 0];
        end
        stimname = stimname(1:(underscores(end-1)-1));
    end
    DS_stimnames{jj} = stimname;

    this_raw_pair = DS{jj};
    this_respfile = this_raw_pair.respfiles;
    this_respfile = this_raw_pair.respfiles;
    [respdir,respname] = fileparts(this_respfile);
    if rawData
        underscores = findstr(respname,'_');
        if isempty(underscores)
            underscores = [(length(respname)+1) 0 0];
        end

        respname = respname(1:(underscores(end-2)-1));
    end
    DS_respnames{jj} = respname;

    this_name = [DS_stimnames{jj} '  ' DS_respnames{jj}];
    DS_names{end+1}  =this_name;
end
