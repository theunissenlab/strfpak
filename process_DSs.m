function [DS_names, predDS_names,finalDS_names,for_validation,for_final,DS,predDS,finalDS] = process_DSs(DS,predDS,finalDS,rawDS);
%  Function used by select_validation_files.m
%  Takes the structs DS, predDS and rawDS and returns "nice" (i.e.
%  user-recognizable) strings of the current data sets (both DS and predDS)
%  as well as a vector of 1s and 0s saying if the corresponding entry in
%  rawDS has been reserved for validation.
for_validation = -999*ones(1,length(rawDS));
for_final = -999*ones(1,length(rawDS));

%  First, let's get the names of the files we're working with.
for jj = 1:length(rawDS)
    this_raw_pair = rawDS{jj};
    this_stimfile = this_raw_pair.stimfiles;
    this_respfile = this_raw_pair.respfiles;
    [stimdir,stimname] = fileparts(this_stimfile);
    stimnames{jj} = stimname;
    [respdir,respname] = fileparts(this_respfile);
    respnames{jj} = respname;
end

%  Next, let's get the names of the files in DS
global rawData
DS_stimnames = {};
DS_respnames = {};

for jj = 1:length(DS)
    this_raw_pair = DS{jj};
    this_stimfile = this_raw_pair.stimfiles;
    this_respfile = this_raw_pair.respfiles;
    [stimdir,stimname] = fileparts(this_stimfile);
    [respdir,respname] = fileparts(this_respfile);
    if rawData
        underscores = findstr(stimname,'_');
        stimname = stimname(1:(underscores(end-1)-1));
    end
    DS_stimnames{jj} = stimname;
    if rawData
        underscores = findstr(respname,'_');
        respname = respname(1:(underscores(end-2)-1));
    end
    DS_respnames{jj} = respname;
end

%  Next, let's get the names of the files in predDS
predDS_stimnames = {};
predDS_respnames = {};
for jj = 1:length(predDS)
    this_raw_pair = predDS{jj};
    this_stimfile = this_raw_pair.stimfiles;
    this_respfile = this_raw_pair.respfiles;
    [stimdir,stimname] = fileparts(this_stimfile);
    [respdir,respname] = fileparts(this_respfile);
    if rawData
        underscores = findstr(stimname,'_');
        stimname = stimname(1:(underscores(end-1)-1));
    end
    predDS_stimnames{jj} = stimname;
    if rawData
        underscores = findstr(respname,'_');
        respname = respname(1:(underscores(end-2)-1));
    end
    predDS_respnames{jj} = respname;
end


%  Finally, let's get the names of the files in finalDS
finalDS_stimnames = {};
finalDS_respnames = {};
for jj = 1:length(finalDS)
    this_raw_pair = finalDS{jj};
    this_stimfile = this_raw_pair.stimfiles;
    this_respfile = this_raw_pair.respfiles;
    [stimdir,stimname] = fileparts(this_stimfile);
    [respdir,respname] = fileparts(this_respfile);
    if rawData
        underscores = findstr(stimname,'_');
        stimname = stimname(1:(underscores(end-1)-1));
    end
    finalDS_stimnames{jj} = stimname;
    if rawData
        underscores = findstr(respname,'_');
        respname = respname(1:(underscores(end-2)-1));
    end
    finalDS_respnames{jj} = respname;

end


DS_names = {};
predDS_names = {};
rawDS_names = {};
finalDS_names = {};
for jj = 1:length(for_validation)
    done = 0;
    this_name = [stimnames{jj} '  ' respnames{jj}];
    rawDS_names{end+1}  =this_name;
    for kk = 1:length(DS)
        if 1
            this_DS_name = [DS_stimnames{kk} '  ' DS_respnames{kk}];
            if strcmp(this_name,this_DS_name)

                for_validation(jj) = 0;
                for_final(jj) = 0;
                DS_names{kk} = this_name;
            end
        end
    end

    for kk = 1:length(predDS)
        if 1
            this_predDS_name = [predDS_stimnames{kk} '  ' predDS_respnames{kk}];
            if strcmp(this_name,this_predDS_name)
                for_validation(jj) = 1;
                for_final(jj) = 0;
                predDS_names{kk} = this_name;
            end
        end
    end
    for kk = 1:length(finalDS)
        if 1
            this_finalDS_name = [finalDS_stimnames{kk} '  ' finalDS_respnames{kk}];
            if strcmp(this_name,this_finalDS_name)
                for_validation(jj) = 0;
                for_final(jj) = 1;
                finalDS_names{kk} = this_name;
            end
        end
    end
end
