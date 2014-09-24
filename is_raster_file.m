function [out,loaded,loaded2] = is_raster_file(filename)
%Returns 1 is the file is likely to be a raster file, 0 otherwise.
%Assumption: no decimals are used in raster files, and decimals are used in
%every spike time file.  Also, not all trials have monotonically-increasing
%firing rates.
[path, name, ext] = fileparts(filename);
switch ext
    case {'.mat'}
        loaded = load(filename);
        flds = fieldnames(loaded);
        rawResp = getfield(loaded, char(flds{1}));
        if iscell(rawResp)
            out = 0;
        else
            out = prod(size(rawResp)) > 1;
        end
    otherwise
        fid = fopen(filename);
        loaded = char(fread(fid,'char')');
        fclose(fid);
        out = 1;
        if any(loaded == 46)
            out = 0;
            %disp('detected a decimal');
        end
        if out
            try
                loaded2 = textread(filename,'','emptyvalue',-999);
            catch
                out = 0;
            end
            if out
                if(any(any(loaded2 == -999)))
                    out = 0;
                else
                    all_increasing = 1;
                    for jj=1:size(loaded2,1)
                        if all_increasing
                            nonzero_parts = loaded2(jj,find(loaded2(jj,:)>0));
                            if length(nonzero_parts) > 1
                                diffs = nonzero_parts(2:end) - nonzero_parts(1:(end-1));
                                if any(diffs < 0)
                                    all_increasing = 0;
                                end
                            end
                        end
                    end
                    if all_increasing
                        out = 0;
                        %disp('All rows of the input file had monotonically-increasing entries, so it's probably a spike-time file.' );
                    end
                end
            end
        end
end