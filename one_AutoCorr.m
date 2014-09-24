function [CS,CS_ns] = one_AutoCorr(stim_env,nband,nlen,stim_avg,ntrials)
running_flag = 1;

    nlen = size(stim_env,2);
    xb = 1;
    stimval = zeros(nband, nlen);
    
%     % Check if input data are chosen properly.
%     thisLength = size(stim_env, 2);
%     
%     if thisLength < nlen
%         answ= questdlg(['Data Error: Please check your input data by clicking ',...
%                 '"Get Files" Button in the main window: The first data file need ',...
%                 'to be stimuli and the second data file need to be its corresponding',...
%                 ' response file. If you made a mistake, please type "clear all" ',...
%                 ' or hit "reset" button first and then choose input data again.',...
%                 ' Otherwise, we will truncate the longer length. Do you want to continue?'],...
%             'Input Data Error','Yes','No','No');
%         switch answ
%             case 'No'
%                 errFlg = 1;
%                 return
%         end
%     end
    
%     nlen = min(nlen, thisLength);

    % subtract mean of stim from each stim 
    for tgood = 1:nlen
        stimval(:, tgood) = stim_env(1:nband, tgood) - stim_avg(1:nband);
    end
    
    % do autocorrelation calculation
    if running_flag == 1
        tempWait = waitbar(0,...
            sprintf('Calculating Auto-correlation for dataset'));
    end
    for ib1 = 1:nband
        
        if running_flag == 1
            waitbar(ib1/nband, tempWait);
        end
        
        for ib2 = ib1:nband
            
            % NEW version of algorithm by using xcorr 
            CS(xb, :) = CS(xb, :) + xcorr(stimval(ib2,:), stimval(ib1, :),...
                twindow(2))*ntrials;
            
            xb = xb +1;
        end            % END of ib2
    end                 % END of ib1
    if running_flag == 1
        close(tempWait)
    end 
    
    % Count the total trials for later normalization 
    lengthVec = ones(1, nlen);
    CS_ns(1, :) = CS_ns(1, :) + DS{fidx}.ntrials *...
        xcorr(lengthVec, lengthVec, twindow(2));
