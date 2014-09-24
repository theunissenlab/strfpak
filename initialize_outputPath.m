if ~exist('outputPath','var')
    global outputPath
    if isempty(outputPath)
        currentPath = pwd;
        prompt={['Please Enter where you want to put your intermediate results:']};
        def = {fullfile(currentPath, 'Output')};
        dlgTitle='Path for intermediate results';
        lineNo=1;

        % picture feature
        AddOpts.Resize='on';
        AddOpts.WindowStyle='normal';
        AddOpts.Interpreter='tex';
        datadir =inputdlg2(prompt,dlgTitle,lineNo,def,AddOpts);

        % Check if user input valid directory
        if isempty(datadir)
            errordlg('Please enter valid directory.','Input Error','modal')
            return
        end

        outdatafiledir = datadir{1};

        if not(exist(outdatafiledir,'dir'))
            disp('Directory not found. Creating new directory.');
            [p, n, e] = fileparts(outdatafiledir);
            if not(exist(p, 'dir'))
                errordlg('The upper directory is not found. So exist.');
                return
            end
            cd (p)
            mkdir(n)
            % Junli: 11/11/2005
        else   % existing directory
            tt = dir(outdatafiledir);
            if length(tt) >2    % more file there

                anw = questdlg({'The directory is not empty. What do you want to do?'},...
                    'Warning Message', 'Overwrite the directory', 'Create new directory',...
                    'Create new directory');
                switch anw
                    case 'Overwrite the directory'
                        current_dir = pwd;
                        cd(outdatafiledir)
                        delete *_*.mat;
                        cd(current_dir);

                    case 'Create new directory'

                        newDir = 0;
                        while newDir == 0
                            currentPath = pwd;
                            prompt={['Please Enter where you want to put your intermediate results:']};
                            def = {fullfile(currentPath, 'Output')};
                            dlgTitle='Path for intermediate results';
                            lineNo=1;

                            % picture feature
                            AddOpts.Resize='on';
                            AddOpts.WindowStyle='normal';
                            AddOpts.Interpreter='tex';
                            datadir =inputdlg2(prompt,dlgTitle,lineNo,def,AddOpts);

                            % Check if user input valid directory
                            if isempty(datadir)
                                errordlg('Please enter valid directory.','Input Error','modal')
                                return
                            end

                            outdatafiledir = datadir{1};

                            if not(exist(outdatafiledir,'dir'))
                                disp('Directory not found. Creating new directory.');
                                [p, n, e] = fileparts(outdatafiledir);
                                if not(exist(p, 'dir'))
                                    errordlg('The upper directory is not found. Exist.');
                                    return
                                end
                                cd (p)
                                mkdir(n)
                                newDir == 1;
                            else
                                tt =warndlg('Directory is exsiting. Please create new one.', 'Warning', 'modal');
                                uiwait(tt);
                            end
                        end
                    otherwise
                        errordlg('Computation is cancelled.', 'modal')
                        return;
                end

            end
        end

        outputPath = outdatafiledir;
    else   % not empty outputPath
        outdatafiledir = outputPath;
        tt = dir(outdatafiledir);
        if length(tt) >2    % more file there

            anw = questdlg({'The directory is not empty. What do you want to do?'},...
                'Warning Message', 'Overwrite the directory', 'Create new directory',...
                'Create new directory');
            switch anw
                case 'Overwrite the directory'
                    current_dir = pwd;
                    cd(outdatafiledir)
                    delete *_*.mat;
                    cd(current_dir);

                case 'Create new directory'
                    newDir = 0;
                    while newDir == 0
                        currentPath = pwd;
                        prompt={['Please Enter where you want to put your intermediate results:']};
                        def = {fullfile(currentPath, 'Output')};
                        dlgTitle='Path for intermediate results';
                        lineNo=1;

                        % picture feature
                        AddOpts.Resize='on';
                        AddOpts.WindowStyle='normal';
                        AddOpts.Interpreter='tex';
                        datadir =inputdlg2(prompt,dlgTitle,lineNo,def,AddOpts);

                        % Check if user input valid directory
                        if isempty(datadir)
                            errordlg('Please enter valid directory.','Input Error','modal')
                            return
                        end

                        outdatafiledir = datadir{1};

                        if not(exist(outdatafiledir,'dir'))
                            disp('Directory not found. Creating new directory.');
                            [p, n, e] = fileparts(outdatafiledir);
                            if not(exist(p, 'dir'))
                                errordlg('The upper directory is not found. Exist...');
                                return
                            end
                            cd (p)
                            mkdir(n)
                            newDir = 1;
                        else
                            tt =warndlg('Directory exists. Please create new one.', 'Warning', 'modal');
                            uiwait(tt);
                        end
                    end
                otherwise
                    errordlg('Computation is cancelled.', 'modal')
                    return;
            end
        end
        outputPath = outdatafiledir;
    end
end
