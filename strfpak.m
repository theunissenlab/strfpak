%
% This matlab script file is first program of STRFPAK. 
% It does:
%      1. First Add the current directory to matlab path
%      2. Run firstGUI of STRF
%
%             STRFPAK: STRF Estimation Software
% Copyright ©2003. The Regents of the University of California (Regents).
% All Rights Reserved.
% Created by Theunissen Lab and Gallant Lab, Department of Psychology, Un
% -iversity of California, Berkeley.
%
% Permission to use, copy, and modify this software and its documentation
% for educational, research, and not-for-profit purposes, without fee and
% without a signed licensing agreement, is hereby granted, provided that
% the above copyright notice, this paragraph and the following two paragr
% -aphs appear in all copies and modifications. Contact The Office of Tec
% -hnology Licensing, UC Berkeley, 2150 Shattuck Avenue, Suite 510,
% Berkeley, CA 94720-1620, (510) 643-7201, for commercial licensing
% opportunities.
%
%IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
%SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
%ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
%REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
%LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
%PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY,
%PROVIDED HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PRO
%-VIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
%  
% Created by JXZ, Mar. 19, 2003.

% Check box for user if the user accept copyright
answ = questdlg({...
    '             STRFPAK: STRF Estimation Software                          ',...
    ' Copyright 2003. The Regents of the University of California (Regents). ',...
    ' All Rights Reserved.                                                   ',...
    ' Created by Theunissen Lab and Gallant Lab, Department of Psychology at ',...
    ' University of California, Berkeley.                                    ',...
    '                                                                        ',...
    ' Permission to use, copy, and modify this software and its documentation',...
    ' for educational, research, and not-for-profit purposes, without fee and',...
    ' without a signed licensing agreement, is hereby granted, provided that ',...
    ' the above copyright notice, this paragraph and the following two       ',...
    ' paragraphs appear in all copies and modifications. Contact The Office  ',...
    ' of Technology Licensing, UC Berkeley, 2150 Shattuck Avenue, Suite 510, ',...
    ' Berkeley, CA 94720-1620, (510) 643-7201, for commercial licensing      ',...
    ' opportunities.                                                         ',...
    '                                                                        ',...
    '   Do you accept the above Copyright issue?                             ',...
    '                                                                        '},...
    'Copyright Message', 'Yes', 'No', 'Yes');
      
switch answ
case 'Yes'
        
    %         
    % Add matlab path
    addpath(pwd);

    clear all;
    % Run strfFirstGUI
    using_cache = allow_user_cache_input;
    
    if ~using_cache
        answ2 = questdlg({['You are not currently set up to automatically save re-usable calculations.' ...
            char(10)  'Would you like to start?' ]},'Start a cache?','Start Cache','Not Now','Start Cache');
        switch answ2
            case 'Start Cache'
                msgbox('First type "help Caching_Readme", then type "edit dir_of_caches.m" to get started with caching.');
                return
            otherwise
                strfFirstGUI;
        end
    end
    strfFirstGUI;
    
case 'No'
    msgbox('Bye Bye.');
    return
end


