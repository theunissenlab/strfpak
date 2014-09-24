function [cc_best, cc_best_low, cc_best_high, cc_constant, cc_constant_low, cc_constant_high, corrval,...
        corrval_low, corrval_high,filterwidth] = cc_plot_est(spike1, spike2, widthvector, widthconstant, smoothflag)
% The function is going to return 10 arguments.
% function [cc_best, cc_constant, corrval, filterwidth] = cc_plot_est(
%          spike1, spike2, widthvector, widthconstant, smoothflag)
%   -- Calculate best correlation coefficients (cc), constant cc and cc
%   -- Find optimal filterwidth
%   Input:
%        spike1 -- one spike train or one psth
%        spike2 -- one spike train or one psth
%        widthvector -- provide filter width ranges
%        widthconstant -- constant width (default is 21)
%        smoothflag -- flag to smooth filter
%   Output:
%        cc_best -- best cc from given filter wide range
%        cc_constant -- constant cc with given constant filter width
%        filterwidth -- best filter width 
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

% Created by FET and ANNE
% Modified by JXZ, 2003
%

%
% Constants used in this function
i = 1;
corrmax = -1.0;
corrmax_low = -1.0;
corrmax_high = -1.0;

%
% Go through all the given width range and calcualte corrval
% and find optimal filterwidth
%
for width=widthvector(1):widthvector(2):widthvector(3)
    clear sest spre
    window_width = width;  
    wind1 = hanning(window_width)./sum(hanning(window_width));
    sest = conv(spike1,wind1);

    if ( smoothflag )  % Smooth both signals
        spre = conv(spike2,wind1);
    else
        % Do not smooth the second signal - the window is a delta function
        if mod(window_width,2)==0
            window_width=window_width+1; 
            wind2 = zeros(window_width,1);
            wind2((window_width+1)/2) = 1;
            wind2 = wind2(1:end-1);
        else
            wind2 = zeros(window_width,1);
            wind2((window_width+1)/2) = 1;
        end
        spre = conv(spike2,wind2);
    end
    [rcenter rp rlow rhigh] = corrcoef(spre,sest);
    
    corrval(i)= diag(rcenter,1);
    
    % Check for init or whatever...
    corrval_low(i) = diag(rlow,1);
    corrval_high(i) = diag(rhigh,1);
    
    % Negative correlation coefficients are sampling noise and set to zero
    if (corrval(i) < 0.0)
        corrval(i) = 0.0;   
        if (corrval_high(i) < 0.0)
            %print warning on significant negative correlation coefficient
            %warndlg('Warning: there is significant negative correlation coefficent.');
        end
    end
              
    
    % Calculate significance - correlation coefficient with p>0.05 are set
    % to zero. fet and jxz 08/17/05    
%     if ( diag(rp,1) < 0.05) 
%         corrval(i)=0.0;
%     end
    corrval_low(i) = max(0.0, corrval_low(i));
    corrval_high(i) = max(0.0, corrval_high(i));
    
    if  corrval(i) > corrmax 
        corrmax = corrval(i);
        corrmax_low = corrval_low(i);
        corrmax_high = corrval_high(i);
        filterwidth = width;
    end
    i = i+1;
end
corrval=corrval.^2;
corrval_low = corrval_low.^2;
corrval_high = corrval_high.^2;


cc_best=corrmax.^2;
cc_best_low = corrmax_low.^2;
cc_best_high = corrmax_high.^2;


%
% Now calculate constant_cc ADD SAME MODS DOWN BELOW
%

if mod(widthconstant,2)==0
    widthconstant=widthconstant-1;
    display('width changed to widthconstant-1');
end
wind1 = hanning(widthconstant)/sum(hanning(widthconstant));
wind2 = zeros(widthconstant,1);
wind2((widthconstant+1)/2) = 1;
sestconstant = conv(spike1,wind1);
spreconstant = conv(spike2,wind2);

% Modified on Dec. 11, 2003 
if ( smoothflag )
    spreconstant = conv(spike2,wind1);
else
    spreconstant = conv(spike2,wind2);
end

[cc_const, p_const, cc_low, cc_high] = corrcoef(sestconstant, spreconstant);
cc_const = max(0.0, diag(cc_const, 1));
cc_low = max(0.0, diag(cc_low,1));
cc_high = max(0.0, diag(cc_high,1));
% cc_constant=max(0.0,cc_constant);
% npts = length(sestconstant);
% tval = cc_constant*sqrt(npts-2)./sqrt(1-cc_constant^2);
% pval = tcdf(tval, npts-2);
% if (pval < 0.95)
%     cc_constant = 0.0;
% end

cc_constant=cc_const.^2;
cc_constant_low = cc_low.^2;
cc_constant_high = cc_high.^2;

% ======================================================
% END of cc_plot_est
% ======================================================


