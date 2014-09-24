% Loads and process the stim - Spike cross-correlation
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

stim_spike = Check_And_Load(crossfile);

disp('Done loading stim_spike cross correlatin.');

global NBAND
global TimeLag

nb = size(stim_spike, 1);
nt = size(stim_spike, 2);
t = (nt -1)/2;
t=-t:t;

f= 1:nb;

minval = min(min(stim_spike));
maxval = max(max(stim_spike));
axisval(1)=t(1);
axisval(2)=t(nt);
axisval(3)=minval;
axisval(4)=maxval;

% two-dimension display
splitX = floor(sqrt(nb));
STA = reshape(stim_spike, splitX, splitX, nt);
amax = max(abs(STA(:)));
amin = -amax;
   
numList = 12;
if numList > nt
    numList = nt
end
    
titlestr = 'Lag=';
% display total 12 images
for fr = 1:numList
      
        h=subplot(2,6,fr);
      
        if amin~=amax,
            imagesc(STA(:,:,fr),[amin,amax]);
        else
            imagesc(zeros(size(STA,1),size(STA,2)));
        end
      
        if size(STA,1)>1 & size(STA,2)>1,
            axis image
        end
      
        set(h,'YTickLabel',[]);
        set(h,'XTickLabel',[]);
      
         
        curTitle = strcat(titlestr, num2str(fr-1));
        title(curTitle);
        colormap(redblue);
end


% ===============================================================
% END of plot_stimspikeCross2d
% ===============================================================
