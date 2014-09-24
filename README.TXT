
STRFPAK (a spatio-temporal receptive field estimation software)
==============================================================

Version 4.1.2 - August, 2006
------------------------

Description
-----------
STRFPAK is a Matlab toolbox for estimating the linear and nonlinear 
(future versions) stimulus-response transfer function of a sensory neuron.
The resulting spatio-temporal receptive field (STRF) provides a quantitative 
description of neural filtering properties that can be used in subsequent 
computational modeling studies. The estimation techniques implemented by 
STRFPAK are quite general. Several algorithms are provided for estimating 
both linear and nonlinear STRFs from responses to either simple or complex 
stimuli, including natural signals.

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
%PROVIDED HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO 
%PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

Installation
------------

1. Unzip STRFPAK-4.1.2 and put all the .m files in the same directory
The unzip methods depends on the platform you are using. See detail
on STRFPAK User's manual.

How to use STRF package
---------------------------------
1. Start to run Matlab first.

2. On the MATLAB command prompt line, type:
   >> strfpak

3. Follow the flow on the main window of STRFPAK to get input data,
   estimate the kernel, predict PSTH using estimated kernel and 
   validate the goodness of fit.

4. All the intermediate results are saved under the directory you
   specified when you click "Calculate" button on the main window.

5. If you have done all the calculation, you want to display them
   from the main window. You need first specify the global variable
   "outputPath" at the matlab command line:
>> global outputPath
>> outputPath = '/home/junli/STRFPAK/Vision_Output';

 
Main intermediate output file list
-----------------
Stim_autocorr.mat               : stimuli autocorrelation matrix
StimResp_crosscorr.mat          : stimuli-response crosscorrelation vector
SR_crosscorrJN.mat              : Jackknifed version crosscorrelation
strfResult.mat                  : estimated STRF results
predResult-EstSpike_Tolx.mat    : predicted psth for tolrence value x 
predResult-avgSpike_Tolx.mat    : original psth used for prediciton at tol1 
info_r_result.mat               : validation results

Main files list
----------------
1. Calculation Functions:
----------------
cal_AVG.m      : Calculate average of stimuli that used for removing noise
cal_AutoCorr.m : Calculate autocorrelation of stimuli
cal_AutoCorrSep.m  : Calculate autocorrelation of stimuli using space-temp 
                     separable algorithm
cal_CrossCorr.m    : Calculate crosscorrelation of stimuli and response
cal_Strf.m         : Calculate strf 
cal_StrfSep.m      : Calculate strf using separable algorithm 
calStrf_script.m   : Calcualte strf and normalization
calStrfSep_script.m: Calcualte strf and normalization using separable alg.
cal_PredStrf.m     : Predict psth using estimated STRF and new dataset
cal_Validate.m     : Validate its prediction

2. Display Functions:
-----------------
displayxxxxx.x    : Display routines for xxxxx

3. Helper Functions:
------------------
rest

Questions and Comments
----------------------
Please contact Patrick Gill for your questions and comments.
E-Mail:   patgill@berkeley.edu 
WWW:      http://strfpak.berkeley.edu

