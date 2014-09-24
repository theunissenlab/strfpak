% This file extracts the CC and info values at best tolerance for a selected data set

clear all
close all

data_files_mld = {
   
        '/auto/fdata/junli/PPdata_Nov5/blublu0916/1_A/' ...
        '/auto/fdata/junli/PPdata_Nov5/oo0108/4_A/' ...
        '/auto/fdata/junli/PPdata_Nov5/blabla0713/3_A/' ...
        '/auto/fdata/junli/PPdata_Nov5/blabla1904/1_A/' ...
        '/auto/fdata/junli/PPdata_Nov5/blabla1904/2_A/' ...
        '/auto/fdata/junli/PPdata_Nov5/blabla1904/5_A/' ...
        '/auto/fdata/junli/PPdata_Nov5/oo1010/3_A/' ...
        '/auto/fdata/junli/PPdata_Nov5/oo0909/4_A/' ...
        '/auto/fdata/junli/PPdata_Nov5/pupu0333/5_A/' ...
};

data_files_fieldl = {
    
        '/auto/fdata/junli/PPdata_Nov5/rr3334/1_B/' ...
        '/auto/fdata/junli/PPdata_Nov5/pupu0606/2_B' ...
        '/auto/fdata/junli/PPdata_Nov5/yy1617/1_B' ...
        '/auto/fdata/junli/PPdata_Nov5/yy1617/2_B' ...
        '/auto/fdata/junli/PPdata_Nov5/yy1617/3_B' ...
        '/auto/fdata/junli/PPdata_Nov5/yy1617/5_B' ... 
        '/auto/fdata/junli/PPdata_Nov5/blabla0713/3_B' ... 
        '/auto/fdata/junli/PPdata_Nov5/blabla0713/4_B' ... 
        '/auto/fdata/junli/PPdata_Nov5/blabla1904/3_B' ... 
        '/auto/fdata/junli/PPdata_Nov5/blublu0809/4_B' ... 
        '/auto/fdata/junli/PPdata_Nov5/pupu0333/1_B' ... 
        '/auto/fdata/junli/PPdata_Nov5/pupu0333/2_B' ...
        
};

start_cd = pwd;
nfiles = length(data_files_mld);
for ifiles=1:nfiles
    cd(data_files_mld{ifiles});
    
    if  exist(fullfile(pwd,'con_strfpak/Output/info_r_result.mat')) == 2
     
        
        filt_allres = read_strfpak_infoResult(fullfile(pwd,'con_strfpak/Output/info_r_result.mat'));
        [infomax indmax]=max([filt_allres.infopre]);
        filt_goodtol_mld_con(ifiles) = filt_allres(indmax);  
        
    else
        disp(['Please redo ',fullfile(pwd,'con_strfpak')]);
    end
 
    
    if  exist(fullfile(pwd,'songrip_strfpak/Output/info_r_result.mat')) == 2
        
       
        filt_allres=read_strfpak_infoResult(fullfile(pwd,...
                       'songrip_strfpak/Output/info_r_result.mat'));
        [infomax indmax]=max([filt_allres.infopre]);
        filt_goodtol_mld_songrip(ifiles) = filt_allres(indmax);   
        
    else
        disp(['Please redo ',fullfile(pwd,'songrip_strfpak')]);
    end
    
    
    
    if  exist(fullfile(pwd,'flatrip_strfpak/Output/info_r_result.mat')) == 2
        
       
        filt_allres=read_strfpak_infoResult(fullfile(pwd,...
                      'flatrip_strfpak/Output/info_r_result.mat'));
        [infomax indmax]=max([filt_allres.infopre]);
        filt_goodtol_mld_flatrip(ifiles) = filt_allres(indmax);   
        
    else
        disp(['Please redo ',fullfile(pwd,'flatrip_strfpak')]);
      
    end
    
end
cd (start_cd);

nfiles = length(data_files_fieldl);
for ifiles=1:nfiles
    cd(data_files_fieldl{ifiles});
    
    if  exist(fullfile(pwd,'con_strfpak/Output/info_r_result.mat')) == 2
     
        
        filt_allres=read_strfpak_infoResult(fullfile(pwd,...
                          'con_strfpak/Output/info_r_result.mat'));
        [infomax indmax]=max([filt_allres.infopre]);
        filt_goodtol_fieldl_con(ifiles) = filt_allres(indmax);  
    else
        disp(['Please run STRFPAK on ',fullfile(pwd,'con_strfpak'),' first!']);
        
    end
    
    
    
    
    if  exist(fullfile(pwd,'songrip_strfpak/Output/info_r_result.mat')) == 2
        
       
        filt_allres=read_strfpak_infoResult(fullfile(pwd,...
                       'songrip_strfpak/Output/info_r_result.mat'));
        [infomax indmax]=max([filt_allres.infopre]);
        filt_goodtol_fieldl_songrip(ifiles) = filt_allres(indmax);   
    else
        disp(['Please redo ',fullfile(pwd,'songrip_strfpak')]);
        
    end
    
    
    
    if  exist(fullfile(pwd,'flatrip_strfpak/Output/info_r_result.mat')) == 2
        
       
        filt_allres=read_strfpak_infoResult(fullfile(pwd,...
                      'flatrip_strfpak/Output/info_r_result.mat'));
        [infomax indmax]=max([filt_allres.infopre]);
        filt_goodtol_fieldl_flatrip(ifiles) = filt_allres(indmax);   
    else
        disp(['Please redo ',fullfile(pwd,'flatrip_strfpak')]);
        
      
    end
   
    
end

cd(start_cd);

fprintf(1, 'STRFPAK: Mean predicted information for MLd Con: %f\n', mean([filt_goodtol_mld_con.infopre]));    
fprintf(1, 'STRFPAK: Mean predicted information for Field L Con: %f\n', mean([filt_goodtol_fieldl_con.infopre]));   

fprintf(1, 'STRFPAK: Mean predicted cc_ratio_max for MLd Con: %f\n', mean([filt_goodtol_mld_con.cc_ratio_max]));
fprintf(1, 'STRFPAK: Mean predicted cc_ratio_max for Field L Con: %f\n', mean([filt_goodtol_fieldl_con.cc_ratio_max])); 

fprintf(1, 'STRFPAK: Mean cc_spike_pre_constant  for MLd Con: %f\n', mean([filt_goodtol_mld_con.cc_spike_pre_constant]));    
fprintf(1, 'STRFPAK: Mean cc_spike_pre_constant for Field L Con: %f\n', mean([filt_goodtol_fieldl_con.cc_spike_pre_constant]));
fprintf(1, 'STRFPAK: Mean cc_two_halves_constant  for MLd Con: %f\n', mean([filt_goodtol_mld_con.cc_two_halves_constant]));    
fprintf(1, 'STRFPAK: Mean cc_two_halves_constant for Field L Con: %f\n', mean([filt_goodtol_fieldl_con.cc_two_halves_constant]));

fprintf(1, 'STRFPAK: Mean cc ratio constant for MLd Con: %f\n', mean([filt_goodtol_mld_con.cc_spike_pre_constant]./[filt_goodtol_mld_con.cc_two_halves_constant]));    
fprintf(1, 'STRFPAK: Mean cc ratio constant for Field L Con: %f\n', mean([filt_goodtol_fieldl_con.cc_spike_pre_constant]./[filt_goodtol_fieldl_con.cc_two_halves_constant]));    
