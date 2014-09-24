function [strf,best_tol_index,best_std_index,savefile] = copy_best_strf;
global outputPath Tol_val Std_val
the_r = load(fullfile(outputPath,'info_r_result.mat'));
infopre = the_r.infopre;
best_info = -999;
best_tol = -999;
best_tol_index = -999;
best_std = -999;
for jj = 1:length(Tol_val)
    for kk = 1:length(Std_val)
        if infopre{jj}{kk} > best_info
            best_info = infopre{jj}{kk};
            best_tol = Tol_val(jj);
            best_tol_index = jj;
            best_std = Std_val(kk);
            best_std_index = kk;
        end
    end
end
loaded = load(fullfile(outputPath,['strfResult_Tol' num2str(best_tol_index) '.mat']));
strf = fast_filter_filter(loaded.STRF_Cell,mean(loaded.STRFJNstd_Cell,3),best_std);
savefile = fullfile(outputPath,'best_strf.mat');
save(savefile,'strf','best_tol','best_std');
disp(['Saving best strf in file ' savefile '.']);
