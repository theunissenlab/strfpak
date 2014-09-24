function [results ] = read_strfpak_infoResult (infofilename)

    % first load that result file
    Rt = load(infofilename);

    totalNum = length(Rt.infopre);
    for ii = 1:totalNum
        
	    results(ii).infopre = Rt.infopre{ii};
	    results(ii).infouppre = Rt.infouppre{ii};
	    results(ii).infodownpre = Rt.infodownpre{ii};
        
	    results(ii).cc_spike_pre_constant = Rt.cc_spike_pre_constant{ii};
        results(ii).cc_spike_pre_best = Rt.cc_spike_pre_best{ii};
        
	    results(ii).cc_two_halves_constant= Rt.cc_two_halves_constant{ii};
        results(ii).cc_two_halves_best = Rt.cc_two_halves_tmax{ii};
        
	    results(ii).cc_ratio_max= Rt.cc_ratio_max{ii};
        
        results(ii).info = Rt.info{ii};
        results(ii).infoup = Rt.infoup{ii};
        results(ii).infodown = Rt.infodown{ii};
        
        results(ii).tmax_pre = Rt.tmax_pre{ii};
        results(ii).tmax_ratio = Rt.tmax_ratio{ii};
     end

 
   
