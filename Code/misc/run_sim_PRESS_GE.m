function run_sim_PRESS_GE
% all_mets = {'Ala','Asc','Asp','bHB','bHG','Cit','Cr','EtOH','GABA','GPC','GSH','Glc','Gln' ...
%     ,'Glu','Gly','H2O','Ins','Lac','NAA','NAAG','PCh','PCr','PE','Phenyl' ...
%     ,'Scyllo','Ser','Tau','Tyros'};

all_mets = {'GABA'};

for kk = 1:length(all_mets)
    [out_2] = sim_PRESS_fast_fastphase_TE30_GE({all_mets{kk}});
end

end

