function [posterior,out]=run_single_subject_sonrisa_vba()

load('/Users/martapecina/Box Sync/PITT/RESEARCH/fMRI_shared/NFB/NFB_response/SON1&2_behav_results/nfball.mat')
nfball = out;
protocol = {'SON1'};%SON1 or SON2
subj_name = {'subj017'}; %ex: subj001
admin = {'admin1'}; %son1 = admin1 or admin2 son2 = adminPlac or adminNalt
vba_input.data = nfball.(protocol{:}).(subj_name{:}).(admin{:});
vba_input.protocol=protocol{:};
vba_input.admin=admin{:};
vba_input.subj_name = subj_name{:};
model_name = 'twoLR_S_fixD_oneK';

[posterior,out]=pavlov_vba_sonrisa(vba_input, model_name);