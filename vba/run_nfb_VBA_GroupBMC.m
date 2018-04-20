function groupBMC_out=run_nfb_VBA_GroupBMC(L_struct)

groupBMC_in.L_SON1_admin1 = [];
groupBMC_in.L_SON1_admin2 = [];
groupBMC_in.L_SON2_adminPlac = [];
groupBMC_in.L_SON2_adminNalt = [];

%Get models inside our vba input structure
models = fieldnames(L_struct);
protocols = {'SON1', 'SON2'};
%protocols = {'SON1'};
admins = {'admin1', 'admin2', 'adminPlac', 'adminNalt'};
%admins = {'admin1', 'admin2'};

%Censor ids here as a cell with the format 'subjxxx'
censor.son1_censor = {};
censor.son2_censor = {};




for m = 1:length(models)
    for protocol = protocols
        for admin = admins
            try
                idx=~cellfun(@isempty,L_struct.(models{m}).(protocol{:}).(admin{:}).subj);
                
                %Censor our specific subjects if any
                if strcmp(protocol{:},'SON1')
                    idx = censor_subjects(idx, L_struct.(models{m}).(protocol{:}).(admin{:}).subj,censor.son1_censor); 
                else
                    idx = censor_subjects(idx, L_struct.(models{m}).(protocol{:}).(admin{:}).subj,censor.son2_censor); 
                end
                
                %Compile input for BMC
                groupBMC_in.(['L_' protocol{:} '_' admin{:}]) = [groupBMC_in.(['L_' protocol{:} '_' admin{:}]);...
                    L_struct.(models{m}).(protocol{:}).(admin{:}).L(idx)'];
                subject_idx.(['L_' protocol{:} '_' admin{:}]) = L_struct.(models{m}).(protocol{:}).(admin{:}).subj(idx);
            catch
                
            end
        end
    end
end

%Set up analysis loop
analysis_fnames = fieldnames(groupBMC_in);

%Censor subjects here

% % % censor_admins = fieldnames(censor);
% % % for censor_admin = censor_admins'
% % %     if isempty(censor.(censor_admin{:}))
% % %         continue %Skip if no censor values
% % %     end
% % %     stop=1;
% % %     ids_to_censor = censor.(censor_admin{:});
% % %     
% % %     
% % %     
% % % end

%Start analysis loop
for analysis_fname = analysis_fnames'
    [post,out]=VBA_groupBMC(groupBMC_in.(analysis_fname{:}));
    groupBMC_out.(analysis_fname{:}).posterior = post;
    groupBMC_out.(analysis_fname{:}).out = out;
end

groupBMC_out.model_order = models;


function idx=censor_subjects(idx,subj_list,censor_list)
if ~isempty(censor_list)
    tmp_cen_idx = cellfun(@(x) strcmp(x,subj_list), censor_list, 'UniformOutput',0);
    censor_idx=sum(horzcat(tmp_cen_idx{:}),2);
    idx = idx & ~logical(censor_idx);
else
    idx = idx;
end