function groupBMC_out=run_nfb_VBA_GroupBMC(L_struct)

groupBMC_in.L_SON1_admin1 = [];
groupBMC_in.L_SON1_admin2 = [];
groupBMC_in.L_SON2_adminPlac = [];
groupBMC_in.L_SON2_adminNalt = [];

%Get models inside our vba input structure
models = fieldnames(L_struct);
protocols = {'SON1', 'SON2'};
admins = {'admin1', 'admin2', 'adminPlac', 'adminNalt'};

for m = 1:length(models)
    for protocol = protocols
        for admin = admins
            try
                idx=~cellfun(@isempty,L_struct.(models{m}).(protocol{:}).(admin{:}).subj);
                groupBMC_in.(['L_' protocol{:} '_' admin{:}]) = [groupBMC_in.(['L_' protocol{:} '_' admin{:}]);...
                    L_struct.(models{m}).(protocol{:}).(admin{:}).L(idx)'];
            catch
                
            end
        end
    end
end

analysis_fnames = fieldnames(groupBMC_in);
for analysis_fname = analysis_fnames'
    [post,out]=VBA_groupBMC(groupBMC_in.(analysis_fname{:}));
    groupBMC_out.(analysis_fname{:}).posterior = post;
    groupBMC_out.(analysis_fname{:}).out = out;
end

groupBMC_out.model_order = models;

