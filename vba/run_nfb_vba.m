 function L_struct=run_nfb_vba
%Parent script to run subjects through VBA data analysis scripts

%% Load in main data structure
if ispc
    load('E:\Box Sync\fMRI_shared/NFB/NFB_response/SON1&2_behav_results/nfball.mat') %Called out
else
    load('/Users/martapecina/Box Sync/PITT/RESEARCH/fMRI_shared/NFB/NFB_response/SON1&2_behav_results/nfball.mat')
end

%Rename nfball to reduce confusion from VBA output, clear original out
nfball=out;
clear out

%% Organize data further
% We may want to take care of this step in the compiler script, but for
% now let's loop over subjects and admins for SON1 and SON 2

%For SON1
%protocol = {'SON1', 'SON2'};
protocol = {'SON1','SON2'};
%model_names = {'null','oneLR_twoK','oneLR_fixD_twoK','twoLR_twoK','twoLR_fixD','twoLR_S_twoK','twoLR_S_fixD_twoK','twoLR_S_fixD_oneK','oneLR_S_twoK','oneLR_S_fixD_twoK','twoLR_S_fixD_m_twoB_twoK','momentum','momentum_oneK'};
model_names = {'null','twoLR_S_fixD_oneK'};


for model_name = model_names
    for kk = 1:length(protocol)
        subj_names = fieldnames(nfball.(protocol{kk}));
        for ii = 1:length(subj_names)
            admins = fieldnames(nfball.(protocol{kk}).(subj_names{ii}));
            for jj = 1:length(admins)
                
                if ~strcmpi(admins{jj},'admin1') && ~strcmpi(admins{jj},'admin2') && ~strcmpi(admins{jj},'adminNalt') && ~strcmpi(admins{jj},'adminPlac')
                    continue
                end
                
                %filter any "bad" subjects
                % % %             if any(ismember(bad_subjs,subj_names{ii}))
                % % %                 continue; %Skip them
                % % %             end
                
                vba_input=struct; %Initialize the input variable for the vba script
                
                %Get and set the data
                vba_input.data = nfball.(protocol{kk}).(subj_names{ii}).(admins{jj});
                
                %vba_input.admin = admins{jj}; %Get administration
                vba_input.subj_name = subj_names{ii}; %Get subj name
                vba_input.admin = admins{jj}; %Get the current administration
                vba_input.protocol = protocol{kk};
                
                %vba_input.subj_name = subj_names{8}; %Get one subj name
                %TODO graphics, multinomial, multisession,ect
                %       try
                
                fprintf('Fitting %s %s from %s using VBA now \n\n',protocol{kk}, vba_input.subj_name, admins{jj})
                
                [posterior,out]=pavlov_vba_sonrisa(vba_input,model_name{:});
                
                %Save the L as a struct which we can then use for model
                %comparisons
                try
                    L_struct.(model_name{:}).(protocol{kk}).(admins{jj});
                catch
                    L_struct.(model_name{:}).(protocol{kk}).(admins{jj}) = table();
                end
                L_struct.(model_name{:}).(protocol{kk}).(admins{jj}).subj(ii,1) = {subj_names{ii}};
                L_struct.(model_name{:}).(protocol{kk}).(admins{jj}).L(ii,1) = out.F;
                %L_struct.(model_name{:}).(protocol{kk}).(admins{jj}) = [L_struct.(model_name{:}).(protocol{kk}).(admins{jj}) out.F]; %Output of log model evidence
                
                
                fprintf('----------------------------------------------\n\n')
                %       catch
                %           fprintf('Subject %s died...\n\n',vba_input.subj_name)
                %           L(ii,jj) = nan;
                %       end
            end
        end
    end
end



%Save the output
time_complete = datetime;
save('vba/L_struc/L_struct.mat','L_struct','time_complete')
