function run_nfb_vba
%Parent script to run subjects through VBA data analysis scripts

%% Load in main data structure
if ispc
    load('E:\Box Sync\fMRI_shared\NFB_response\nfball.mat') %Called out
else
    load('/Users/martapecina/Box Sync/PITT/RESEARCH/fMRI_shared/NFB/NFB_response/SON1&2_behav_results/nfball.mat')
end

%Rename nfball to reduce confusion from VBA output, clear original out 
nfball=out;
clear out

%% Organize data further
% We may want to take care of this step in the compiler script, but for
% now let's loop over subjects and admins for SON1, SON 2 will be another
% day

%For SON1
subj_names = fieldnames(nfball.SON1);
for ii = 1:length(subj_names)
   admins = fieldnames(nfball.SON1.(subj_names{ii}));
   for jj = 1:length(admins)
       
       if ~strcmpi(admins{jj},'admin1') && ~strcmpi(admins{jj},'admin2')
           continue
       end
       
      vba_input=struct; %Initialize the input variable for the vba script      
      
      %Get and set the data
      data = nfball.SON1.(subj_names{ii}).(admins{jj}); 
      vba_input.data = data;      
      
      vba_input.admin = admins{jj}; %Get administration
      vba_input.subj_name = subj_names{ii}; %Get subj name
      %TODO graphics, multinomial, multisession,ect
      try
        [posterior,out]=pavlov_vba_sonsira(vba_input);
        L(ii,jj) = out.F;
      catch
          fprintf('Subject %s died...\n\n',vba_input.subj_name)
          L(ii,jj) = nan;
      end
   end
end
