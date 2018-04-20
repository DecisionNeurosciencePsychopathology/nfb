function out=compile_nfb_data_1_fsl(regressor_flag)
%Funciton will organize and create simple plots of nfb task data. For more
%information see README at https://github.com/heffjos/nfb.
%
%To create regressors for all subjects in both studies use the command:
%out=compile_nfb_data_1(1)

%For now...
global gen_concat_reg final_recorded_time

gen_concat_reg = 0; %1==AFNI 0==FSL

%%Create regreossors or not
try regressor_flag; catch, regressor_flag=0; end

%% Datalocation -- set as needed
if ispc
    data_path = 'E:\Box Sync\fMRI_Shared\NFB\NFB_response\'; %Change as needed
    vba_data_path = ''
    %For the time being because depending on the system we might not have a
    %'/' in front of the string
    str_direction=1;
else
    data_path = '/Users/martapecina/Box Sync/PITT/RESEARCH/fMRI_shared/NFB/NFB_response/';
    vba_data_path = '/Users/martapecina/GitHub/Nfb_sonrisa/vba_data/';
    str_direction=-1;
end

%% Initialize storgae variables
out.SON1 = struct;
out.SON2 = struct;
out.SON1_MDF_1 = [];
out.SON1_MDF_2 = [];
out.SON2_MDF_Nalt = [];
out.SON2_MDF_Plac = [];

%%Additional options
options.rflag=regressor_flag;
options.str_direction = str_direction;
% Compile loop

%Grab parent directories
experiment_paths=glob([data_path 'SON*'])';

for experiment_path = experiment_paths
    
%     if strfind(experiment_path{:},'SONRISA1')
%         experiment_path = experiment_path{:};
%     end
    
    %Grab dir list from datalocation
    expression = '.*\d{3,9}.*';
    data_dirs = regexp(glob([experiment_path{:} '/SON*']), expression, 'match'); %Think of a way to just get the uniue ids here , probably with cellfun
    data_dirs = [data_dirs{:}];
    
    for data_dir = data_dirs
        %Put everything into one table
        out=compile_single_subject_data(out,data_dir,options);
    end
end

%% Compile all sonrisa data into two csvs one for SON1 and one for SON2

%Load in lookup table for SON2 to SON1 migration
[~,~,lkup_tbale]=xlsread([data_path 'SON1&2_behav_results/son2_to_son1_lkup_table.xlsx']);
out = merge_son2_plac_to_son1(out,lkup_tbale);

%Create the two seperate csvs
%Compile
son1_all_tmp = [out.SON1_MDF_1; out.SON1_MDF_2];
son2_all_tmp = [out.SON2_MDF_Nalt; out.SON2_MDF_Plac];

%Sort
son1_all_tmp = sortrows(son1_all_tmp,'Participant','ascend');
son2_all_tmp = sortrows(son2_all_tmp,'Participant','ascend');
sortedNames_son1 = sort(son1_all_tmp.Properties.VariableNames(2:end));
sortedNames_son2 = sort(son2_all_tmp.Properties.VariableNames(2:end));
out.son1_all = [son1_all_tmp(:,1) son1_all_tmp(:,sortedNames_son1)];
out.son2_all = [son2_all_tmp(:,1) son2_all_tmp(:,sortedNames_son2)];

%Write
writetable(out.son1_all,[data_path 'SON1&2_behav_results/son1_all.csv']);
writetable(out.son2_all,[data_path 'SON1&2_behav_results/son2_all.csv']);

%% Save final output
save([data_path 'SON1&2_behav_results/nfball'],'out')


function out=compile_single_subject_data(out,data_dir,options)
%Function will go though each directory in the predefined data location
%and compile single subjects into one csv, then using the single subject
%csvs, compile a group level csv for each task (SON1-2) & administration
%(baseline (1) vs 8 weeks later (2))

%Get the subj id
expression = '\d{3,9}';
id=regexp(data_dir,expression,'match');
id=id{:};

%Get the administration code from the tail of the dir - Null is SON2
expression = '_(\d)\\|\/$';
admin_code=regexp(data_dir,expression); %Not sure why the grouping didn't work...
admin_code=data_dir{:}(admin_code{:}+options.str_direction);

%For SON-2 designate the type of run (Nalt or Plac)
expression = '([A-Za-z]{4})(?:[\\|\/])?$';
run_str=regexp(data_dir,expression,'tokens'); 
%Change into a function later
while iscell(run_str) && ~isempty(run_str)
    run_str=run_str{:};
end 

%Get the csv files in the directory (they should start with SON1 or SON2)
expression = '(?!.*[\\|\/])SON.*Run.*csv';
files=regexp(glob([data_dir{:} '*csv']),expression,'match');
files=[files{:}]; %May have to sort these if it doesn't do it automatically

%Check if all runs are present if not create missing data run files
% if length(files)~=4
%     files=file_checker(data_dir{:},files);
% end

%Initialize table var
T=[];

%File loop
for file = files
    T=[T; readtable([data_dir{:} file{:}])];
end

%Fill any missing runs with nans for completion purposes
if length(files)~=4
    %We are assuming 4 runs here
    runs=1:4;
    missing_runs = runs(~ismember(runs,unique(T.Run))); %Which missing runs
    vars=T.Properties.VariableNames(10:end); %Keep some vars
    for missing_run=missing_runs
        tmp_T=readtable([data_dir{:} file{:}]);
        tmp_T.Run=ones(length(tmp_T.Run),1).*missing_run;
        for var=vars
            if iscell(tmp_T.(var{:}))
                tmp_T.(var{:})=repmat({'Nan'},length(tmp_T.Run),1);
            else
                tmp_T.(var{:})=nan(length(tmp_T.Run),1);
            end
        end
        %update the main subject table
        T=[T; tmp_T];
    end
end

%Numerate responses
T=numerate_responses(T);

%Decide which struct to use & write the data to the subj lvl table
if strfind(data_dir{:},'SON1')
    out.SON1.(['subj' id{:}]).(['admin' admin_code]) = T;
    proto={'SON1'};
    suffix_codes={'a','b'};
    reg_suffix = suffix_codes{str2double(admin_code)};
    str = 'SON1_MDF_';
elseif strfind(data_dir{:},'SON2')
    admin_code = run_str; %Override admin code
    out.SON2.(['subj' id{:}]).(['admin' admin_code]) = T;
    proto={'SON2'};
    reg_suffix=run_str;
    str = 'SON2_MDF_';
else
    error('Something went wrong')
end

%Save the administration
T.administration = repmat(admin_code,height(T),1);

%Parse out the data for easier regressor creation
%out.(proto{:}).(['subj' id{:}]).(['admin' admin_code '_data'])=parse_single_subj_data(out.(proto{:}).(['subj' id{:}]).(['admin' admin_code]),proto{:},id,reg_suffix);
tmp_data=parse_single_subj_data(out.(proto{:}).(['subj' id{:}]).(['admin' admin_code]),proto{:},id,reg_suffix);
out.(proto{:}).(['subj' id{:}]).(['admin' admin_code '_reg_id']) = tmp_data.reg_id;

%Make sure there aren't any trailing letters on the subj ids
T.Participant = repmat({[proto{:} '_' id{:}]},height(T),1);
id_tmp=regexp(id{:},'[1-9]+.*','match');
id_tmp=id_tmp{:};
T.subject_id = repmat({id_tmp},height(T),1);

%TODO: fix this to make it more generic allowing for differnt models to be
%easily chosen, also spruce up the loading of vba parameters
model = 'twoLR_S_fixD_oneK';
vba_file = glob(sprintf('/Users/martapecina/GitHub/Nfb_sonrisa/vba_data/%s/subj_subj%s_admin%s_vba_data_*.mat',model,id{:},admin_code));
%vba_file = glob(sprintf('/Users/martapecina/GitHub/Nfb_sonrisa/vba_data/exp_yoked/subj_subj%s_admin%s_vba_data_*.mat',id{:},admin_code));
if ~isempty(vba_file)
    vba_data = load(vba_file{:},'out','posterior');
    close all; %In case figures got saved too...
    [r,c] = size(vba_data.posterior.muX);
    value_logical_idx = zeros(r,c);
    for i = 1:length(value_logical_idx)
        current_stim = vba_data.out.u(4,i);
        value_logical_idx(current_stim,i) = 1;
    end
    T.value_per_stim = vba_data.posterior.muX(logical(value_logical_idx));
    T.PE = vba_data.posterior.muX(end,:)';
end


%Update master list
try
    out.([str admin_code]) = [out.([str admin_code]); T];
catch
    %In case the columns don't match up replace missing vars with nans, I
    %don't think this should overwrite any columns but be careful here
    T_colmissing = setdiff(out.([str admin_code]).Properties.VariableNames, T.Properties.VariableNames);
    out_colmissing = setdiff(T.Properties.VariableNames, out.([str admin_code]).Properties.VariableNames);
    T = [T array2table(nan(height(T), numel(T_colmissing)), 'VariableNames', T_colmissing)];
    out.([str admin_code]) = [out.([str admin_code]) array2table(nan(height(out.([str admin_code])), numel(out_colmissing)), 'VariableNames', out_colmissing)];
    
    %If the data doesn't want to merge becasue if data type
    for colname = T_colmissing
        if iscell(out.([str admin_code]).(colname{1}))
            T.(colname{1}) = cell(height(T), 1);
        end
    end
    for colname = out_colmissing
        if iscell(T.(colname{1}))
            out.([str admin_code]).(colname{1}) = cell(height(out.([str admin_code])), 1);
        end
    end
    
    %Now update
    out.([str admin_code]) = [out.([str admin_code]); T];
end


%Hack be wary
out.(proto{:}).(['subj' id{:}]).(['admin' admin_code]) = T;

%Save new file in data location
writetable(T,[data_dir{:} sprintf('subj_%s_all_runs.csv',id{:})])

%Create regressors if needed -- this could be moved out of here if need be
if options.rflag
    %write_nfb_regressors_2(out.(proto{:}).(['subj' id{:}]).(['admin' admin_code]),out.(proto{:}).(['subj' id{:}]).(['admin' admin_code '_reg_id']))
    write_nfb_regressors_2(T,out.(proto{:}).(['subj' id{:}]).(['admin' admin_code '_reg_id']))
end



function foo=file_checker(data_dir,files)

function ss_data=parse_single_subj_data(T,proto,id,suffix_code)

%Grab the column names
col_names=T.Properties.VariableNames;

%Set subject id name as we need it to look for the afni models Ex:
%SON1_001_a for subject 001 admin 1
ss_data.reg_id = sprintf('%s_%s_%s',proto,id{:},suffix_code);

%Parse out the table data into a struct
% for i=2:length(col_names)
%     ss_data.(col_names{i})=T.(col_names{i});
% end


function T=numerate_responses(T)

T.WillImpRespBin=double(cellfun(@(x) strcmpi(x,'yes'), T.WillImpRespText));
T.ImprovedRespBin=double(cellfun(@(x) strcmpi(x,'yes'), T.ImprovedRespText));

%Take care of nans
T.WillImpRespBin(cellfun(@(x) strcmpi(x,'NaN'), T.WillImpRespText))=nan;
T.ImprovedRespBin(cellfun(@(x) strcmpi(x,'NaN'), T.ImprovedRespText))=nan;

function data = merge_son2_plac_to_son1(data,lkup_table)
son_2_subjs = fieldnames(data.SON2);
for son_2_subj = son_2_subjs'
    
    %Pull the numerical id
    son_2_subj=regexp(son_2_subj{:},'\d+','match');
    
    %If subject from SON2 is found in the lookup table rip the placebo data
    %as SON1 administration 
    %N.B. Nans may creep into the lookup table remove them manually or
    %write a snippet.
    [~, lkup_idx]=ismember(son_2_subj,lkup_table(:,1));
    if lkup_idx>0
        %Set SON1 id
        son_1_id = lkup_table{lkup_idx,2};
        
        %Grab the SON@ subj's placebo data 
        tmp_data=data.SON2.(['subj' son_2_subj{:}]).adminPlac;
        tmp_data.Participant = repmat({['SON1_' son_1_id]},height(tmp_data),1); %Replace id (must be cell?)
        tmp_data.administration = repmat('1',height(tmp_data),1); %Replace admin code -- this should only be the first administration
        tmp_id = regexp(son_1_id,'[1-9]+.*','match');
        tmp_id=tmp_id{:};
        tmp_data.subject_id = repmat({tmp_id},height(tmp_data),1);
        
        %Add to SON1 datastruct
        data.SON1.(['subj' son_1_id]).admin1 = tmp_data;
        data.SON1.(['subj' son_1_id]).admin1_reg_id = ['SON1_' son_1_id '_a']; %add reg id
        
        %Add to SON1 admin 1 dataframe -- This is ugly fix it up
        try
            data.SON1_MDF_1 = [data.SON1_MDF_1; tmp_data];
        catch
            tmp_data_colmissing = setdiff(data.SON1_MDF_1.Properties.VariableNames, tmp_data.Properties.VariableNames);
            data.SON1_MDF_1_colmissing = setdiff(tmp_data.Properties.VariableNames, data.SON1_MDF_1.Properties.VariableNames);
            tmp_data = [tmp_data array2table(nan(height(tmp_data), numel(tmp_data_colmissing)), 'VariableNames', tmp_data_colmissing)];
            data.SON1_MDF_1 = [data.SON1_MDF_1 array2table(nan(height(data.SON1_MDF_1), numel(data.SON1_MDF_1_colmissing)), 'VariableNames', data.SON1_MDF_1_colmissing)];
            data.SON1_MDF_1 = [data.SON1_MDF_1; tmp_data];
        end
    end
    
end




% % % %Grab field names
% % % fnames = fieldnames(data);
% % % for i = 1:length(fnames)
% % %     if strfind(fnames{i},'MDF')>0 %If it is a master data file
% % %         data.(fnames{i}).WillImpRespNum=double(cellfun(@(x) strcmpi(x,'yes'), data.(fnames{i}).WillImpRespText));
% % %         data.(fnames{i}).ImprovedRespNum=double(cellfun(@(x) strcmpi(x,'yes'), data.(fnames{i}).ImprovedRespText));
% % %         
% % %         %Take care of nans
% % %         data.(fnames{i}).WillImpRespNum(cellfun(@(x) strcmpi(x,'NaN'), data.(fnames{i}).WillImpRespText))=nan;
% % %         data.(fnames{i}).ImprovedRespNum(cellfun(@(x) strcmpi(x,'NaN'), data.(fnames{i}).ImprovedRespText))=nan;
% % %     end
% % %     
% % % end
% % %         
        



function write_nfb_regressors_2(data,id)

global final_recorded_time

%% Datalocation -- set as needed
if ispc
    reg_file_dest = 'C:/kod/nfb/analysis/regs/'; %Change as needed
else
    reg_file_dest = '/Users/martapecina/GitHub/Nfb_sonrisa/regs/';
end


%Pull the final recorded times
runs = unique(data.Run);
for run = 1:max(runs)
   final_recorded_time(run) = data.J2Onset(max(data.TrialNum(data.Run==run))==data.TrialNum & data.Run==run);
end
final_recorded_time = cumsum(final_recorded_time);


%% Pull variables from dataset (maybe do this in parent script)
%Infusion (1) no infusion (-1)
inf_no_inf = double(cellfun(@(x) strcmp(x,'A'), data.Infusion) | cellfun(@(x) strcmp(x,'B'), data.Infusion));
inf_no_inf_contrast = inf_no_inf;
inf_no_inf_contrast(inf_no_inf==0) = -1;
%inf_no_inf(inf_no_inf==0) = -1;

%Feedback (1) no feedback (-1) is signal vs baseline
feedback = double(cellfun(@(x) strcmpi(x,'Signal'), data.Feedback));
feedback_contrast = feedback;
feedback_contrast(feedback==0) = -1;

%Will improve responses yes (1) no (-1)
willImpResp = double(cellfun(@(x) strcmpi(x,'yes'), data.WillImpRespText));
willImpResp(willImpResp==0) = -1;

%Will improve responses yes (1) no (-1)
impResp = double(cellfun(@(x) strcmpi(x,'yes'), data.ImprovedRespText));
impResp(impResp==0) = -1;

%% Trial
write3Ddeconv_startTimes(reg_file_dest,data.InfOnset,data.J2Onset,sprintf('%s_trial',id),true(length(data.TrialNum),1),data);

%% Motor -- odd runs display yes/no even runs display no/yes
willImpLeft = data.WillImpRespNum==7;
%willImpLeft=double(willImpLeft); %This is for contrasts
%willImpLeft(willImpLeft==0)=-1; %This is for contrasts
impLeft = data.ImprovedRespNum==7;
%impLeft=double(impLeft); %This is for contrasts
%impLeft(impLeft==0)=-1; %This is for contrasts

%Initialize full motor arrays
motor_left = zeros(length(willImpLeft)+length(impLeft),1);
motor_onset = zeros(length(willImpLeft)+length(impLeft),1);
motor_offset = zeros(length(willImpLeft)+length(impLeft),1);

%Set 1 as left -1 as right
motor_left(1:2:end)=willImpLeft;
motor_left(2:2:end)=impLeft;
%motor_left(motor_left==0)=-1; 
%motor_right = motor_left.*-1; %This is for contrasts
motor_right = motor_left==0;
motor_right = double(motor_right);

%Set event timings
motor_onset(1:2:end)=data.WillImpOnset; %This is for contrasts
motor_onset(2:2:end)=data.ImprovedOnset; %This is for contrasts
motor_offset(1:2:end)=data.WillImpOnset+data.WillImpRt; %This is for contrasts
motor_offset(2:2:end)=data.ImprovedOnset+data.ImprovedRt; %This is for contrasts

 %Create new binary regressors
left_only = motor_left;
left_only(left_only==0)=nan;

right_only = motor_right;
right_only(right_only==0)=nan;

%Override trial index
motor_trial_index = [1 65 129 193];

write3Ddeconv_startTimes(reg_file_dest,motor_onset,motor_offset,sprintf('%s_left_only_all',id),left_only,data,motor_trial_index);
write3Ddeconv_startTimes(reg_file_dest,motor_onset,motor_offset,sprintf('%s_right_only_all',id),right_only,data,motor_trial_index);


%% Create (if any) vba regressors
try
    write3Ddeconv_startTimes(reg_file_dest,data.InfOnset,data.InfOnset+1,sprintf('%s_value_per_stim_stick',id),data.value_per_stim,data);
    write3Ddeconv_startTimes(reg_file_dest,data.Feed2Onset,data.Feed2Onset+1,sprintf('%s_PE_stick',id),data.PE,data);
catch
    if ~exist('T.PE', 'var')
        fprintf('Subject %s does not have VBA data to create regressors with!\n\n',id)
    else
        fprintf('Subject %s VBA regressors were no made successfully, please DEBUG further\n\n',id)    
    end
end

%% Infusion 1 sec stick
write3Ddeconv_startTimes(reg_file_dest,data.InfOnset,data.InfOnset+1,sprintf('%s_infusion_onset',id),true(length(data.TrialNum),1),data);

%Inf on to offset (4 sec)
write3Ddeconv_startTimes(reg_file_dest,data.InfOnset,data.WillImpOnset,sprintf('%s_infusion_onset_exact',id),true(length(data.TrialNum),1),data);

% Infusion aligned from onset to end of willImp RT
write3Ddeconv_startTimes(reg_file_dest,data.InfOnset,data.WillImpOnset+data.WillImpRt,sprintf('%s_infusion_onsetToRT',id),true(length(data.TrialNum),1),data);


%% Will Improve Responses EVENT TIMES ONLY - RT convolv
write3Ddeconv_startTimes(reg_file_dest,data.WillImpOnset,data.WillImpOnset+data.WillImpRt,sprintf('%s_willImpRespEvent',id),true(length(data.TrialNum),1),data);
write3Ddeconv_startTimes(reg_file_dest,data.WillImpOnset,data.WillImpOnset+1,sprintf('%s_willImpRespEventStick',id),true(length(data.TrialNum),1),data);

%% Will Improve Responses ACTUAL RESPONSES AS PARAMETERIC MOD
write3Ddeconv_startTimes(reg_file_dest,data.WillImpOnset,data.WillImpOnset+data.WillImpRt,sprintf('%s_willImpRespMod',id),willImpResp,data);



%IR IRo
%write3Ddeconv_startTimes(reg_file_dest,data.WillImpOnset,data.J1Onset,sprintf('%s_IR',id),inf_no_inf.*willImpResp,data);
%write3Ddeconv_startTimes(reg_file_dest,data.WillImpOnset,data.J1Onset,sprintf('%s_IRo',id),~inf_no_inf.*willImpResp,data);

write3Ddeconv_startTimes(reg_file_dest,data.WillImpOnset,data.J1Onset,sprintf('%s_IRp',id),inf_no_inf.*data.WillImpRespBin,data);
write3Ddeconv_startTimes(reg_file_dest,data.WillImpOnset,data.J1Onset,sprintf('%s_IRp0',id),~inf_no_inf.*data.WillImpRespBin,data);

%% Infusion / no infusion 
write3Ddeconv_startTimes(reg_file_dest,data.InfOnset,data.WillImpOnset,sprintf('%s_infusion_no_infusion_exact',id),inf_no_inf_contrast,data);
write3Ddeconv_startTimes(reg_file_dest,data.InfOnset,data.WillImpOnset,sprintf('%s_infusion_exact',id),inf_no_inf,data);
write3Ddeconv_startTimes(reg_file_dest,data.InfOnset,data.WillImpOnset,sprintf('%s_no_infusion_exact',id),~inf_no_inf,data);

%Align with infusion time (4sec)
%write3Ddeconv_startTimes(reg_file_dest,data.InfOnset,data.InfOnset+4,sprintf('%s_infusion_no_infusion_',id),inf_no_inf,data);

%Align inf no inf with RT
write3Ddeconv_startTimes(reg_file_dest,data.WillImpOnset,data.WillImpOnset+data.WillImpRt,sprintf('%s_infusion_no_infusion_RTconvolv',id),inf_no_inf_contrast,data);

%Align inf no inf from onset to end of RT
write3Ddeconv_startTimes(reg_file_dest,data.InfOnset,data.WillImpOnset+data.WillImpRt,sprintf('%s_infusion_no_infusion_onsetToRT',id),inf_no_inf_contrast,data);

%% Feedback 1 sec stick - start at inflection of feedback signal
%Feed 1 onset is the start of the feedback screen
%Feed 2 onset is the start of inflection (1.5-3.3 secs)
%Feed 3 osnet is the point when the inflection reaches its peak

%Feedback 1 second stick
write3Ddeconv_startTimes(reg_file_dest,data.Feed2Onset,data.Feed2Onset+1,sprintf('%s_feedback2_onset',id),true(length(data.TrialNum),1),data);

%Align feedback from feedback onset to end of RT
write3Ddeconv_startTimes(reg_file_dest,data.Feed2Onset,data.ImprovedOnset+data.ImprovedRt,sprintf('%s_feedback2_onsetToRT',id),true(length(data.TrialNum),1),data);

%Align with exact time immp onset
write3Ddeconv_startTimes(reg_file_dest,data.Feed2Onset,data.ImprovedOnset,sprintf('%s_feedback2_onset_exact',id),true(length(data.TrialNum),1),data);


%write3Ddeconv_startTimes(reg_file_dest,data.Feed1Onset,data.ImprovedOnset,sprintf('%s_feedback_all',id),true(length(data.TrialNum),1),data);
%%Delete? ^^^

%% Feedback / no Feedback -- determine which feedback gives the best signal -- it is fb2
write3Ddeconv_startTimes(reg_file_dest,data.Feed2Onset,data.Feed2Onset+1,sprintf('%s_feedback_no_feedback',id),feedback,data);

%Align as exact time difference 
write3Ddeconv_startTimes(reg_file_dest,data.Feed2Onset,data.ImprovedOnset,sprintf('%s_feedback_no_feedback_exact',id),feedback_contrast,data);
write3Ddeconv_startTimes(reg_file_dest,data.Feed2Onset,data.ImprovedOnset,sprintf('%s_feedback_exact',id),feedback,data);
write3Ddeconv_startTimes(reg_file_dest,data.Feed2Onset,data.ImprovedOnset,sprintf('%s_no_feedback_exact',id),~feedback,data);

%Align feedback no feedback with RT
write3Ddeconv_startTimes(reg_file_dest,data.ImprovedOnset,data.ImprovedOnset+data.ImprovedRt,sprintf('%s_feedback_no_feedback_RTconvolv',id),feedback,data);

%Align feedback no feedback from feedback onset to the end of RT
write3Ddeconv_startTimes(reg_file_dest,data.Feed2Onset,data.ImprovedOnset+data.ImprovedRt,sprintf('%s_feedback_no_feedback_onsetToRT',id),feedback,data);


%% Improve Responses EVENT TIMES ONLY - RT convolv
write3Ddeconv_startTimes(reg_file_dest,data.ImprovedOnset,data.ImprovedOnset+data.ImprovedRt,sprintf('%s_impRespEvent',id),true(length(data.TrialNum),1),data);
write3Ddeconv_startTimes(reg_file_dest,data.ImprovedOnset,data.ImprovedOnset+1,sprintf('%s_impRespEventStick',id),true(length(data.TrialNum),1),data);

%% Improve Responses ACTUAL RESPONSES AS PARAMETERIC MOD
write3Ddeconv_startTimes(reg_file_dest,data.ImprovedOnset,data.ImprovedOnset+data.ImprovedRt,sprintf('%s_impRespMod',id),impResp,data);


%write3Ddeconv_startTimes(reg_file_dest,data.ImprovedOnset,data.ImprovedOnset+data.J2Onset,sprintf('%s_FR',id),feedback.*impResp,data);
%write3Ddeconv_startTimes(reg_file_dest,data.ImprovedOnset,data.ImprovedOnset+data.J2Onset,sprintf('%s_FRo',id),~feedback.*impResp,data);

write3Ddeconv_startTimes(reg_file_dest,data.ImprovedOnset,data.J2Onset,sprintf('%s_FRp',id),feedback.*data.ImprovedRespBin,data);
write3Ddeconv_startTimes(reg_file_dest,data.ImprovedOnset,data.J2Onset,sprintf('%s_FRp0',id),~feedback.*data.ImprovedRespBin,data);

%% Parametric interactions
write3Ddeconv_startTimes(reg_file_dest,data.ImprovedOnset,data.ImprovedOnset+data.ImprovedRt,sprintf('%s_feedbackXinf_no_inf',id),feedback.*inf_no_inf,data);
write3Ddeconv_startTimes(reg_file_dest,data.Feed2Onset,data.ImprovedOnset+data.ImprovedRt,sprintf('%s_feedbackXinf_no_inf_fb_onset',id),feedback.*inf_no_inf,data);


%% Censor
%Only sonrisa 1 right nbow
%if strfind(id,'SON1')
try
    create_censor_file_2_1(data,id)
catch
    fprintf('\n%s could not construct Censor file! Check proc dir...\n', num2str(id))
end
%end


function [x,y]=write3Ddeconv_startTimes(file_loc,event_beg,event_end,fname,modulator,data,override)
% Function will write FSL styled regressors in dat files for fMRI analysis
% Inputs:
% file_loc: file location (str)
% event_beg: the time in miliseconds of the event beginning
% event_end: the time in milliseconds of the event ending
% fname: the file names
% censor: the censor vector or parametric value vector depending on the regressor
% trial_index: the position of when a new block starts (trialwise)

global gen_concat_reg final_recorded_time
 
format long

%Get the number of runs
n_runs = max(data.Run);

%If there is no trial override given
try
    override;
catch
    override=[];
end

%args should be cols
if ~iscolumn(event_beg)
    event_beg=event_beg';
end

if ~iscolumn(event_end)
    event_end=event_end';
end

x(:,1) = event_beg;
x(:,2) = event_end-event_beg;
%x=x./1000; %Convert to seconds (not for clock already in seconds)
x(:,3) = ones(length(x),1).*modulator; %originally was modulator'
%write the -stim_times_FSL


%Save to regs folder
%As we transition to FSL make this an options of we want to write single
%text files for each regressor for each run or one concatenated regressor
%for AFNI

if gen_concat_reg
    c = asterisk(x,data,override); %Add in asterisks and clean up data
    dlmcell([file_loc fname '.dat'],c,'delimiter','\t')
else
    
    if ~exist([file_loc '/fsl_regs/'],'dir')
        mkdir([file_loc '/fsl_regs/'])
    end
    
    tmp_run_data = [];
    
    for i = 1:n_runs
        run_data.(['run_' num2str(i)]) = x(i==data.Run,:);
        nan_data = sum(isnan(run_data.(['run_' num2str(i)])),2)>=1;
        run_data.(['run_' num2str(i)])(nan_data,:)=[];
        run_data.(['run_' num2str(i)])=array2table(run_data.(['run_' num2str(i)]));
        if i > 1
           run_data.(['run_' num2str(i)])(:,1).Var1 = run_data.(['run_' num2str(i)])(:,1).Var1 + final_recorded_time(i-1);
        end
        tmp_run_data = [tmp_run_data; run_data.(['run_' num2str(i)])];
        %writetable(run_data,[file_loc '/fsl_regs/' fname '_run_' num2str(i) '.dat'],'Delimiter','\t','WriteVariableNames',0)
    end
end

writetable(tmp_run_data,[file_loc '/fsl_regs_helmet/' fname '.dat'],'Delimiter','\t','WriteVariableNames',0)

y=0;

return

function c = asterisk(x,b,trial_index)
%adding asterisk to existing .dat files also removes any nans present

c=[];
ast = {'*', '*', '*'};

%For now manually override trial index if you want to
if isempty(trial_index)
    trial_index=find(b.TrialNum==1);
end

%Short while loop to place asterisks for regressors
while ~isempty(trial_index)
%     if sum(isnan(b.Jitter1(trial_index(1):trial_index(1)+31,:)))==length(trial_index(1):trial_index(1)+31)
%         %pass if no data was recorded
    if 1<length(trial_index)
        block_data = num2cell(x(trial_index(1): trial_index(2)-1,:));
        c = [c; block_data; ast];
    else
        block_data = num2cell(x(trial_index(1): end,:));
        c = [c; block_data;];
    end
    trial_index(1)=[];
end

%DO NOTE SON1_008_a has duplicate * * *'s due to missing one block, do we
%need to account for this? or will afni figure it out?


%clean up any nans
%fh = @(y) all(isnan(y(:)));
c = c(~any(cellfun(@isnan,c),2),:);
%c(cellfun(fh, c)) = [];
%Check on c!

return

function [reg_out_onset, reg_out_offset]=create_motor_reg(data,motor_bin)

%Initialize
motor_onset = zeros(length(data.WillImpOnset)+length(data.ImprovedOnset),1);
motor_offset = motor_onset;

motor_onset(1:2:end)=data.WillImpOnset; 
motor_onset(2:2:end)=data.ImprovedOnset; 
motor_offset(1:2:end)=data.WillImpOnset+data.WillImpRt; 
motor_offset(2:2:end)=data.ImprovedOnset+data.ImprovedRt; 

%Just in case it's not a logical
motor_bin = logical(motor_bin);

%Onsets
motor_onset = motor_onset(motor_bin);
reg_out_onset = zeros(length(motor_onset_willImp)+length(motor_onset_imp),1);
reg_out_onset(1:2:end)=motor_onset_willImp; 
reg_out_onset(2:2:end)=motor_onset_imp; 

%Offsets
motor_offset_willImp = data.WillImpOnset(motor_bin) + data.WillImpRt(motor_bin);
motor_offset_imp = data.ImprovedOnset(motor_bin) + data.ImprovedRt(motor_bin);
reg_out_offset = zeros(length(motor_offset_willImp)+length(motor_offset_imp),1);
reg_out_offset(1:2:end)=motor_onset_willImp; 
reg_out_offset(2:2:end)=motor_onset_imp; 


        
