function out=compile_nfb_data_1(regressor_flag)
%Funciton will organize and create simple plots of nfb task data. For more
%information see README at https://github.com/heffjos/nfb.

%%Create regreossors or not
try regressor_flag; catch, regressor_flag=0; end

%% Datalocation -- set as needed
if ispc
    data_path = 'E:\Box Sync\fMRI_Shared\NFB\NFB_response\'; %Change as needed
    %For the time being because depending on the system we might not have a
    %'/' in front of the string
    str_direction=1;
else
    data_path = '/Users/martapecina/Box Sync/PITT/RESEARCH/fMRI_shared/NFB/NFB_response/';
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
sortedNames = sort(son1_all_tmp.Properties.VariableNames(2:end));
out.son1_all = [son1_all_tmp(:,1) son1_all_tmp(:,sortedNames)];
out.son2_all = [son2_all_tmp(:,1) son2_all_tmp(:,sortedNames)];

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


%Save new file in data location
writetable(T,[data_dir{:} sprintf('subj_%s_all_runs.csv',id{:})])

%Create regressors if needed -- this could be moved out of here if need be
if options.rflag
    write_nfb_regressors_2(out.(proto{:}).(['subj' id{:}]).(['admin' admin_code]),out.(proto{:}).(['subj' id{:}]).(['admin' admin_code '_reg_id']))
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
        
        %Add to SON1 admin 1 dataframe
        data.SON1_MDF_1 = [data.SON1_MDF_1; tmp_data]; 
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
        
        
