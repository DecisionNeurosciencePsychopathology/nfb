function out=compile_nfb_data(regressor_flag)
%Funciton will organize and create simple plots of nfb task data. For more
%information see README at https://github.com/heffjos/nfb.

%%Create regreossors or not
try regressor_flag; catch, regressor_flag=0; end

%% Datalocation -- set as needed
if ispc
    data_path = 'E:\Box Sync\fMRI_Shared\NFB\NFB_response\'; %Change as needed
else
    data_path = '/Users/martapecina/Box Sync/PITT/RESEARCH/fMRI_shared/NFB/NFB_response/';
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

%% Compile loop

%Grab parent directories
experiment_paths=glob([data_path 'SON*'])';

for experiment_path = experiment_paths
    
    if strfind(experiment_path{:},'SONRISA1')
        experiment_path = {[experiment_path{:} 'CorrectResponse']};
    end
    
    %Grab dir list from datalocation
    expression = '.*\d{3,9}.*';
    data_dirs = regexp(glob([experiment_path{:} '/SON*']), expression, 'match'); %Think of a way to just get the uniue ids here , probably with cellfun
    data_dirs = [data_dirs{:}];
    
    for data_dir = data_dirs
        %Put everything into one table
        out=compile_single_subject_data(out,data_dir,options);
    end
end

%% Save final output
save([data_path '/nfball'],'out')


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
admin_code=data_dir{:}(admin_code{:}-1);

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

%Parse out the data for easier regressor creation
%out.(proto{:}).(['subj' id{:}]).(['admin' admin_code '_data'])=parse_single_subj_data(out.(proto{:}).(['subj' id{:}]).(['admin' admin_code]),proto{:},id,reg_suffix);
tmp_data=parse_single_subj_data(out.(proto{:}).(['subj' id{:}]).(['admin' admin_code]),proto{:},id,reg_suffix);
out.(proto{:}).(['subj' id{:}]).(['admin' admin_code '_reg_id']) = tmp_data.reg_id;

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

%Create regressors if needed
if options.rflag
    write_nfb_regressors(out.(proto{:}).(['subj' id{:}]).(['admin' admin_code]),out.(proto{:}).(['subj' id{:}]).(['admin' admin_code '_reg_id']))
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



