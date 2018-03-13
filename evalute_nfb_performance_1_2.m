function evalute_nfb_performance_1_2
%Function will evalute nfb subjects responses, place resutls into graphics
%for ease of viewing and assisting with evaluation of subejcts perofmrance 
%on a case by case basis as well as experiment wide performance.

%Load i nthe data automatically
if ispc
    load('E:\Box Sync\fMRI_shared\NFB\NFB_response\SON1&2_behav_results\nfball.mat')
else
    load('/Users/martapecina/Box Sync/PITT/RESEARCH/fMRI_shared/NFB/NFB_response/SON1&2_behav_results/nfball.mat')
end


%Grab field names -- we only want to use Master Data Files
fnames = fieldnames(out);

for i = 1:length(fnames)
   if strfind(fnames{i},'MDF')>0 
       
       %Set up options
       options.average_response_plot_fig_num = 10+i;
       options.dfname = fnames{i};
       
       %MDF data
       data=out.(fnames{i});

       %Expermient wide plot
       options.bar_plot_fig_number = i+100;
       plot_behavior(data,options);
       
       %TODO
       %Grab the indivual subjects here (regrep(Participant numbers, unqiue
       %it?, cellfun(@(x) strcmp...) feed to
       %bar_plot_behavior(data,graph_options)
       
       %Individual subject loop
%        subjs=regexp(out.(fnames{i}).Participant,'\d{3}','match');
%        subjs=unique([subjs{:}]');
%        for j = 1:length(subjs)
%            subjs_data_idx = cellfun(@(x) ~isempty(strfind(x,subjs{j})),data.Participant,'UniformOutput',0);
%            subjs_data = data(cell2mat(subjs_data_idx),:);
%            plot_behavior(subjs_data,options);
%        end
       
   end
end


function plot_behavior(data,options)
%TODO
%Make this the generic plot funciton with a second arguement graph_options
%graph_options will have things like Infusion_title, subject(s)_name,
%not_infusion_title, type_of_plot(bar line ... )ect...
infusion_trials = cellfun(@(x) strcmpi(x,'A') | strcmpi(x,'B') ,data.Infusion);
will_imp_resp_inf = data.WillImpRespBin(infusion_trials);
will_imp_resp_no_inf = data.WillImpRespBin(~infusion_trials);
figure(options.bar_plot_fig_number)
subplot(1,2,1)
bar([sum(will_imp_resp_inf==1),sum(will_imp_resp_inf==0), sum(isnan(will_imp_resp_inf))])
title([options.dfname ' Infusion will improve'])
ax=gca;
ax.XTickLabel = {'Yes', 'No', 'Nan'};
subplot(1,2,2)
bar([sum(will_imp_resp_no_inf==1),sum(will_imp_resp_no_inf==0), sum(isnan(will_imp_resp_no_inf))])
title('No Infusion will improve')
ax=gca;
ax.XTickLabel = {'Yes', 'No', 'Nan'};

imp_resp_inf = data.ImprovedRespBin(infusion_trials);
imp_resp_no_inf = data.ImprovedRespBin(~infusion_trials);
figure(options.bar_plot_fig_number*2)
subplot(1,2,1)
bar([sum(imp_resp_inf==1),sum(imp_resp_inf==0), sum(isnan(imp_resp_inf))])
title([options.dfname ' Infusion improve'])
ax=gca;
ax.XTickLabel = {'Yes', 'No', 'Nan'};
subplot(1,2,2)
bar([sum(imp_resp_no_inf==1),sum(imp_resp_no_inf==0), sum(isnan(imp_resp_no_inf))])
title('No Infusion improve')
ax=gca;
ax.XTickLabel = {'Yes', 'No', 'Nan'};


%Plot the mean responses for son 1 (with son2 plac subjects) & 2 and son 2
%Plac and Nalt
group_specs=[];
group_specs=compile_group_data(data,group_specs);
group_specs.name = options.dfname;
plot_group_average_response(data,options.average_response_plot_fig_num,group_specs)






%--------------------------------------------------------------------------
function group_specs=compile_group_data(data,group_specs)
%Inf/no-inf
group_specs.signal_a_or_b = double(cellfun(@(x) strcmp(x,'A'), data.Infusion) | cellfun(@(x) strcmp(x,'B'), data.Infusion));

%Feedback
group_specs.signal_or_baseline = double(cellfun(@(x) strcmp(x,'Signal'), data.Feedback));

%Subj responses for will improve
group_specs.will_imp=cellfun(@(x) str2double(x),cellfun(@count_subj_response, data.WillImpRespText, 'UniformOutput', false));

%Subj responses for improved
group_specs.imp=cellfun(@(x) str2double(x),cellfun(@count_subj_response, data.ImprovedRespText, 'UniformOutput', false));

conditions = {'a','b','c','d'};

%%we will use version as out key instead of subject becasue of naming
%%convention
%This works becasue there are only 4 versions! not valid if this changes!
versions = unique(data.Version)';
for version = versions
    ver = data.Version==version;
    for condition = conditions
        %Find the indices of the condition per version
        group_specs.(['cond_' condition{:} '_index']) = find(cellfun(@(x) strcmpi(condition{:},x), data.Infusion) & ver);
        
        %Find the subj responses
        tmp_willImp = group_specs.will_imp(group_specs.(['cond_' condition{:} '_index']));
        tmp_imp = group_specs.imp(group_specs.(['cond_' condition{:} '_index']));
        
        %Reshape
        cols=length(tmp_willImp)/max(data.TrialNum);
        group_specs.willImpResp.(condition{:}).(['cond_' condition{:} '_willImpResponse_' num2str(version)]) = reshape(tmp_willImp,max(data.TrialNum),cols);
        group_specs.ImpResp.(condition{:}).(['cond_' condition{:} '_ImpResponse_' num2str(version)]) = reshape(tmp_imp,max(data.TrialNum),cols);
    end
end

%Merge the structs into one array
for condition  = conditions
    tmp_willImp=[];
    tmp_imp=[];
    fnames_will=fieldnames(group_specs.willImpResp.(condition{:}))';
    fnames_imp=fieldnames(group_specs.ImpResp.(condition{:}))';
    for fname = fnames_will
        tmp_willImp=[tmp_willImp group_specs.willImpResp.(condition{:}).(fname{:})];
    end
    
    for fname = fnames_imp
        tmp_imp=[tmp_imp group_specs.ImpResp.(condition{:}).(fname{:})];
    end
    
    %Take the mean
    group_specs.(['cond_' condition{:} '_willImpMeanResponse']) = nanmean(tmp_willImp,2)';
    group_specs.(['cond_' condition{:} '_ImpMeanResponse']) = nanmean(tmp_imp,2)';
    
end

% group_specs.name = 'Son1';
  group_specs.current_admin='NA';


%--------------------------------------------------------------------------
function y = count_subj_response(str)
if strcmpi('yes',str)
    y = '1';
elseif strcmpi('no',str)
    y = '0';
else
    y = str;
end

%--------------------------------------------------------------------------
function fig_num = plot_group_average_response(data,fig_num,group_specs)
conditions = unique(data.Infusion);
for i = 1:length(conditions)
    figure(fig_num)
    subplot(2,1,1)
    plot(smooth(group_specs.(['cond_' lower(conditions{i}) '_willImpMeanResponse'])))
    hold on
    title(sprintf('Subject %s Will improve by stimulus %s',group_specs.name, group_specs.current_admin))
    axis([0 35 0 1])
    
    subplot(2,1,2)
    plot(smooth(group_specs.(['cond_' lower(conditions{i}) '_ImpMeanResponse'])))
    hold on
    title(sprintf('Subject %s Improved by stimulus %s',group_specs.name, group_specs.current_admin))
    axis([0 35 0 1])
end
legend('A', 'B', 'C', 'D', 'Location','best')
fig_num = fig_num+1;