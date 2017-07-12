function plot_nfb_data
%Function to make a bunch of plots for the nfb data

%Clear up any open graphs
close all

%Load in main data structure  (out)
%load('/Users/martapecina/Google Drive/fMRI_Responses/NFB_response/nfball.mat')
%load('E:\Box Sync\fMRI_shared\NFB_response\nfball.mat')
load('/Users/martapecina/Box Sync/PITT/RESEARCH/fMRI_shared/NFB/NFB_response/nfball.mat')
protocols = {'SON1', 'SON2'};
fig_ct=1;
for protocol = protocols
    %Plot group level data here
    subjs=fieldnames(out.(protocol{:}))';
    son1_and_2 = [out.SON1_MDF_1;out.SON1_MDF_2];
    group_specs=[];
    group_specs=compile_group_data(son1_and_2,group_specs);
    fig_ct=plot_group_average_response(son1_and_2,fig_ct,group_specs);
    for subj = subjs
        admins=fieldnames(out.(protocol{:}).(subj{:}))';
        for admin = admins
            %Pull single subject data
            data=out.(protocol{:}).(subj{:}).(admin{:});
            
            %Get subject name
            subj_specs.name=data.Participant{1};
            subj_specs.current_admin=admin{:};
            
            %Pull all plotting data from single subject
            subj_specs=compile_single_subject_data(data,subj_specs);
            
            %Set up simple plots
            %fig_ct=plot_subj_responses_per_trial(fig_ct,subj_specs);
            %fig_ct=gen_rt_histograms(data,fig_ct,subj_specs);
            %fig_ct=plot_subj_responses_by_infusion_no_infusion(data,fig_ct,subj_specs);
            %fig_ct=plot_learning_curves(fig_ct,subj_specs);
            fig_ct=plot_one_graph_resp_per_condition(data,fig_ct,subj_specs);
            
        end
    end
end


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
function fig_num=plot_subj_responses_per_trial(fig_num,subj_specs)
figure(fig_num)
subplot(2,1,1)
%yes_vs_no=cellfun(@(x) str2double(x),cellfun(@count_subj_response, data.WillImpRespText, 'UniformOutput', false));
plot(subj_specs.will_imp)
title(sprintf('Subject %s will improve responses per trial %s',subj_specs.name,subj_specs.current_admin))
axis([0 length(subj_specs.will_imp) -.5 1.5])
subplot(2,1,2)
%yes_vs_no=cellfun(@(x) str2double(x),cellfun(@count_subj_response, data.ImprovedRespText, 'UniformOutput', false));
plot(subj_specs.imp, 'r')
axis([0 length(subj_specs.imp) -.5 1.5])
title(sprintf('Subject %s improved responses per trial %s',subj_specs.name,subj_specs.current_admin))
fig_num = fig_num+1;

%--------------------------------------------------------------------------
function fig_num=gen_rt_histograms(data,fig_num,subj_specs)
figure(fig_num)
subplot(2,1,1)
histogram(data.WillImpRt)
title(sprintf('Subject %s RTs of %s',subj_specs.name,subj_specs.current_admin))
subplot(2,1,2)
histogram(data.ImprovedRt)
title(sprintf('Subject %s RTs of %s',subj_specs.name,subj_specs.current_admin))
fig_num = fig_num+1;

%--------------------------------------------------------------------------
function fig_num=plot_subj_responses_by_infusion_no_infusion(data,fig_num,subj_specs)
conditions = unique(data.Infusion);
subplot_idx=[1:2; 3:4; 5:6; 7:8];

for i = 1:length(conditions)
    cond_idx_infusion = cellfun(@(x) strcmp(conditions{i},x), data.Infusion) & cellfun(@(x) strcmp('Signal',x), data.Feedback);
    cond_idx_baseline = cellfun(@(x) strcmp(conditions{i},x), data.Infusion) & cellfun(@(x) strcmp('Baseline',x), data.Feedback);
    
    inf_response_willImp = cellfun(@(x) str2double(x),cellfun(@count_subj_response, data.WillImpRespText(cond_idx_infusion), 'UniformOutput', false));
    baseline_response_willImp = cellfun(@(x) str2double(x),cellfun(@count_subj_response, data.WillImpRespText(cond_idx_baseline), 'UniformOutput', false));
    
    inf_response_Imp = cellfun(@(x) str2double(x),cellfun(@count_subj_response, data.ImprovedRespText(cond_idx_infusion), 'UniformOutput', false));
    baseline_response_Imp = cellfun(@(x) str2double(x),cellfun(@count_subj_response, data.ImprovedRespText(cond_idx_baseline), 'UniformOutput', false));
    
    %Create the figure
    figure(fig_num)
    subplot(4,2,subplot_idx(i,1))
    %%Signal%%
    plot(inf_response_willImp,'o','MarkerSize',15,'MarkerFaceColor',[0,0,1]) %Add 2 to seperate the two responses [2 - 3]
    hold on
    plot(inf_response_Imp,'o','MarkerFaceColor',[1,0,0])
    ylim([-.5 1.5])
    
    ylabel(sprintf('Condition %s',conditions{i}),'FontSize',12)
    if strcmp(conditions{i}, 'A')
        title('Signal')
        legend('Will Imp. Resp.', 'Imp. Resp.','Location','best')
    end
    
    %%Baseline%%
    subplot(4,2,subplot_idx(i,2))
    plot(baseline_response_willImp,'o','MarkerSize',15,'MarkerFaceColor',[0,0,1]) %Add 2 to seperate the two responses [2 - 3]
    hold on
    plot(baseline_response_Imp,'o','MarkerFaceColor',[1,0,0])
    ylim([-.5 1.5])
    
    if strcmp(conditions{i}, 'A')
        title('Baseline')
        legend('Will Imp. Resp.', 'Imp. Resp.','Location','best')
    end
end

fig_num = fig_num+1;

savefig(sprintf('Subj_%s_Inf_vs_noInf_per_condition',subj_specs.name))

%--------------------------------------------------------------------------
function fig_num=plot_learning_curves(fig_num,subj_specs)
%Plot "learning curves" ie will improve (1 vs 0) and condition A or B
%present and improve (1 vs 0) if signal or baseline. All conditions
%included

figure(fig_num)
%Plot will improve answers based on if signal was to be expected
subplot(2,1,1)
plot(smooth(subj_specs.signal_a_or_b),'r--')
hold on
plot(smooth(subj_specs.will_imp))
title(sprintf('Subject %s Protocol/Calibration vs Will Imp. Resp %s',subj_specs.name, subj_specs.current_admin))
legend('Contingency', 'Response','Location','best')
%Plot improved answers based on if signal or baseline
subplot(2,1,2)
plot(smooth(subj_specs.signal_or_baseline),'r--')
hold on
plot(smooth(subj_specs.imp))
title(sprintf('Subject %s Inf/noInf vs Imp. Resp %s',subj_specs.name, subj_specs.current_admin))

%savefig(sprintf('Subj_%s_learning_curves',subj_specs.name))

fig_num = fig_num+1;

%--------------------------------------------------------------------------
function subj_specs=compile_single_subject_data(data,subj_specs)
%Protocol/Calibration
subj_specs.signal_a_or_b = double(cellfun(@(x) strcmp(x,'A'), data.Infusion) | cellfun(@(x) strcmp(x,'B'), data.Infusion));

%Inf/no-inf
subj_specs.signal_or_baseline = double(cellfun(@(x) strcmp(x,'Signal'), data.Feedback));

%Subj responses for will improve
subj_specs.will_imp=cellfun(@(x) str2double(x),cellfun(@count_subj_response, data.WillImpRespText, 'UniformOutput', false));

%Subj responses for improved
subj_specs.imp=cellfun(@(x) str2double(x),cellfun(@count_subj_response, data.ImprovedRespText, 'UniformOutput', false));

subj_specs.cond_a_index = find(cellfun(@(x) strcmp('A',x), data.Infusion));
subj_specs.cond_b_index = find(cellfun(@(x) strcmp('B',x), data.Infusion));
subj_specs.cond_c_index = find(cellfun(@(x) strcmp('C',x), data.Infusion));
subj_specs.cond_d_index = find(cellfun(@(x) strcmp('D',x), data.Infusion));

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

group_specs.name = 'Son1 admin 1 and 2';
group_specs.current_admin='NA';

%--------------------------------------------------------------------------
function fig_num=plot_one_graph_resp_per_condition(data,fig_num,subj_specs)
conditions = unique(data.Infusion);
% colors = {[0,0,1], [0,0,.9];
%           [1,0,0], [.9,0,0];
%           [0,1,0], [0,.9,0];
%           [1,0,1], [.7,0,.7];};
colors = {'r','b','g','k'};
subplot_idx=[1:2; 3:4; 5:6; 7:8];
for i = 1:length(conditions)
    %Initialize for nans
    plot_data_will_imp = nan(length(subj_specs.will_imp),1);
    plot_data_imp = nan(length(subj_specs.imp),1);
    cond_index=find(cellfun(@(x) strcmp(conditions{i},x), data.Infusion));
    
    %Set the data
    plot_data_will_imp = subj_specs.will_imp(cond_index);
    %plot_data_will_imp(cond_index) = subj_specs.will_imp(cond_index);
    plot_data_imp = subj_specs.imp(cond_index);
    %plot_data_imp(cond_index) = subj_specs.imp(cond_index);
    
    
    %%JUST OVER LAY THE CONTIONGENCY NOW!! to have all condition on one
    %%plot use subplot(2,1,1) and comment out the contingencies
    
    figure(fig_num)
    subplot(2,1,1)
    plot(smooth(plot_data_will_imp))
    ylim([-.5 1.5])
    hold on
    title(sprintf('Subject %s Will improve by stimulus %s',subj_specs.name, subj_specs.current_admin))
    
    subplot(2,1,2)
    plot(smooth(plot_data_imp))
    ylim([-.5 1.5])
    hold on
    title(sprintf('Subject %s Improved by stimulus %s',subj_specs.name, subj_specs.current_admin))
    
    
    %subplot(2,1,1)%     subplot(4,2,subplot_idx(i,1))
    %     hold on
    %     plot(plot_data_will_imp,colors{i})
    %     title(sprintf('Subject %s Protocol/Calibration vs Will Imp. Resp per condition %s',subj_specs.name, subj_specs.current_admin))
    %     %subplot(2,1,2)
    %     subplot(4,2,subplot_idx(i,2))
    %     hold on
    %     plot(plot_data_imp,colors{i})
    %     title(sprintf('Subject %s Inf/noInf vs Imp. Resp per condition %s',subj_specs.name, subj_specs.current_admin))
    
    
    % %     plot(plot_data_will_imp,'o','MarkerSize',15,'MarkerFaceColor',colors{i,1}) %Add 2 to seperate the two responses [2 - 3]
    % %     hold on
    % %     plot(plot_data_imp,'o','MarkerFaceColor',colors{i,2})
    % %     ylim([-.5 1.5])
end
%savefig(sprintf('Subj_%s_graph_resp_per_condition%s_admin',subj_specs.name,subj_specs.current_admin))

legend('A', 'B', 'C', 'D', 'Location','best')
fig_num = fig_num+1;


function fig_num = plot_group_average_response(data,fig_num,group_specs)
conditions = unique(data.Infusion);
for i = 1:length(conditions)
    figure(fig_num)
    subplot(2,1,1)
    plot(smooth(group_specs.(['cond_' lower(conditions{i}) '_willImpMeanResponse'])))
    hold on
    title(sprintf('Subject %s Will improve by stimulus %s',group_specs.name, group_specs.current_admin))
    
    subplot(2,1,2)
    plot(smooth(group_specs.(['cond_' lower(conditions{i}) '_ImpMeanResponse'])))
    hold on
    title(sprintf('Subject %s Improved by stimulus %s',group_specs.name, group_specs.current_admin))
end
legend('A', 'B', 'C', 'D', 'Location','best')
fig_num = fig_num+1;