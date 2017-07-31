function write_nfb_regressors(data,id)

%% Datalocation -- set as needed
if ispc
    reg_file_dest = 'C:/kod/nfb/analysis/regs/'; %Change as needed
else
    reg_file_dest = '/Users/martapecina/GitHub/nfb_analysis/regs/';
end


%% Pull variables from dataset (maybe do this in parent script)
%Infusion (1) no infusion (0)
inf_no_inf = double(cellfun(@(x) strcmp(x,'A'), data.Infusion) | cellfun(@(x) strcmp(x,'B'), data.Infusion));

%Feedback (1) no feedback (0) is signal vs baseline
signal_vs_baseline = double(cellfun(@(x) strcmpi(x,'Signal'), data.Feedback) | cellfun(@(x) strcmp(x,'Baseline'), data.Feedback));

%Will improve resposnes yes (1) no (0)
willImpResp = data.WillImpRespNum==7;

%Will improve resposnes yes (1) no (0)
impResp = data.ImprovedRespNum==7;

%% Trial
write3Ddeconv_startTimes(reg_file_dest,data.InfOnset,data.J2Onset,sprintf('%s_trial',id),true(length(data.TrialNum),1),0,data);

%% Infusion 1 sec stick
write3Ddeconv_startTimes(reg_file_dest,data.InfOnset,data.InfOnset+1,sprintf('%s_infusion_onset',id),true(length(data.TrialNum),1),0,data);

%% Infusion / no infusion 
write3Ddeconv_startTimes(reg_file_dest,data.InfOnset,data.WillImpOnset,sprintf('%s_infusion_no_infusion',id),inf_no_inf,0,data);

%% Feedback / no Feedback
write3Ddeconv_startTimes(reg_file_dest,data.Feed1Onset,data.Feed3Onset,sprintf('%s_feedback_no_feedback',id),signal_vs_baseline,0,data);

%% Will Improve Responses EVENT TIMES ONLY - RT convolv
write3Ddeconv_startTimes(reg_file_dest,data.WillImpOnset,data.WillImpOnset+data.WillImpRt,sprintf('%s_willImpRespEvent',id),true(length(data.TrialNum),1),0,data);

%% Will Improve Responses ACTUAL RESPONSES AS PARAMETERIC MOD
write3Ddeconv_startTimes(reg_file_dest,data.WillImpOnset,data.WillImpOnset+data.WillImpRt,sprintf('%s_willImpRespMod',id),willImpResp,0,data);

%% Feedback 1 sec stick
write3Ddeconv_startTimes(reg_file_dest,data.Feed1Onset,data.Feed1Onset+1,sprintf('%s_feedback_onset',id),true(length(data.TrialNum),1),0,data);

%% Improve Responses EVENT TIMES ONLY - RT convolv
write3Ddeconv_startTimes(reg_file_dest,data.ImprovedOnset,data.ImprovedOnset+data.ImprovedRt,sprintf('%s_impRespEvent',id),true(length(data.TrialNum),1),0,data);

%% Improve Responses ACTUAL RESPONSES AS PARAMETERIC MOD
write3Ddeconv_startTimes(reg_file_dest,data.ImprovedOnset,data.ImprovedOnset+data.ImprovedRt,sprintf('%s_impRespMod',id),impResp,0,data);

%% Censor
%Worry about this later, ge the model done first
%create_censor_file(data,id)


function [x,y]=write3Ddeconv_startTimes(file_loc,event_beg,event_end,fname,modulator,noFSL,b)
% Function will write FSL styled regressors in dat files for fMRI analysis
% Inputs:
% file_loc: file location (str)
% event_beg: the time in miliseconds of the event beginning
% event_end: the time in milliseconds of the event ending
% fname: the file names
% censor: the censor vector or parametric value vector depending on the regressor
% noFSL either to write a FSL file or a different single line version (see 3dDeconvolve help for more info)
% trial_index: the position of when a new block starts (trialwise)

if nargin <6
    %censor = 1;
    noFSL=0;
end
format long

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

if ~noFSL
    %Save to regs folder
    c = asterisk(x,b); %Add in asterisks and clean up data
    dlmcell([file_loc fname '.dat'],c,'delimiter','\t')
    y=0;
else
    %write the -stim_times file
    fname = [fname '_noFSL'];
    y = x(logical(x(:,3)),1)';
    %Quick fix hack for just first ten trials troubleshoot SPMG2
    %y = y(1:10);
    dlmwrite([file_loc fname '.dat'],y,'delimiter','\t','precision','%.6f')
end
return

function c = asterisk(x,b)
%adding asterisk to existing .dat files also removes any nans present

c=[];
ast = {'*', '*', '*'};
trial_index=find(b.TrialNum==1);

%Short while loop to place asterisks for regressors
while ~isempty(trial_index)    
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
