function write_nfb_regressors(data,id)

reg_file_dest = '/Users/martapecina/GitHub/nfb_analysis/regs/';

%% Trial
write3Ddeconv_startTimes(reg_file_dest,data.InfOnset,data.J2Onset,sprintf('%s_trial',id),true(length(data.TrialNum),1),0,data);

%% Infusion
write3Ddeconv_startTimes(reg_file_dest,data.InfOnset,data.WillImpOnset,sprintf('%s_infusion_onset',id),true(length(data.TrialNum),1),0,data);

%% Will Improve Responses
write3Ddeconv_startTimes(reg_file_dest,data.WillImpOnset,data.WillImpOnset+data.WillImpRt,sprintf('%s_willImpRT',id),true(length(data.TrialNum),1),0,data);

%% Feedback
write3Ddeconv_startTimes(reg_file_dest,data.Feed1Onset,data.Feed3Onset,sprintf('%s_feedback_onset',id),true(length(data.TrialNum),1),0,data);

%% Improve Responses
write3Ddeconv_startTimes(reg_file_dest,data.ImprovedOnset,data.ImprovedOnset+data.ImprovedRt,sprintf('%s_impRT',id),true(length(data.TrialNum),1),0,data);


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
