function write_nfb_regressors_2(data,id)

%% Datalocation -- set as needed
if ispc
    reg_file_dest = 'C:/kod/nfb/analysis/regs/'; %Change as needed
else
    reg_file_dest = '/Users/martapecina/GitHub/Nfb_sonrisa/regs/';
end


%% Pull variables from dataset (maybe do this in parent script)
%Infusion (1) no infusion (-1)
inf_no_inf = double(cellfun(@(x) strcmp(x,'A'), data.Infusion) | cellfun(@(x) strcmp(x,'B'), data.Infusion));
inf_no_inf(inf_no_inf==0) = -1;

%Feedback (1) no feedback (-1) is signal vs baseline
feedback = double(cellfun(@(x) strcmpi(x,'Signal'), data.Feedback));
feedback(feedback==0) = -1;

%Will improve resposnes yes (1) no (-1)
willImpResp = double(cellfun(@(x) strcmpi(x,'yes'), data.WillImpRespText));
willImpResp(willImpResp==0) = -1;

%Will improve resposnes yes (1) no (-1)
impResp = double(cellfun(@(x) strcmpi(x,'yes'), data.ImprovedRespText));
impResp(impResp==0) = -1;

%% Trial
write3Ddeconv_startTimes(reg_file_dest,data.InfOnset,data.J2Onset,sprintf('%s_trial',id),true(length(data.TrialNum),1),0,data);

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
% 
motor_offset(1:2:end)=data.WillImpOnset+data.WillImpRt; %This is for contrasts
motor_offset(2:2:end)=data.ImprovedOnset+data.ImprovedRt; %This is for contrasts

 %[motor_left_onset,motor_left_offset]=create_motor_reg(data,motor_left);
 %[motor_right_onset,motor_right_onset]=create_motor_reg(data,motor_right);

 %Create new binary regressors
left_only = motor_left;
left_only(left_only==0)=nan;

right_only = motor_right;
right_only(right_only==0)=nan;

%Override trial index
motor_trial_index = [1 65 129 193];

write3Ddeconv_startTimes(reg_file_dest,motor_onset,motor_offset,sprintf('%s_left_only_all',id),left_only,0,data,motor_trial_index);
write3Ddeconv_startTimes(reg_file_dest,motor_onset,motor_offset,sprintf('%s_right_only_all',id),right_only,0,data,motor_trial_index);
 
%write3Ddeconv_startTimes(reg_file_dest,motor_left_onset,motor_left_offset,sprintf('%s_left_all',id),motor_left,0,data);
%write3Ddeconv_startTimes(reg_file_dest,motor_left_onset,motor_left_offset,sprintf('%s_right_all',id),motor_right,0,data);

%Strip it down use only one motor event 
%write3Ddeconv_startTimes(reg_file_dest,data.ImprovedOnset,data.ImprovedOnset+data.ImprovedRt,sprintf('%s_left',id),impLeft,0,data);
%write3Ddeconv_startTimes(reg_file_dest,data.ImprovedOnset,data.ImprovedOnset+data.ImprovedRt,sprintf('%s_right',id),~impLeft,0,data);


%% Infusion 1 sec stick
write3Ddeconv_startTimes(reg_file_dest,data.InfOnset,data.InfOnset+1,sprintf('%s_infusion_onset',id),true(length(data.TrialNum),1),0,data);

% Infusion aligned from onset to end of willImp RT
write3Ddeconv_startTimes(reg_file_dest,data.InfOnset,data.WillImpOnset+data.WillImpRt,sprintf('%s_infusion_onsetToRT',id),true(length(data.TrialNum),1),0,data);

%% Infusion / no infusion 
write3Ddeconv_startTimes(reg_file_dest,data.InfOnset,data.WillImpOnset,sprintf('%s_infusion_no_infusion',id),inf_no_inf,0,data);

%Align inf no inf with RT
write3Ddeconv_startTimes(reg_file_dest,data.WillImpOnset,data.WillImpOnset+data.WillImpRt,sprintf('%s_infusion_no_infusion_RTconvolv',id),inf_no_inf,0,data);

%Align inf no inf from onset to end of RT
write3Ddeconv_startTimes(reg_file_dest,data.InfOnset,data.WillImpOnset+data.WillImpRt,sprintf('%s_infusion_no_infusion_onsetToRT',id),inf_no_inf,0,data);

%% Feedback / no Feedback -- determine which feedback gives the best signal -- it is fb2
write3Ddeconv_startTimes(reg_file_dest,data.Feed2Onset,data.Feed2Onset+1,sprintf('%s_feedback_no_feedback',id),feedback,0,data);

%Align feedback no feedback with RT
write3Ddeconv_startTimes(reg_file_dest,data.ImprovedOnset,data.ImprovedOnset+data.ImprovedRt,sprintf('%s_feedback_no_feedback_RTconvolv',id),feedback,0,data);

%Align feedback no feedback from feedback onset to the end of RT
write3Ddeconv_startTimes(reg_file_dest,data.Feed2Onset,data.ImprovedOnset+data.ImprovedRt,sprintf('%s_feedback_no_feedback_onsetToRT',id),feedback,0,data);

%% Will Improve Responses EVENT TIMES ONLY - RT convolv
write3Ddeconv_startTimes(reg_file_dest,data.WillImpOnset,data.WillImpOnset+data.WillImpRt,sprintf('%s_willImpRespEvent',id),true(length(data.TrialNum),1),0,data);
write3Ddeconv_startTimes(reg_file_dest,data.WillImpOnset,data.WillImpOnset+1,sprintf('%s_willImpRespEventStick',id),true(length(data.TrialNum),1),0,data);

%% Will Improve Responses ACTUAL RESPONSES AS PARAMETERIC MOD
write3Ddeconv_startTimes(reg_file_dest,data.WillImpOnset,data.WillImpOnset+data.WillImpRt,sprintf('%s_willImpRespMod',id),willImpResp,0,data);

%% Feedback 1 sec stick - start at inflection of feedback signal
%Feed 1 onset is the start of the feedback screen
%Feed 2 onset is the start of inflection (1.5-3.3 secs)
%Feed 3 osnet is the point when the inflection reaches its peak
write3Ddeconv_startTimes(reg_file_dest,data.Feed1Onset,data.Feed1Onset+1,sprintf('%s_feedback1_onset',id),true(length(data.TrialNum),1),0,data);
write3Ddeconv_startTimes(reg_file_dest,data.Feed2Onset,data.Feed2Onset+1,sprintf('%s_feedback2_onset',id),true(length(data.TrialNum),1),0,data);
write3Ddeconv_startTimes(reg_file_dest,data.Feed3Onset,data.Feed3Onset+1,sprintf('%s_feedback3_onset',id),true(length(data.TrialNum),1),0,data);
write3Ddeconv_startTimes(reg_file_dest,data.Feed1Onset,data.ImprovedOnset,sprintf('%s_feedback_all',id),true(length(data.TrialNum),1),0,data);

%Align feedback from feedback onset to end of RT
write3Ddeconv_startTimes(reg_file_dest,data.Feed2Onset,data.ImprovedOnset+data.ImprovedRt,sprintf('%s_feedback2_onsetToRT',id),true(length(data.TrialNum),1),0,data);

%% Improve Responses EVENT TIMES ONLY - RT convolv
write3Ddeconv_startTimes(reg_file_dest,data.ImprovedOnset,data.ImprovedOnset+data.ImprovedRt,sprintf('%s_impRespEvent',id),true(length(data.TrialNum),1),0,data);
write3Ddeconv_startTimes(reg_file_dest,data.ImprovedOnset,data.ImprovedOnset+1,sprintf('%s_impRespEventStick',id),true(length(data.TrialNum),1),0,data);

%% Improve Responses ACTUAL RESPONSES AS PARAMETERIC MOD
write3Ddeconv_startTimes(reg_file_dest,data.ImprovedOnset,data.ImprovedOnset+data.ImprovedRt,sprintf('%s_impRespMod',id),impResp,0,data);

%% Censor
%Only sonrisa 1 right nbow
%if strfind(id,'SON1')
try
    create_censor_file_2_1(data,id)
catch
    fprintf('\n%s could not construct Censor file! Check proc dir...\n')
end
%end


function [x,y]=write3Ddeconv_startTimes(file_loc,event_beg,event_end,fname,modulator,noFSL,b,override)
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

if ~noFSL
    %Save to regs folder
    c = asterisk(x,b,override); %Add in asterisks and clean up data
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

