function create_censor_file_2_1(data,id)
%This is a file of 1s and 0s, indicating which
%time points are to be included (1) and which are
%to be excluded (0).

%Load in volume file -- have to fix this for sonrisa 2 later !!!!
if strfind(id,'SON1')>0
    vol_T = readtable('vol_info/fd_vol_output_sonrisa1.csv');
else
    vol_T = readtable('vol_info/fd_vol_output_sonrisa2.csv');
end

%Set some constants grab the block lenths
frequency_scale_hz = 10;
scan_tr = 1; %According to 3dinfo -verb!

%Try to keep vol lookup table updated
if any(strcmp(vol_T.Subjects,id))
    block_lengths = vol_T(strcmp(vol_T.Subjects,id),:);
    block_lengths = table2array(block_lengths(1,2:end));
else
    %However if subjects are not in vol_info we assume all runs are
    %completed with all volumes intact
    block_lengths = [672 672 672 672];
end


%Set the onset time
stim_OnsetTime = data.InfOnset*1000; %Convert to ms
stim_NextOnsetTime = []; %Hacky work around to get matrices to play nice

%Set logical trial censor
%data.trials_to_censor = ~data.subjects.(data.id).run_censor_full; %DID I HAVE THIS BACKWARDS!? there was no '~' before 11/4

% % %I was actually censoring the 0 trials :( let's try what we originally
% % %thought was going on
% % data.trials_to_censor = zeros(length(data.subjects.(data.id).run_censor_full),1);


%Censor only the NANs
missed_trials = zeros(length(data.WillImpRespNum)+length(data.ImprovedRespNum),1);
missed_trials(1:2:end)=data.WillImpRespNum;
missed_trials(2:2:end)=data.ImprovedRespNum;
missed_trials = isnan(missed_trials);

%Need to loop over block lengths as they are all different for each run.
num_blocks = unique(data.Run);

for block_n = 1:length(num_blocks)
    %Set up trial ranges
    trial_idx = data.Run==block_n;
    %     trial_index_1 = data.trial_index(block_n);
    %     trial_index_2 = trial_index_1 + data.trials_per_block-1;
    
    stim_OnsetTime_tmp = stim_OnsetTime(trial_idx);
    stim_NextOnsetTime = [0; stim_OnsetTime_tmp(2:end)];
    %stim_NextOnsetTime = [stim_NextOnsetTime; stim_OnsetTime(logical([0; trial_idx(2:end)]))];
    
    bin_size = 1/frequency_scale_hz*1000; % convert Hz to msec
    %epoch_window = data.stim_OnsetTime(trial_index_1:trial_index_2):bin_size:data.stim_OnsetTime(trial_index_1:trial_index_2)+scan_tr*block_lengths*1000;
    epoch_window = stim_OnsetTime_tmp:bin_size:stim_OnsetTime_tmp+scan_tr*block_lengths(block_n)*1000;
    event_beg = stim_OnsetTime_tmp; event_end = stim_NextOnsetTime;
    
    tmp_reg.(['regressors' num2str(block_n)]).to_censor = ...
        createSimpleRegressor(event_beg, event_end, epoch_window, missed_trials(trial_idx));
    tmp_reg.(['regressors' num2str(block_n)]).to_censor = ones(size(tmp_reg.(['regressors' num2str(block_n)]).to_censor)) - tmp_reg.(['regressors' num2str(block_n)]).to_censor;
    
    
    % NB: the first 5s are censored because they capture HRF to events
    % preceding the first trial
    tmp_reg.(['hrfreg' num2str(block_n)]).to_censor = ...
        gsresample( ...
        [zeros(50,1)' tmp_reg.(['regressors' num2str(block_n)]).to_censor(1:end-51)], ...
        10,1./scan_tr);
    
end


fnm = fieldnames(tmp_reg.regressors1)';
ct=1:length(fnm);
hrf_regs.to_censor = [];
for i = 1:length(num_blocks)
    try
        hrf_regs.to_censor = [hrf_regs.to_censor tmp_reg.(['hrfreg' num2str(i)]).to_censor];
    catch
        error('Censor is bad new bears, check the number of runs')
    end
end
hrf_regs.to_censor = 1-(ceil(hrf_regs.to_censor));
hrf_regs.to_censor = ~hrf_regs.to_censor;

%Write it to file
gdlmwrite(['regs/' sprintf('%s_CensorOnly.regs',id)],hrf_regs.to_censor');


function foo = createSimpleRegressor(event_begin,event_end,epoch_window,conditional_trials)
% this was not a problem earlier, but for some reason it is now: find indices that would
% result in a negative value and set them both to 0
qbz = ( event_begin == 0 ); qez = ( event_end == 0 );
event_begin( qbz | qez ) = 0; event_end( qbz | qez ) = 0;

% check if optional censoring variable was used
if(~exist('conditional_trials','var') || isempty(conditional_trials))
    conditional_trials = true(length(event_begin),1);
elseif(~islogical(conditional_trials))
    % needs to be logical format to index cells
    conditional_trials = logical(conditional_trials);
end

% this only happened recently, but it's weird
if(any((event_end(conditional_trials)-event_begin(conditional_trials)) < 0))
    error('MATLAB:bandit_fmri:time_travel','feedback is apparently received before RT');
end

% create epoch windows for each trial
epoch = arrayfun(@(a,b) a:b,event_begin,event_end,'UniformOutput',false);
foo = zeros(size(epoch_window));

for n = 1:numel(epoch)
    if(conditional_trials(n))
        foo = logical(foo + histc(epoch{n},epoch_window));
    end
end

return