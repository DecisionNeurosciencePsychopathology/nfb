function [posterior,out]=pavlov_vba_sonrisa(vba_input,model_name)
%NFB placebo vba modeling code

%TODO set up all the variables and defaults, is no struct is given just
%assume its the data?

close all; %Get ride of all figures

%To save or not to save results
save_results=1;
save_figure=0;
%model_name = 'exp';

%To plot or not
graphics=0;

%If we want to yoke the stimuli to have the same starting parameters
yoke_stimuli = 0;
yoked_muX = 0;
one_hidden_state=0;
use_expectancy_as_prior = 1;

%% read in data
try
    data = vba_input.data;
catch
    data = vba_input;
end

%Set up simple test case for now
% load('vba/test_subj.mat')
% data = vba_test_subj;

%% construct the u and y

%Infusion or no infusion

cs = strcmp(data.Infusion,'A')' | strcmp(data.Infusion,'B')';
%cs = (cellfun(@(x) strcmp(x,'A'), data.Infusion) | cellfun(@(x) strcmp(x,'B'), data.Infusion))'; %Inf(1) or noInf (0)

%Subj responses for will improve
will_imp=cellfun(@str2double,cellfun(@count_subj_response, data.WillImpRespText, 'UniformOutput', false));

%Subj feedback
fb = strcmp(data.Feedback,'Signal')';
%fb = (cellfun(@(x) strcmp(x,'Signal'), data.Feedback))'; %Feedback = Signal or Baseline

%Subj responses for improved
imp=cellfun(@(x) str2double(x),cellfun(@count_subj_response, data.ImprovedRespText, 'UniformOutput', false));

%Concatenate responses
resp = [will_imp'; imp'];

%Set up y -- clean this up
%  y = [data.WillImpRespBin(2:end); NaN]';
%y = data.ImprovedRespBin'; % looks like you don't need to shift the y, Jon to check
y = data.WillImpRespBin';

%Create the hidden state index
cond = [strcmp(data.Infusion,'A')';
    strcmp(data.Infusion,'B')';
    strcmp(data.Infusion,'C')';
    strcmp(data.Infusion,'D')';];

%Create condition array for update indexing
condition = double(cond);
condition(2,:)=condition(2,:).*2;
condition(3,:)=condition(3,:).*3;
condition(4,:)=condition(4,:).*4;
condition=sum(condition);

%just for expected ratings for now
%  y = [y(1,:) & will_imp'==1;
%      y(2,:)& will_imp'==1;
%      y(3,:)& will_imp'==1;
%      y(4,:)& will_imp'==1;];

%Set up u - second half is shifted akin to bandit VBA!
u = [cs 0; ...  % 1 infusion on current trial
    [resp [0;0]]; ...  % 2,3 willImp resp; Imp resp
    condition 0; ... %4 Which stimulus
    fb 0; ... %5 Feedback
    0 cs; ...% 6 infusion on previous trial (leading to the current, pre-update value)
    [[0;0] resp]; ... % 7,8 willImp resp; Imp resp previous trial
    0 condition; %9 Which stim for previous trial
    0 fb; %10 feedback on previous trial
    fb(2:end) [0 0]; %11 feedback on next trial
    ];


%Delete later check if exp aligns with rating
% u(10,4) = 10;
% u(10,10) = 10;

%% define f and g functions
switch model_name
    case 'null'
        f_fname = @f_pavlov_null; % evolution function
        g_fname=@g_sigmoid_null; %observation function
        num_hidden_states=4;
        n_theta=0;
        n_phi=0;
    case 'oneLR_twoK'
        f_fname = @f_pavlov_1LR; % evolution function
        g_fname=@g_sigmoid_1B_2K; %observation function
        num_hidden_states=4;
        n_theta=2;
        n_phi=3;
    case 'oneLR_fixD_twoK'
        f_fname = @f_pavlov_1LR_fixD; % evolution function
        g_fname=@g_sigmoid_1B_2K; %observation function
        num_hidden_states=4;
        n_theta=1;
        n_phi=3;
    case 'twoLR_twoK'
        f_fname = @f_pavlov_2LR; % evolution function
        g_fname=@g_sigmoid_1B_2K; %observation function
        num_hidden_states=4;
        n_theta=3;
        n_phi=3;
    case 'twoLR_fixD'
        f_fname = @f_pavlov_2LR_fixD; % evolution function
        g_fname=@g_sigmoid_1B_2K; %observation function
        num_hidden_states=4;
        n_theta=2;
        n_phi=3;
    case 'twoLR_S_twoK'
        f_fname = @f_pavlov_2LR_S; % evolution function
        g_fname=@g_sigmoid_1B_2K; %observation function
        num_hidden_states=4;
        n_theta=4;
        n_phi=3;
    case 'twoLR_S_fixD_twoK'
        f_fname = @f_pavlov_2LR_S_fixD; % evolution function
        g_fname=@g_sigmoid_1B_2K; %observation function
        num_hidden_states=4;
        n_theta=3;
        n_phi=3;
   case 'twoLR_S_fixD_oneK'
        f_fname = @f_pavlov_2LR_S_fixD; % evolution function
        g_fname=@g_sigmoid_1B_1K; %observation function
        num_hidden_states=4;
        n_theta=3;
        n_phi=2;
    case 'oneLR_S_twoK'
        f_fname = @f_pavlov_1LR_S; % evolution function
        g_fname=@g_sigmoid_1B_2K; %observation function
        num_hidden_states=4;
        n_theta=3;
        n_phi=3;
    case 'oneLR_S_fixD_twoK'
        f_fname = @f_pavlov_1LR_S_fixD; % evolution function
        g_fname=@g_sigmoid_1B_2K; %observation function
        num_hidden_states=4;
        n_theta=2;
        n_phi=3;
    case 'twoLR_S_fixD_m_twoB_twoK'
        f_fname = @f_pavlov_2LR_S_fixD_m; % evolution function
        g_fname=@g_sigmoid_2B_2K_m; %observation function
        num_hidden_states=5;
        n_theta=4;
        n_phi=4;
    case 'momentum'
        f_fname = @f_pavlov_m; % evolution function
        g_fname=@g_sigmoid_2B_2K_m; %observation function
        num_hidden_states=5;
        n_theta=4;
        n_phi=4;
    case 'momentum_oneK'
        f_fname = @f_pavlov_m; % evolution function
        g_fname=@g_sigmoid_2B_1K_m; %observation function
        num_hidden_states=5;
        n_theta=4;
        n_phi=3;
    otherwise
        error('Improper model selected')
end
%% Set up options
options.sources(1).type=1;
options.sources(1).out  = 1;

%options.isYout      = zeros(size(y)) ;
options.isYout=isnan(y) ;
%willImpMissedTrials= [isnan(u(2,:)) | isnan(u(6,:))];
%options.isYout(:,willImpMissedTrials')=1;
options.isYout(:,1) = 1;
%options.inF.binomial_override = 1;

% skip first trial
options.skipf = zeros(1,length(y));
options.skipf(1) = 1;

%Turn off measurement noise
options.binomial = 1 ;

% dim = struct('n',4,'n_theta',1,'n_phi',1); % no decay version
if one_hidden_state
    num_hidden_states=1;
    model_name='twoLR_S_fixD_m';
end


%We want to track PE
num_hidden_states = num_hidden_states + 1;


dim = struct('n',num_hidden_states,'n_theta',n_theta,'n_phi',n_phi);
priors.muTheta = zeros(dim.n_theta,1);
priors.SigmaTheta = 1e1*eye(dim.n_theta);
%priors.muTheta(1) = -1.3801; %% fix LR at .2


priors.muPhi = zeros(dim.n_phi,1);
priors.SigmaPhi = 1e1*eye(dim.n_phi);

priors.a_alpha = Inf;
priors.b_alpha = 0;

%Turn graphics on or off

options.DisplayWin = graphics;
%options.GnFigs = graphics;

%Set prior covariance of x_0
%priors.SigmaX0 = zeros*eye(dim.n);
priors.muX0 = zeros(dim.n,1);

if yoke_stimuli
    model_name = 'yoked';
    priors.SigmaX0 = 1e1*eye(dim.n);  %% try fitting prior expectancy
    
    %A and B
    priors.SigmaX0(2,1) = priors.SigmaX0(1,1);
    priors.SigmaX0(1,2) = priors.SigmaX0(1,1);
    
    %C and D
    priors.SigmaX0(3,4) = priors.SigmaX0(1,1);
    priors.SigmaX0(4,3) = priors.SigmaX0(1,1);
elseif yoked_muX
    %     model_name = 'yoked_muX';
    priors.muX0 = [0 0 0 0 0]';
    priors.SigmaX0 = zeros(dim.n);
else
    priors.SigmaX0 = 1e1*eye(dim.n); %% try fitting prior expectancy
end


if use_expectancy_as_prior
    try
        %Which protocol are we in - todo: make path generic
        if strcmp(vba_input.protocol,'SON1')
            expectancies = readtable('/Users/martapecina/Box Sync/PITT/RESEARCH/fMRI_shared/NFB/NFB_response/SON1&2_behav_results/Expectancies.csv');
        else
            expectancies = readtable('/Users/martapecina/Box Sync/PITT/RESEARCH/fMRI_shared/NFB/NFB_response/SON1&2_behav_results/Expectancies.csv');
        end
        
        %Normalize the exp. ratings if the DNE simply make them 0
        expectancies.expectations_q1 = expectancies.expectations_q1./3;
        expectancies.expectations_q2 = expectancies.expectations_q2./3;
        expectancies.expectations_q1(isnan(expectancies.expectations_q1)) = 0;
        expectancies.expectations_q2(isnan(expectancies.expectations_q2)) = 0;
        
        %Pull the subj id to match the csv
        rating_id = [vba_input.protocol, '_', vba_input.subj_name(end-2:end)];
        
        %Create logical arrays to index exp ratings csv
        subj_idx = strcmp(expectancies.subject_id,rating_id);
        if iscell(expectancies.admin)
            admin_idx = ~cellfun(@isempty,(strfind(expectancies.admin,vba_input.admin(end))));
        else
            admin_idx = expectancies.admin==str2num(vba_input.admin(end));
        end
        
        %Compile data
        a_b_x0 = expectancies.expectations_q1(subj_idx & admin_idx);
        c_d_x0 = expectancies.expectations_q2(subj_idx & admin_idx);
        
        if isempty(a_b_x0) || isempty(c_d_x0)
            error('Ratings returned empty check!')
        end
            
        %Apply to priors (A B C D PE)
        additional_hidden_states = num_hidden_states - 4; 
        priors_to_add = zeros(1,additional_hidden_states);
        priors.muX0 = [a_b_x0 a_b_x0 c_d_x0 c_d_x0 priors_to_add]';
        
        %Inf. percision priors for sigma
        priors.SigmaX0 = zeros(dim.n);
    catch
        error('EXP matching to priors failed')
    end
    
end

options.priors = priors;

%Set prior mean for x_0
%priors.muX0 = zeros(dim.n,1);

%Set up options for evolution and observation functions
options.inF.num_hidden_states = num_hidden_states;
options.inG.num_hidden_states = num_hidden_states;

%% skip first trial
options.skipf = zeros(1,length(y));
options.skipf(1) = 1;

%% Run main f(x)
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

%L(ct) = out.F;

if save_results
    %TODO - make this save path generic
    save_path = ['/Users/martapecina/GitHub/Nfb_sonrisa/vba_data/', model_name];
    if ~exist(save_path)
        mkdir(save_path)
    end
    
    save([save_path sprintf('/subj_%s_%s_vba_data_%s.mat',vba_input.subj_name, vba_input.admin, model_name)],'out','posterior')
end

if save_figure
    h1 = figure(1);
    savefig(h1,sprintf('vba_data/subj_%s_%s_vba_data_%s.fig',vba_input.subj_name, vba_input.admin, model_name))
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



