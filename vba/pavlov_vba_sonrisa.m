function [posterior,out]=pavlov_vba_sonrisa(vba_input)
%NFB placebo vba modeling code

%TODO set up all the variables and defaults, is no struct is given just
%assume its the data?


%To save or not to save results
save_results=1;
model_name = 'vanilla';


%To plot or not
graphics=1;

%If we want to yoke the stimuli to have the same starting parameters
yoke_stimuli = 0;
yoke_muX = 0;
one_hidden_state=1;

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
cs = (cellfun(@(x) strcmp(x,'A'), data.Infusion) | cellfun(@(x) strcmp(x,'B'), data.Infusion))'; %Inf(1) or noInf (0)

%Subj responses for will improve
will_imp=cellfun(@(x) str2double(x),cellfun(@count_subj_response, data.WillImpRespText, 'UniformOutput', false));

%Subj responses for improved
imp=cellfun(@(x) str2double(x),cellfun(@count_subj_response, data.ImprovedRespText, 'UniformOutput', false));

%Subj feedback 
fb = (cellfun(@(x) strcmp(x,'Signal'), data.Feedback))'; %Feedback = Signal or Baseline

resp = [will_imp'; imp'];
%Set up y -- clean this up 
%  y = [data.WillImpRespBin(2:end); NaN]';
 y = data.WillImpRespBin'; % looks like you don't need to shift the y, Jon to check
 
 
 %Create the hidden state index
 cond = [(cellfun(@(x) strcmp(x,'A'), data.Infusion))';
     (cellfun(@(x) strcmp(x,'B'), data.Infusion))';
     (cellfun(@(x) strcmp(x,'C'), data.Infusion))';
     (cellfun(@(x) strcmp(x,'D'), data.Infusion))';];
 
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
    condition 0; ...
    fb 0; ...
    0 cs; ...% 4 infusion on previous trial (leading to the current, pre-update value)
    [[0;0] resp]; ... % 5,6 willImp resp; Imp resp previous trial
    0 condition;
    0 fb; %feedback on previous trial
    fb(2:end) [0 0]; % feedback on next trial
    ]; 
 

%Delete later check if exp aligns with rating
% u(10,4) = 10;
% u(10,10) = 10;

%% define f and g functions
f_fname = @f_pavlov; % evolution function 
g_fname=@g_sigmoid; %observation function

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
    model_name='one_hidden_state';
else
    num_hidden_states=4;
end

dim = struct('n',num_hidden_states,'n_theta',2,'n_phi',2);
priors.muTheta = zeros(dim.n_theta,1);
priors.SigmaTheta = 1e1*eye(dim.n_theta);
%priors.muTheta(1) = -1.3801; %% fix LR at .2


priors.muPhi = zeros(dim.n_phi,1);
priors.SigmaPhi = 1e1*eye(dim.n_phi);

priors.a_alpha = Inf;
priors.b_alpha = 0;
options.priors = priors;

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
elseif yoke_muX
    model_name = 'yoked_muX';
    priors.muX0 = [0.5 0.5 0 0]';
    priors.SigmaX0 = zeros(dim.n);
else
    priors.SigmaX0 = 1e1*eye(dim.n); %% try fitting prior expectancy
end
%Set prior mean for x_0
%priors.muX0 = zeros(dim.n,1);

%Set up options for evolution and observation functions
options.inF.num_hidden_states = num_hidden_states;
options.inG.num_hidden_states = num_hidden_states;



[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

%L(ct) = out.F;
 
% if save_results
%    save(sprintf('vba_data/subj_%s_%s_vba_data.mat',vba_input.subj_name, vba_input.admin),'out','posterior')
%    
% end

save(sprintf('vba_data/subj_%s_%s_vba_data_%s.mat',vba_input.subj_name, vba_input.admin, model_name),'out','posterior')

%--------------------------------------------------------------------------
function y = count_subj_response(str)
if strcmpi('yes',str)
    y = '1';
elseif strcmpi('no',str)
    y = '0';
else
    y = str;
end



