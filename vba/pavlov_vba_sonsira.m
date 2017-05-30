function [posterior,out]=pavlov_vba_sonsira(input_struct)
%NFB placebo vba modeling code

%TODO set up all the variables and defaults, is no struct is given just
%assume its the data?


%To save or not to save results
save_results=1;

%% read in data
try
    data = input_struct.data;
catch
    data = input_struct;
end

%% construct the u and y

%Infusion or no infusion
cs = (cellfun(@(x) strcmp(x,'A'), data.Infusion) | cellfun(@(x) strcmp(x,'B'), data.Infusion))'; %Inf(1) or noInf (0)

%Subj responses for will improve
will_imp=cellfun(@(x) str2double(x),cellfun(@count_subj_response, data.WillImpRespText, 'UniformOutput', false));

%Subj responses for improved
imp=cellfun(@(x) str2double(x),cellfun(@count_subj_response, data.ImprovedRespText, 'UniformOutput', false));

us = [will_imp'; imp'];

%Set up u
u = [cs 0; ...  % 1 infusion on current trial
    [us [0;0]]; ...  % 2 feedback on current trial
    0 cs; ...% 3 infusion on previous trial (leading to the current, pre-update value)
    [[0;0] us]; ... %4 feedback on previous trial
    ];

%Set up y
y = zeros(4,length(imp));
for i = 1:length(imp)
    %For will improve rating
    try
        if will_imp(i)
            y(1,i)=1;
        else
            y(2,i)=1;
        end
    catch
        %In case there is a nan just leave it as 0? Skip the nans?
    end
    
    %For improve rating
    try
        if imp(i)
            y(3,i)=1;
        else
            y(4,i)=1;
        end
    catch
        %In case there is a nan just leave it as 0? Skip the nans?
    end
end



%% define f and g functions
f_fname = @f_pavlov; % evolution function 
g_fname=@g_feedback_ratings; %observation function

%% Set up options
options.inF.noCS = 0;
options.inG.noCS = 0;

%options.isYout      = zeros(size(y)) ;
options.isYout=isnan(y) ;
options.isYout(1) = 1;
options.inF.binomial_override = 1;
options.inF.decay = decay;
options.inG.decay = decay;
options.inG.biases = biases;
options.inG.infusion_expectancy = infusion_expectancy;
options.inF.track_pe = track_pe;
options.inG.track_pe = track_pe;



dim = struct('n',2,'n_theta',1,'n_phi',2);
priors.muTheta = zeros(dim.n_theta,1);
priors.SigmaTheta = 1e1*eye(dim.n_theta);
priors.muTheta(1) = -1.3801; %% fix LR at .2


priors.muPhi = zeros(dim.n_phi,1);
priors.SigmaPhi = 1e1*eye(dim.n_phi);

priors.a_alpha = Inf;
priors.b_alpha = 0;
options.priors = priors;

%Turn graphics on or off

options.DisplayWin = graphics;
options.GnFigs = graphics;

%Set prior covariance of x_0
priors.SigmaX0 = zeros*eye(dim.n);
% priors.SigmaX0 = 1e1*eye(dim.n); %% try fitting prior expectancy

%Set prior mean for x_0
priors.muX0 = zeros(dim.n,1);


[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

%L(ct) = out.F;

if save_results
    save(sprintf('vba_data/subj_%d_vba_data.mat',tmp_id),'out','posterior')
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



