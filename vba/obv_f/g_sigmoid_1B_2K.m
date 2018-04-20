function [gx] = g_sigmoid_1B_2K(V_t,Phi,u_t,inG)
% Identity observation mapping (partially observable)
% omega = Phi(1); %Free paramter no contraints?
% gx = sigmoid((Vst+omega))';

% INPUT
% - x : Q-values (2x1)
% - P : inverse temperature (1x1)
% - u : [useless]
% - in : [useless]
% OUTPUT
% - gx : P(a=1|x)

%PARAMETERS
beta = exp(Phi(1)); % inverse temperature (stochasticity/decisiveness of ratings)
%kappa = Phi(2)./10; % subject's choice bias
kappa1 = Phi(2)./10;  %placebo bias parameter
kappa2 = Phi(3)./10; %calibration bias parameter

choice_vector = u_t(4);

x = V_t(choice_vector);

%PROBABILITY FUNCTION: gx = sigm(x,[],phi); Probability of chosing one
%response after updating the value of the stimulus.
placebo = u_t(1);
if placebo==1
gx = sig(x*beta + kappa1);
else
gx = sig(x*beta + kappa2);
end
%Account for nans and infs
if any(isnan(gx))>0 || any(isinf(gx))
    stop=1;
end