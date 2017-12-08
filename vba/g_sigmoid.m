function [gx] = g_sigmoid(V_t,Phi,u_t,inG)
% Identity observation mapping (partially observable)

% omega = Phi(1); %Free paramter no contraints?
% 
% gx = sigmoid((Vst+omega))';


% INPUT
% - x : Q-values (2x1)
% - P : inverse temperature (1x1)
% - u : [useless]
% - in : [useless]
% OUTPUT
% - gx : P(a=1|x)

beta = exp(Phi(1)); % inverse temperature (stochasticity of ratings)
kappa = Phi(2)./10; % subject's choice bias

if inG.num_hidden_states==1
    current_stim=1; %Hacky...
else
    current_stim = u_t(4);
end

x = V_t(current_stim);

% gx = sigm(x,[],phi); 
gx = sig(x*beta + kappa);
%p_choice = (exp((x-max(x))/beta)) / (sum(exp((x-max(x))/beta))); %Divide by temperature
%gx = p_choice';

%Account for nans and infs
if any(isnan(gx))>0 || any(isinf(gx))
    stop=1;
end