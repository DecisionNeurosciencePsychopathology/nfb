function [gx] = g_sigmoid_inf_no_inf(V_t,Phi,u_t,inG)
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


if u_t(1) == 1
    vt_to_update = 1; %Index of which condition to update;
else
    vt_to_update = 2;
end


x = V_t(vt_to_update);

% gx = sigm(x,[],phi); 
gx = sig(x*beta + x*kappa);
%p_choice = (exp((x-max(x))/beta)) / (sum(exp((x-max(x))/beta))); %Divide by temperature
%gx = p_choice';

%Account for nans and infs
if any(isnan(gx))>0 || any(isinf(gx))
    stop=1;
end