function [gx] = g_sigmoid(Vst,Phi,u_t,inG)
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

beta = exp(Phi(1));

x = Vst(u_t(11)); %Hardcoded = bad! Feedback on previous trial

gx(u_t(11)) = sigm(x,[],beta); %Update only the previous stimuli?

%p_choice = (exp((x-max(x))/beta)) / (sum(exp((x-max(x))/beta))); %Divide by temperature
%gx = p_choice';

if any(isnan(gx))>0 || any(isinf(gx))
    stop=1;
end