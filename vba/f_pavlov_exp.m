
function  [V_t] = f_pavlov_exp(V_tminus1,theta,u_t,inF )
% function  [ fx,dfdx,dfdP ] = f_pavlov( x_t,P,u_t,in )
% IN:
% - x_t : Q-values (2*1)
% - P : learning rate (1*1)
% - u_t : previous action and feedback
% - in : []

%STIMULUS TO UPDATE
vt_to_update = u_t(9); %Index of which condition to update;

%Update the proper stimuli presented
update_vector = zeros(inF.num_hidden_states,1);
update_vector(vt_to_update) = 1;
update_vector = logical(update_vector); %Needs to be a logical

%PARAMETERS
alpha = 1./(1+exp(-theta(1)));           % learning rate
% alpha_base = 1./(1+exp(-theta(2)));          % learning rate
% s_pos = 1./(1+exp(-theta(3)));               % sensitivity to feedback
gamma = 1./(1+exp(-theta(2)));               % decay



%Decay the other conditions, if there are other conditions
if inF.num_hidden_states > 1
   V_t(~update_vector) = gamma.*V_tminus1(~update_vector); 
end

%REWARD on previous trial
fb = u_t(10); 

%Compute prediction error
delta = fb-V_tminus1(update_vector);

%Q-LEARNING FUNCTION: Vt = Vt-1 + ? (s*r ? Vt-1)
    V_t(update_vector) = V_tminus1(update_vector) + alpha*(delta);

%Account for nans and infs
if any(isnan(V_t)) || any(isinf(V_t))
    V_t=V_t;
end

%Assume the last hidden state is the prediction error value
V_t(end) = delta;


