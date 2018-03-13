
function  [V_t] = f_pavlov(V_tminus1,theta,u_t,inF )
% function  [ fx,dfdx,dfdP ] = f_pavlov( x_t,P,u_t,in )
% IN:
% - x_t : Q-values (2*1)
% - P : learning rate (1*1)
% - u_t : previous action and feedback
% - in : []

%STIMULUS TO UPDATE
vt_to_update = u_t(4); %Index of which condition to update;

%Update the proper stimuli presented
update_vector = zeros(inF.num_hidden_states,1);
update_vector(vt_to_update) = 1;
update_vector = logical(update_vector); %Needs to be a logical

%PARAMETERS
alpha_pos = 1./(1+exp(-theta(1)));           % learning rate
alpha_base = 1./(1+exp(-theta(2)));          % learning rate
% s_pos = 1./(1+exp(-theta(3)));               % sensitivity to feedback
gamma = 1./(1+exp(-theta(4)));               % decay

s_pos = 1;

%Decay the other conditions, if there are other conditions
if inF.num_hidden_states > 1
   V_t(~update_vector) = gamma.*V_tminus1(~update_vector); 
end

%REWARD
fb = u_t(5); 

%Q-LEARNING FUNCTION: Vt = Vt-1 + ? (s*r ? Vt-1)
if fb>0
    V_t(update_vector) = V_tminus1(update_vector) + alpha_pos*(s_pos*fb-V_tminus1(update_vector));
elseif fb==0
    V_t(update_vector) = V_tminus1(update_vector) + alpha_base*(fb-V_tminus1(update_vector));
end

%Account for nans and infs
if any(isnan(V_t)) || any(isinf(V_t))
    V_t=V_t;
end


