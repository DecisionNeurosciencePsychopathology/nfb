
function  [V_tplus1] = f_pavlov(V_t,theta,u_t,inF )
% function  [ fx,dfdx,dfdP ] = f_pavlov( x_t,P,u_t,in )
% IN:
% - x_t : Q-values (2*1)
% - P : learning rate (1*1)
% - u_t : previous action and feedback
% - in : []

alpha_pos = 1./(1+exp(-theta(1)));          % learning rate
alpha_base = 1./(1+exp(-theta(2)));          % learning rate
% gamma = 1./(1+exp(-theta(2)));          % decay
gamma = 0.99;                                 % decay

V_tplus1= zeros(1,inF.num_hidden_states);


% expectation = u_t(6); %will improve response

%If we are running the one hidden state model
if inF.num_hidden_states==1
    vt_to_update=1; %Hacky...
else
    vt_to_update = u_t(9); %Index of which condition to update;
end


%Grab the feedback
fb = u_t(10); 

%Update the proper stimuli presented
update_vector = zeros(inF.num_hidden_states,1);
update_vector(vt_to_update) = 1;
update_vector = logical(update_vector); %Needs to be a logical

%Update the current condition present to the subject
if fb>0
    V_tplus1(update_vector) = V_t(update_vector) + alpha_pos*(fb-V_t(update_vector));
elseif fb==0
    V_tplus1(update_vector) = V_t(update_vector) + alpha_base*(fb-V_t(update_vector));
end


%Decay the other conditions, if there are other conditions
if inF.num_hidden_states > 1
    V_tplus1(~update_vector) = gamma.*V_t(~update_vector); 
end

%Account for nans and infs
if any(isnan(V_tplus1)) || any(isinf(V_tplus1))
    V_tplus1=V_t;
end


