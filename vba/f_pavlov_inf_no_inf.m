
function  [V_tplus1] = f_pavlov_inf_no_inf(V_t,theta,u_t,inF )
% function  [ fx,dfdx,dfdP ] = f_pavlov( x_t,P,u_t,in )
% IN:
% - x_t : Q-values (2*1)
% - P : learning rate (1*1)
% - u_t : previous action and feedback
% - in : []

alpha_pos = 1./(1+exp(-theta(1)));          % learning rate
alpha_base = 1./(1+exp(-theta(2)));          % learning rate
s_pos = 1./(1+exp(-theta(3)));                   % sensitivity to feedback
s_base = 1./(1+exp(-theta(4)));             % sensitivity to feedback
gamma = 1./(1+exp(-theta(5)));               % decay

V_tplus1= zeros(1,inF.num_hidden_states);

% expectation = u_t(6); %will improve response

%If we are running the one hidden state model
if u_t(1) == 1
    vt_to_update = 1; %Index of which condition to update;
else
    vt_to_update = 2;
end

%Grab the feedback
fb = u_t(11); 

% %Update the proper stimuli presented
% update_vector = zeros(inF.num_hidden_states,1);
% update_vector(vt_to_update) = 1;
% update_vector = logical(update_vector); %Needs to be a logical

%Update the current condition present to the subject
if fb>0
    V_tplus1(vt_to_update) = V_t(vt_to_update) + alpha_pos*(s_pos*fb-V_t(vt_to_update));
elseif fb==0
    V_tplus1(vt_to_update) = V_t(vt_to_update) + alpha_base*(s_base*fb-V_t(vt_to_update));
end

% fx(1) = x(1) + alpha*pe(1)*u(1);
% fx(2) = x(2) + alpha*pe(2)*(1-u(1));

%Decay the other conditions, if there are other conditions
if inF.num_hidden_states > 1
    V_tplus1(vt_to_update) = gamma.*V_t(vt_to_update); 
end

%Account for nans and infs
if any(isnan(V_tplus1)) || any(isinf(V_tplus1))
    V_tplus1=V_t;
end


