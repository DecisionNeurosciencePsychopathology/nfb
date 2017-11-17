function  [Vst] = f_pavlov(Vst_1,theta,u_t,in )
% function  [ fx,dfdx,dfdP ] = f_pavlov( x_t,P,u_t,in )
% IN:
% - x_t : Q-values (2*1)
% - P : learning rate (1*1)
% - u_t : previous action and feedback
% - in : []

alpha = 1./(1+exp(-theta(1)));          % learning rate
gamma = 1./(1+exp(-theta(2)));          % decay
Vst = zeros(1,4);


% expectation = u_t(6); %will improve response
vt_to_update = u_t(9); %Index of which condition to update;


%Grab the feedback
fb = u_t(10); 

%Update the proper stimuli presented
update_vector = zeros(4,1);
update_vector(vt_to_update)=1;
update_vector = logical(update_vector); %Needs to be a logical

%Update the current condition present to the subject
Vst(update_vector) = Vst_1(update_vector) + alpha*(fb-Vst_1(update_vector));

%Decay the other conditions
Vst(~update_vector) = gamma.*Vst_1(~update_vector); 

if any(isnan(Vst)) || any(isinf(Vst))
    Vst=Vst_1;
end


