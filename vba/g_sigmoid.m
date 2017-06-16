function [gx] = g_sigmoid(Vst,Phi,u_t,inG)
% Identity observation mapping (partially observable)

omega = Phi(1); %Free paramter no contraints?

gx = sigmoid((Vst+omega))';