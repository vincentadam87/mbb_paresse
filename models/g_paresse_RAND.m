function [ gx ] = g_paresse_RAND(x,P,u,in)
% INPUT
% - x : Qvalues
% - P : parameters
%     1 : bias
% - u : input
%     1 : choice (left = 0, right = 1)
%     2 : reward
%     3 : effort

side_bias = P(1);
cond_bias = P(2);                      
sess_bias = P(3);                      

gx = 1/(1+exp(side_bias+cond_bias+sess_bias)); % bias through sigmoid mapping

end