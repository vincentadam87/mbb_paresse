function [ fx ] = f_paresse_MATHYS(x,P,u,in)
% INPUT
% - x : Qvalues
%     1:5 : pRleft
%     6:10 : pRright
%     11:15 : pEleft
%     16:20 : pEright
% - P : parameters
%     1/2/3/4 : kappa/omegaR/omegaE/theta
% - u : input
%     1 : choice (left = 0, right = 1)
%     2 : reward
%     3 : effort

fx = x; % copy all

a = u(1); % action that was selected on the previous trial
R = u(2); % reward obtained ( = reward category )
E = u(3); % effort done ( = effort category )
    
fx = x;


% REWARD
ur = [5*a+1:5*(a+1)];
fx(ur) = f_VBvolatile_1p(x(ur),P([1,2,4]),R,in);

% EFFORT
ue = 10+[5*a+1:5*(a+1)];
fx(ue) = f_VBvolatile_1p(x(ue),P([1,3,4]),E,in);


end

