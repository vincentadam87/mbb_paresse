function [ gx ] = g_paresse_RL(x,P,u,in)
% INPUT
% - x : Qvalues
%     1 : QRleft
%     2 : QRright
%     3 : QEleft
%     4 : QEright
% - P : parameters
%     1 : inverse temperature in softmax
%     2 : effort discount factor left
%     3 : bias
% - u : input
%     1 : choice (left = 0, right = 1)
%     2 : reward
%     3 : effort

beta = exp(P(1)); % inverse temperature in softmax positive
bias = P(3);                     
    
    if isequal(in.DECIDING,'LINEAR')
        d = (P(2)); % any value
        Vleft = x(1)-d*x(3);
        Vright = x(2)-d*x(4);
    elseif isequal(in.DECIDING,'HYPERBOLIC')
        k = (P(2)); % any value 
        Vleft = x(1)./(1+k.*x(3));
        Vright = x(2)./(1+k.*x(4));
    end

 dV =Vright-Vleft;
 gx = 1/(1+exp(-beta*dV+bias));


end