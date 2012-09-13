function [ gx ] = g_paresse_MATHYS(x,P,u,in)
% INPUT
% - x : Qvalues
%     1:5 : pleft_R
%     6:10 : pright_R
%     11:15 : pleft_E
%     16:20 : pright_E
% - P : parameters
%     1 : inverse temperature in softmax
%     2 : effort discount factor left
%     3 : bias
% - u : input
%     1 : choice (left = 0, right = 1)
%     2 : reward
%     3 : effort

beta = exp(P(1)); % inverse temperature in softmax positive
bias = P(3); % only used for the BIAS models otherwise set to N(0,0). It changes for the BIAS_COND where there is bias per condition


if isequal(in.DECIDING,'LINEAR')
    
    
    % Approximation first order
    pL_R = sig( x(2));
    pR_R = sig( x(12));
    pL_E = sig( x(7));
    pR_E = sig( x(17));
    d = (P(2)); % any value
    
    % fixed form approximation of the expected probability
    a = 0.368;
    pL_R = sig( x(2)./sqrt(1+a*x(3)));
    pR_R = sig( x(12)./sqrt(1+a*x(13)));
    pL_E = sig( x(7)./sqrt(1+a*x(8)));
    pR_E = sig( x(17)./sqrt(1+a*x(18)));
    
    Vleft = pL_R-d*pL_E;
    Vright = pR_R-d*pR_E;
    
    
    
elseif isequal(in.DECIDING,'HYPERBOLIC')
    % Approximation first order
    pL_R = sig( x(2));
    pR_R = sig( x(12));
    pL_E = sig( x(7));
    pR_E = sig( x(17));
    d = (P(2)); % any value
    
    % fixed form approximation of the expected probability
    a = 0.368;
    pL_R = sig( x(2)./sqrt(1+a*x(3)));
    pR_R = sig( x(12)./sqrt(1+a*x(13)));
    pL_E = sig( x(7)./sqrt(1+a*x(8)));
    pR_E = sig( x(17)./sqrt(1+a*x(18)));
    
    
    k = (P(2)); % any value
    Vleft = pL_R./(1+k.*pL_E);
    Vright = pR_R./(1+k.*pR_E);
    
    
end



% value = sum of the best outcomes / total number of choices ...



dV =Vright-Vleft;
gx = 1/(1+exp(-beta*dV+bias));


end


function y=sig(x)
y = 1/(1+exp(-x));
y(y<eps) = eps;
y(y>1-eps) = 1-eps;
end