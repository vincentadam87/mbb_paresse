function [ gx ] = g_paresse_STABLE_BAYES(x,P,u,in)
% Model decides based on a softmax on expected values.
% Expected values are computed under the assumption that contingencies are
% stable.
% INPUT
% - x : Irrelevant
% - P : parameters
%     1 : inverse temperature in softmax
%     2 : effort discount factor left
%     3 : bias
% - u : input
%     1 total right
%     2 total left
%     3 win reward right
%     4 win reward left
%     5 win effort right
%     6 win effort left

beta = exp(P(1)); % inverse temperature in softmax positive
bias = P(3); % BIAS                   
counts = u;
     
    % ---- observation params P(4:11) sets prior beliefs for rewards and efforts
    % per hand. 
    alpha_p_L_E =(exp(P(4)))+1 + counts(6);
    beta_p_L_E = (exp(P(5)))+1 + counts(2)-counts(6);
    p_LEH = alpha_p_L_E/(alpha_p_L_E+beta_p_L_E); % posterior for small efforts given left choices 
    
    alpha_p_L_R =(exp(P(6)))+1 + counts(4);
    beta_p_L_R = (exp(P(7)))+1 + counts(2)-counts(4);
    p_LRH = alpha_p_L_R/(alpha_p_L_R+beta_p_L_R);
    
    alpha_p_R_E =(exp(P(8)))+1 + counts(5);
    beta_p_R_E = (exp(P(9)))+1 + counts(1)-counts(5);
    p_REH = alpha_p_R_E/(alpha_p_R_E+beta_p_R_E);
    
    alpha_p_R_R =(exp(P(10)))+1 + counts(3);
    beta_p_R_R =(exp(P(11)))+1 +counts(1)-counts(3);
    p_RRH = alpha_p_R_R/(alpha_p_R_R+beta_p_R_R);
    
    
    if isequal(in.DECIDING,'LINEAR')
        
        d = (P(2)); % any value for linear discount  
        Vleft = p_LRH-d*p_LEH; % expected value
        Vright = p_RRH-d*p_REH;
        
    elseif isequal(in.DECIDING,'HYPERBOLIC')
        
        k = (P(2)); % any value hyp discount
        Vleft = p_LRH*(1-p_LEH*k/(1+k)); % expected value
        Vright = p_RRH*(1-p_REH*k/(1+k));
    end

 dV =Vright-Vleft;
 gx = 1/(1+exp(-beta*dV+bias));


end