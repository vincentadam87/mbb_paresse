function [ gx ] = g_paresse_mixed(x,P,u,in)
% INPUT
% - x : Qvalues
%     1 : QRleft
%     2 : QRright
%     3 : QEleft
%     4 : QEright
% - P : parameters
%     1 : inverse temperature in softmax
%     2 : effort discount factor left
% - u : input
%     1 : choice (left = 0, right = 1)
%     2 : reward
%     3 : effort

beta = exp(P(1)); % inverse temperature in softmax positive

% bias = P(3); % only used for the BIAS models otherwise set to N(0,0). It changes for the BIAS_COND where there is bias per condition 
                   
counts = u(4:9); 
    % 1 total right
    % 2 total left
    % 3 win reward right
    % 4 win reward left
    % 5 win effort right
    % 6 win effort left
    
    % ---- observation params P(4:11) sets prior beliefs for rewards and efforts
    % per hand. Only used in the Ideal Bayes. 

% if isequal(in.TYPE_REW,'IDEAL') % && ~isequal(in.TYPE_EFF,'IDEAL') 
%    
%     alpha_p_L_R =(exp(P(3)))+1 + counts(4);
%     beta_p_L_R = (exp(P(4)))+1 + counts(2)-counts(4);
%     p_LRH = alpha_p_L_R/(alpha_p_L_R+beta_p_L_R);
%    
%     alpha_p_R_R =(exp(P(5)))+1 + counts(3);
%     beta_p_R_R =(exp(P(6)))+1+ counts(1)-counts(3);
%     p_RRH = alpha_p_R_R/(alpha_p_R_R+beta_p_R_R);
%     
%     reward_R=p_RRH;
%     reward_L=p_LRH; 
%     
% elseif isequal(in.TYPE_EFF,'IDEAL') % && ~isequal(in.TYPE_REW,'IDEAL') % one if the learning rules is Bayesian 
%    
%     alpha_p_L_E =(exp(P(3)))+1+counts(6);
%     beta_p_L_E = (exp(P(4)))+1+counts(2)-counts(6);
%     p_LEH = alpha_p_L_E/(alpha_p_L_E+beta_p_L_E); % posterior for small efforts given left choices 
%     
%     alpha_p_R_E =(exp(P(5)))+1+counts(5);
%     beta_p_R_E = (exp(P(6)))+1+counts(1)-counts(5);
%     p_REH = alpha_p_R_E/(alpha_p_R_E+beta_p_R_E);
%     
%     effort_R=p_REH;
%     effort_L=p_LEH; 
    
if isequal(in.TYPE_REW,'IDEAL') && isequal(in.TYPE_EFF,'IDEAL') % FULL BIAS 
    alpha_p_L_E =(exp(P(3)))+1 + counts(6);
    beta_p_L_E = (exp(P(4)))+1+counts(2)-counts(6);
    p_LEH = alpha_p_L_E/(alpha_p_L_E+beta_p_L_E); % posterior for small efforts given left choices 
    
    alpha_p_L_R =(exp(P(5)))+1 + counts(4);
    beta_p_L_R = (exp(P(6)))+1 + counts(2)-counts(4);
    p_LRH = alpha_p_L_R/(alpha_p_L_R+beta_p_L_R);
    
    alpha_p_R_E =(exp(P(7)))+1 + counts(5);
    beta_p_R_E = (exp(P(8)))+1+counts(1)-counts(5);
    p_REH = alpha_p_R_E/(alpha_p_R_E+beta_p_R_E);
    
    alpha_p_R_R =(exp(P(9)))+1 + counts(3);
    beta_p_R_R =(exp(P(10)))+1+ counts(1)-counts(3);
    p_RRH = alpha_p_R_R/(alpha_p_R_R+beta_p_R_R);
    
    reward_R=p_RRH;
    reward_L=p_LRH; 
    effort_R=p_REH;
    effort_L=p_LEH; 
    
end
    



if isequal(in.TYPE_REW,'QL') || isequal(in.TYPE_REW,'WSLS')
    reward_R = x(2); % the estimated Q-values taken from the evolution function f
    reward_L = x(1);
    
    if isequal(in.TYPE_EFF,'QL') || isequal(in.TYPE_EFF,'WSLS')
        
        effort_R = x(4); % the estimated Q-values taken from the evolution function f
        effort_L = x(3);
        
    elseif isequal(in.TYPE_EFF,'IDEAL')
        
        alpha_p_L_E =(exp(P(3)))+1+counts(6);
        beta_p_L_E = (exp(P(4)))+1+counts(2)-counts(6);
        p_LEH = alpha_p_L_E/(alpha_p_L_E+beta_p_L_E); % posterior for small efforts given left choices
        
        alpha_p_R_E =(exp(P(5)))+1+counts(5);
        beta_p_R_E = (exp(P(6)))+1+counts(1)-counts(5);
        p_REH = alpha_p_R_E/(alpha_p_R_E+beta_p_R_E);
        
        effort_R=p_REH;
        effort_L=p_LEH;
        
        
    end
        
elseif isequal(in.TYPE_EFF,'QL') || isequal(in.TYPE_EFF,'WSLS') 
    effort_R = x(4); % the estimated Q-values taken from the evolution function f
    effort_L = x(3);
    
    if isequal(in.TYPE_REW,'QL') || isequal(in.TYPE_REW,'WSLS')
        
        reward_R = x(2); % the estimated Q-values taken from the evolution function f
        reward_L = x(1);
       
        
    elseif isequal(in.TYPE_REW,'IDEAL')
        
        alpha_p_L_R =(exp(P(3)))+1 + counts(4);
        beta_p_L_R = (exp(P(4)))+1 + counts(2)-counts(4);
        p_LRH = alpha_p_L_R/(alpha_p_L_R+beta_p_L_R);
        
        alpha_p_R_R =(exp(P(5)))+1 + counts(3);
        beta_p_R_R =(exp(P(6)))+1+ counts(1)-counts(3);
        p_RRH = alpha_p_R_R/(alpha_p_R_R+beta_p_R_R);
        
        reward_R=p_RRH;
        reward_L=p_LRH;
        
    end
    
    
end
    


% the discoutn function 
    if isequal(in.DECIDING,'LINEAR')
        
        d =(P(2)); % any value for linear discount
        
        Vleft = reward_L-d.*effort_L;
        Vright =reward_R-d.*effort_R;
        
%         Vleft = p_LRH-d*p_LEH;
%         Vright = p_RRH-d*p_REH;
        
    elseif isequal(in.DECIDING,'HYPERBOLIC') % HOW IS IT COMPUTED A MIXED between VALUES ET probabilities? 
        
        k = (P(2)); % any value hyp discount
        
        Vleft = reward_L./(1+effort_L.*k); % THIS IS INCORRECT FOR THE COMBINATION WITH THE IDEAL OBSERVER 
        Vright =reward_R./(1+effort_R.*k); 
        
%         Vleft = p_LRH*(1-p_LEH*k/(1+k));
%         Vright = p_RRH*(1-p_REH*k/(1+k));
    end


% elseif isequal(in.TYPE,'BIAS_COND')
%   
%      Vleft=0; 
%      Vright=0; 
%      
%      
% elseif isequal(in.TYPE,'BIAS_HAND') % BIAS per CONDITION or BIAS per hand  
%   
%      Vleft=0; 
%      Vright=0; 
%      
%         
% else  
%     
    
%     if isequal(in.DECIDING,'LINEAR')
%         d = (P(2)); % any value
%         Vleft = x(1)-d*x(3);
%         Vright = x(2)-d*x(4);
%     elseif isequal(in.DECIDING,'HYPERBOLIC')
%         k = (P(2)); % any value 
%         Vleft = x(1)./(1+k.*x(3));
%         Vright = x(2)./(1+k.*x(4));
%     end


    
    % value = sum of the best outcomes / total number of choices ... 
    
%end

% the SOFTMAX 

 dV =Vright-Vleft;
 gx = 1/(1+exp(-beta*dV));


end