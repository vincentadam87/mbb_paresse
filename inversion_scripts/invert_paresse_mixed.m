function [posterior,out] = invert_paresse_mixed(i_subject,SESS,COND,TYPE_REW,TYPE_EFF,DECIDING,UPDATE)
% This script performs the inversion of the task 'paresse' for given
% - subject
% - sessions & conditions
% - model specifications


% INPUT
% - i_subject : index of the subject
% - SESS : indices of the sessions considered (1 to 3)
% - COND : indices of the conditions considered (1 to 4 for A to D)
% - TYPE_REW: - learning rule used for the reward learning 
%   TYPE_EFF: learning rule used for E learning
% - DECIDING : linear or hyperbolic discount 
% - UPDATE : shared or different: usedn only for the FULL QL model 

% OUTPUT : standard inversion output
% - posterior
% - out

[Y,U,IsYout,cond] = load_data_paresse(i_subject,SESS,COND);

Ntrials = size(Y,2);

try UPDATE; % for QL models only 
catch e;
    UPDATE = [];
end
try TYPE_REW; 
catch e;
    TYPE_REW = [];
end

try TYPE_EFF; 
catch e;
    TYPE_EFF = [];
end

try DECIDING;
catch e;
    DECIDING = [];
end

%----------------------------------------------------------
% Specifying models

% ------------------------------------
% -- Declare initial model

dim_output = 1; % choices 
dim_data = 9; % choice/reward/effort/total right choices/total left choices/wins R|right/wins R|left/wins E|right/wins E|left 

% structure for single session 

dim = struct('n',4,... % 4Q hidden states 
    'p',dim_output,... % total output dimension
    'n_theta',2,... % evolution parameters 1 learning rate 
    'n_phi',2,... % observation parameters beta/discount/bias 
    'u',dim_data,... % data dimension per session
    'n_t',Ntrials); %


options.DisplayWin = 0; % display inversion
options.GnFigs = 0;
options.binomial = 1; % Dealing with binary data
options.isYout = zeros(1,Ntrials); % Excluding data points
options.dim = dim;
options.inG.DECIDING =DECIDING; % optional information used in the observation and evomution functions 
options.inG.TYPE_REW = TYPE_REW;
options.inG.TYPE_EFF = TYPE_EFF;
options.inF.UPDATE = UPDATE;
options.inF.TYPE_REW = TYPE_REW;
options.inF.TYPE_EFF = TYPE_EFF;

% Feedback necessity -- optional to calculate the PE after 
% options.skipf = zeros(1,Ntrials);
% options.skipf(1) = 1;

U = [zeros(size(U,1),1),U(:,1:end-1)];

% ------------------------------------
% -- Extend to multiple sessions
Nsessions = length(SESS)*length(COND); % 3 x 4 = 12 sessions per subject 

% future dimensions of the extended model

if (isequal(TYPE_REW,'IDEAL') && ~isequal(TYPE_EFF,'IDEAL')) || (isequal(TYPE_EFF,'IDEAL') && ~isequal(TYPE_REW,'IDEAL'))
    dim_e = struct('n_theta',2,... % specify parameter space
        'n_phi',6); % 4 params to update priors + discount + beta
    
elseif isequal(TYPE_REW,'IDEAL') && isequal(TYPE_EFF,'IDEAL')
    dim_e = struct('n_theta',2,... % specify parameter space
        'n_phi',10); % 8 params to update priors + beta + discount
    
    % elseif isequal(TYPE,'BIAS_COND')
    %      dim_e = struct('n_theta',1,... %
    %         'n_phi',4*3); % will keep the three params because bias = P(3) in g_paresse_last script for consistency
else
    dim_e = struct('n_theta',2,... % specify parameter space
        'n_phi',2);% beta/discount/
end

           
% Information for sessions

in_sessions = struct();
in_sessions.n_sess = Nsessions; % number of sessions
in_sessions.ind.theta = ones(Nsessions,1)*[1:2]; % specify parameter use for each session
in_sessions.ind.phi =ones(Nsessions,1)*[1:2]; % specify parameter use for each session


if isequal(TYPE_REW,'QL') && isequal(TYPE_EFF,'QL') % learning rate is only in the Q models 
% 
% % ----specify evolution params for each session/model whether you have shared or different learning rates within the Q models 
% 
if  isequal(UPDATE,'SHARED')
    in_sessions.ind.theta = ones(Nsessions,1)*[1,1]; % one learning rate for R and E 
elseif isequal(UPDATE,'DIFFERENT')
    in_sessions.ind.theta = ones(Nsessions,1)*[1:2]; % two learning rates 
end

end

in_sessions.f_fname =@f_paresse_mixed; % handle of the shared evolution function
in_sessions.g_fname =@g_paresse_mixed; % handle of the shared observation function
in_sessions.dim_e = dim_e; % specify extended model's parameter space
in_sessions.binomial = 1;
in_sessions.DisplayWin = 1;

% ---- optional information for evolution and observation functions 

in_sessions.inF.TYPE_REW =TYPE_REW;
in_sessions.inF.TYPE_EFF =TYPE_EFF;
in_sessions.inF.UPDATE =UPDATE;
in_sessions.inG.DECIDING=DECIDING;
in_sessions.inG.TYPE_REW=TYPE_REW;
in_sessions.inG.TYPE_EFF=TYPE_EFF;


% ----- specify parameter use for each session
if (isequal(TYPE_REW,'IDEAL') && ~isequal(TYPE_EFF,'IDEAL')) || (isequal(TYPE_EFF,'IDEAL') && ~isequal(TYPE_REW,'IDEAL'))
    in_sessions.ind.phi = ones(Nsessions,1)*(1:6); % 4 hyperparameters + beta temperature 
elseif isequal(TYPE_REW,'IDEAL') && isequal(TYPE_EFF,'IDEAL')
    in_sessions.ind.phi = ones(Nsessions,1)*(1:10); % 8 hyperparameters + beta softmax temperature 
end

% elseif isequal(TYPE,'BIAS_COND')
% in_sessions.ind.phi = in_sessions.ind.phi +3*(cond-1)*ones(1,3); % 3*(cond-1)*ones(1,3);
%      
% end

% ------ MAKE EXTENDED MODEL ACROSS SESSIONS 

[f_fname_e,g_fname_e,dim_e,options_e] = makeExtendedModel(dim,options,in_sessions);



%----- Defining Priors on parameters (mean and Covariance matrix)
% ----- Evolution Params 
priors.muPhi = zeros(dim_e.n_phi,1);
priors.SigmaPhi = 1*eye(dim_e.n_phi);

priors.muTheta = zeros(dim_e.n_theta,1);
priors.SigmaTheta = 1*eye(dim_e.n_theta);


if isequal(TYPE_REW,'WSLS') && ~isequal(TYPE_EFF,'WSLS')
   
    priors.muTheta(1,1) = 1;
    priors.SigmaTheta(1,1)=0; %  1*eye(dim_e.n_theta); % set variance to 10. 
    % priors.SigmaTheta(1,1) = 1*eye(dim_e.n_theta);
    priors.muPhi(1) = 0; % beta
    priors.SigmaPhi(1,1) = 1;
    
    priors.muPhi(2) = 0; % discount
    priors.SigmaPhi(2,2) = 1; %
    
elseif isequal(TYPE_EFF,'WSLS') && ~isequal(TYPE_REW,'WSLS')
    
    priors.muTheta(2,1) = 1; % zeros(dim_e.n_theta,1);
    priors.SigmaTheta(2,2) =0; %  1*eye(dim_e.n_theta); % set variance to 10. 
    
    priors.muPhi(1) = 0; % beta
    priors.SigmaPhi(1,1) = 1;
    
    priors.muPhi(2) = 0; % discount
    priors.SigmaPhi(2,2) = 1; %
    
elseif isequal(TYPE_EFF,'WSLS') && isequal(TYPE_REW,'WSLS') % both learnng rates at 1 with 0 valueance 
    priors.muTheta = ones(dim_e.n_theta,1); % set LR to 1 
    priors.SigmaTheta = 0*eye(dim_e.n_theta); % set variance to 0. 
    
    priors.muPhi(1) = 0; % beta
    priors.SigmaPhi(1,1) = 1;
    
    priors.muPhi(2) = 0; % discount
    priors.SigmaPhi(2,2) = 1; %
    
elseif isequal(TYPE_REW,'QL') && ~isequal(TYPE_EFF,'QL')
   
    priors.muTheta(1,1) = 0;
    priors.SigmaTheta(1,1)=1; %  1*eye(dim_e.n_theta); % set variance to 10. 
    % priors.SigmaTheta(1,1) = 1*eye(dim_e.n_theta);
    
    priors.muPhi(1) = 0; % beta
    priors.SigmaPhi(1,1) = 1;
    
    priors.muPhi(2) = 0; % discount
    priors.SigmaPhi(2,2) = 1; %
    
elseif isequal(TYPE_EFF,'QL') && ~isequal(TYPE_REW,'QL')
    
    priors.muTheta(2,1) = 0; % zeros(dim_e.n_theta,1);
    priors.SigmaTheta(2,2) =1; %  1*eye(dim_e.n_theta); % set variance to 10. 
    
    priors.muPhi(1) = 0; % beta
    priors.SigmaPhi(1,1) = 1;
    
    priors.muPhi(2) = 0; % discount
    priors.SigmaPhi(2,2) = 1; %
    
elseif isequal(TYPE_EFF,'QL') && isequal(TYPE_REW,'QL') % both 
    priors.muTheta = zeros(dim_e.n_theta,1);
    priors.SigmaTheta =1*eye(dim_e.n_theta); 
    
    priors.muPhi = zeros(dim_e.n_theta,1); % beta
    priors.SigmaPhi=1*eye(dim_e.n_theta); 
    
    if isequal(UPDATE,'SHARED') % one learning rate
        priors.SigmaTheta(2,2)=0; % do not infer 2nd learning rate if shared
    end

    
elseif isequal(TYPE_REW,'IDEAL') && ~isequal(TYPE_EFF,'IDEAL') % both learning rates at 0 xith no variance 
       priors.muTheta(1,1) = 0;
       priors.SigmaTheta(1,1)=0;
       
       priors.muPhi(1) = 0; % beta
       priors.SigmaPhi(1,1) = 1;
       
       priors.muPhi(2) = 0; % discount
       priors.SigmaPhi(2,2) = 1; %
       
       priors.muPhi(3) = 0; % prior on low discount
       priors.SigmaPhi(3,3) = 1; % NO BIAS prior on low discount
       
       priors.muPhi(4) = 0; % prior on low discount
       priors.SigmaPhi(4,4) = 1; % NO BIAS prior on low discount
       
       priors.muPhi(5) = 0; % prior on low discount
       priors.SigmaPhi(5,5) = 1; % NO BIAS prior on low discount
       
       priors.muPhi(6) = 0; % prior on low discount
       priors.SigmaPhi(6,6) = 1; % NO BIAS prior on low discount
       
elseif isequal(TYPE_EFF,'IDEAL') && ~isequal(TYPE_REW,'IDEAL')
       priors.muTheta(2,1) = 0;
       priors.SigmaTheta(2,2)=0;
       
       priors.muPhi(1) = 0; % beta
       priors.SigmaPhi(1,1) = 1;
       
       priors.muPhi(2) = 0; % discount
       priors.SigmaPhi(2,2) = 1; %
       
       priors.muPhi(3) = 0; % prior on low discount
       priors.SigmaPhi(3,3) = 1; % NO BIAS prior on low discount
       
       priors.muPhi(4) = 0; % prior on low discount
       priors.SigmaPhi(4,4) = 1; % NO BIAS prior on low discount
       
       priors.muPhi(5) = 0; % prior on low discount
       priors.SigmaPhi(5,5) = 1; % NO BIAS prior on low discount
       
       priors.muPhi(6) = 0; % prior on low discount
       priors.SigmaPhi(6,6) = 1; % NO BIAS prior on low discount
       
elseif isequal(TYPE_REW,'IDEAL') && isequal(TYPE_EFF,'IDEAL')
        priors.muTheta = zeros(dim_e.n_theta,1);
        priors.SigmaTheta = 0*eye(dim_e.n_theta);
        
       priors.muPhi(1) = 0; % beta 
       priors.SigmaPhi(1,1) = 1;
       
       priors.muPhi(2) = 0; % discount
       priors.SigmaPhi(2,2) = 1; %
       
       priors.muPhi(3) = 0; %priors on beleifs about proportions of wins and failures for each hand / each dimension
       priors.SigmaPhi(3,3) = 1; % NO BIAS prior on low discount
       
       priors.muPhi(4) = 0;
       priors.SigmaPhi(4,4) = 1;
       
       priors.muPhi(5) = 0; %
       priors.SigmaPhi(5,5) = 1; %
       
       priors.muPhi(6) = 0; %
       priors.SigmaPhi(6,6) = 1; %
       
       priors.muPhi(7) = 0; %
       priors.SigmaPhi(7,7) = 1; %
       
       priors.muPhi(8) = 0; %
       priors.SigmaPhi(8,8) = 1; %
       
       priors.muPhi(9) = 0; %
       priors.SigmaPhi(9,9) = 1; %
       
       priors.muPhi(10) = 0; %
       priors.SigmaPhi(10,10) = 1; %

end

%     if isequal(UPDATE,'SHARED') % one learning rate
%         priors.SigmaTheta(2,2)=0; % do not infer 2nd learning rate if shared
%     end


% ---- Observation params 

% priors.muPhi(1) = 0; % beta
% priors.SigmaPhi(1,1) = 1;
% 
% priors.muPhi(2) = 0; % discount 
% priors.SigmaPhi(2,2) = 1; % 

% if isequal(TYPE,'BIAS_COND')
%     
%     priors.muPhi(1) = 0; % beta
%     priors.SigmaPhi(1,1) = 0; %
%     
%     priors.muPhi(2) = 0; % discount
%     priors.SigmaPhi(2,2) = 0; %
%     
%      priors.muPhi(3) = 0; % bias
%     priors.SigmaPhi(3,3) = 1; %
    
% if (isequal(TYPE_REW,'IDEAL') && ~isequal(TYPE_EFF,'IDEAL')) || (~isequal(TYPE_REW,'IDEAL') && isequal(TYPE_EFF,'IDEAL'))
%     priors.muPhi(1) = 0; % beta 
%     priors.SigmaPhi(1,1) = 0; %
%     
%     priors.muPhi(2) = 0; % discount 
%     priors.SigmaPhi(2,2) = 1; % prior on low discount
%     
%     priors.muPhi(3) = 0; % prior on low discount
%     priors.SigmaPhi(3,3) = 1; % NO BIAS prior on low discount
% 
%     priors.muPhi(4) = 0; % prior on low discount
%     priors.SigmaPhi(4,4) = 1; % NO BIAS prior on low discount
%     
%     priors.muPhi(5) = 0; % prior on low discount
%     priors.SigmaPhi(5,5) = 1; % NO BIAS prior on low discount
%     
%     priors.muPhi(6) = 0; % prior on low discount
%     priors.SigmaPhi(6,6) = 1; % NO BIAS prior on low discount
%     
% elseif (isequal(TYPE_REW,'IDEAL') && isequal(TYPE_EFF,'IDEAL')) % FULL Ideal model 
%     priors.muPhi(1) = 0; % beta 
%     priors.SigmaPhi(1,1) = 0; %
%     
%     priors.muPhi(2) = 0; % discount 
%     priors.SigmaPhi(2,2) = 1; % prior on low discount
%     
%     priors.muPhi(3) = 0; %priors on beleifs about proportions of wins and failures for each hand / each dimension 
%     priors.SigmaPhi(3,3) = 1; % NO BIAS prior on low discount
% 
%     priors.muPhi(4) = 0;
%     priors.SigmaPhi(4,4) = 1; 
%     
%     priors.muPhi(5) = 0; % 
%     priors.SigmaPhi(5,5) = 1; %
%     
%     priors.muPhi(6) = 0; % 
%     priors.SigmaPhi(6,6) = 1; % 
%     
%     priors.muPhi(7) = 0; % 
%     priors.SigmaPhi(7,7) = 1; % 
%     
%     priors.muPhi(8) = 0; % 
%     priors.SigmaPhi(8,8) = 1; % 
%     
%     priors.muPhi(9) = 0; % 
%     priors.SigmaPhi(9,9) = 1; % 
%     
%     priors.muPhi(10) = 0; % 
%     priors.SigmaPhi(10,10) = 1; % 
% end

% if isequal(BIAS,'BIAS OFF')
% priors.SigmaPhi(3,3) = 0; % NO BIAS in the model default is BIAS ON 
% end

% Priors on initial hidden states (values) or initial beleifs of rewards
% and effort frequencies for the right and left hands 
x0 = [0.5;0.5;0.5;0.5];

priors.muX0 = repmat(x0,Nsessions,1);
priors.SigmaX0 = 0*eye(dim_e.n); % NO VARIANCE add variance on the priors for initial hidden states 

% No state noise for deterministic update rules
priors.a_alpha = Inf;
priors.b_alpha = 0;

options_e.priors = priors;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MODEL INVERSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[posterior,out] = VBA_NLStateSpaceModel(Y,U,f_fname_e,g_fname_e,dim_e,options_e);
