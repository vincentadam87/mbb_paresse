function [posterior,out] = invert_paresse_STABLE_BAYES(i_subject,SESS,COND,DECIDING,BIAS,DISPLAY)
% This script performs the inversion of the task 'paresse' for given
% - subject
% - sessions & conditions
% - model specifications
% DECIDING :- hyperbolic or linear integration in the softmax 
% UPDATE - shared or different learning rates - for the softmax only.
% BIAS - whether there is bias in learning for the hand; default is bias on, 

% INPUT
% - i_subject : index of the subject
% - SESS : indices of the sessions considered (1 to 3)
% - COND : indices of the conditions considered (1 to 4 for A to D)
% - TYPE : model type 
% - DECIDING : linear or hyperbolic discount 

% OUTPUT : standard inversion output
% - posterior
% - out

try DISPLAY;
catch e; DISPLAY=0;
end

%-------------------------------------------
% TEST :
if i_subject == 0
    i_subject = 7;
   SESS= [1,2,3];
   COND=[1,2,3,4];
   DECIDING = 'LINEAR';
   BIAS = 'BIAS OFF';
end
%-------------------------------------------

[Y,O,IsYout,cond,sess,U] = load_data_paresse(i_subject,SESS,COND);

Ntrials = size(Y,2);


% Check arguments
try DECIDING;catch e; DECIDING = [];end
try BIAS;    catch e; BIAS =[];end
%----------------------------------------------------------
% Specifying models

% ------------------------------------
% -- Declare initial model (for single session)

dim_output = 1; % choices 
dim_data = 6; % 
% - choice
% - reward
% - effort
% - total right choices
% - total left choices
% - wins R|right ( winning = getting best output out of 2 )
% - wins R|left
% - wins E|right
% - wins E|left 

% structure for single session 
dim = struct('n',0,... % no hidden states 
    'p',dim_output,... % total output dimension
    'n_theta',0,... % evolution parameters 2 learning rates 
    'n_phi',11,... % observation parameters beta/discount/bias 
    'u',dim_data,... % data dimension per session
    'n_t',Ntrials); %


options.DisplayWin = DISPLAY; % display inversion
options.GnFigs = 0;
options.binomial = 1; % Dealing with binary data
options.isYout = zeros(1,Ntrials); % Excluding data points
options.dim = dim;

% Feedback necessity -- optional to calculate the PE after 
options.skipf = zeros(1,Ntrials);
options.skipf(1) = 1;
U = [zeros(size(U,1),1),U(:,1:end-1)]; % shift data U one to the right [u1,u2,...,uN]-> [~,u1,u2,...,uN-1]

% ------------------------------------
% -- Extend to multiple sessions
Nsessions = length(SESS)*length(COND); % 3 x 4 = 12 sessions per subject 

% -- Dimensions of the extended model to account for multiple sessions

    dim_e = struct('n_theta',0,... % specify parameter space
        'n_phi',11); % 1+8 params to update priors

           
% Information for sessions

in_sessions = struct();
in_sessions.n_sess = Nsessions; % number of sessions

in_sessions.f_fname = []; % handle of the shared evolution function
in_sessions.g_fname = @g_paresse_STABLE_BAYES_last; % handle of the shared observation function

in_sessions.dim_e = dim_e; % specify extended model's parameter space
in_sessions.binomial = 1;
in_sessions.DisplayWin = DISPLAY;

% ---- optional information for evolution and observation functions 

in_sessions.inG.DECIDING=DECIDING;
in_sessions.inF = [];

% ----- specify parameter use for each session
in_sessions.ind.phi = ones(Nsessions,1)*(1:11);
in_sessions.ind.theta = [];


[ f_fname_e,g_fname_e,dim_e,options_e ] = makeExtendedModel(dim,options,in_sessions);



%----- Defining Priors on parameters (mean and Covariance matrix)
% ----- Evolution Params 
priors.muPhi = zeros(dim_e.n_phi,1);
priors.SigmaPhi = 1*eye(dim_e.n_phi);

priors.muTheta = zeros(dim_e.n_theta,1); % no     
priors.SigmaTheta =0*eye(dim_e.n_theta); % no  

% ---- Observation params 


    
    priors.muPhi(1) = 0; % to avoid numbers < 1 CHECK !!! 
    priors.SigmaPhi(1,1) = 1; %
    
    priors.muPhi(2) = 0; % NO discount 
    priors.SigmaPhi(2,2) = 1; % prior on low discount
    
 if isequal(BIAS,'BIAS OFF') % default is BIAS ON
    priors.muPhi(3) = 0; %
    priors.SigmaPhi(3,3) = 0; %
else
    priors.muPhi(3) = 0; %
    priors.SigmaPhi(3,3) = 1; %
end


    priors.muPhi(4) = 0; % prior on low discount
    priors.SigmaPhi(4,4) = 1; % NO BIAS prior on low discount
    
       priors.muPhi(5) = 0; % prior on low discount
    priors.SigmaPhi(5,5) = 1; % NO BIAS prior on low discount
    
    priors.muPhi(6) = 0; % prior on low discount
    priors.SigmaPhi(6,6) = 1; % NO BIAS prior on low discount
    
       priors.muPhi(7) = 0; % prior on low discount
    priors.SigmaPhi(7,7) = 1; % NO BIAS prior on low discount
    
       priors.muPhi(8) = 0; % prior on low discount
    priors.SigmaPhi(8,8) = 1; % NO BIAS prior on low discount

     priors.muPhi(9) = 0; % prior on low discount
    priors.SigmaPhi(9,9) = 1; % NO BIAS prior on low discount
    
       priors.muPhi(10) = 0; % prior on low discount
    priors.SigmaPhi(10,10) = 1; % NO BIAS prior on low discount

       priors.muPhi(11) = 0; % prior on low discount
    priors.SigmaPhi(11,11) = 1; % NO BIAS prior on low discount

% Priors on initial hidden states (values) or initial beleifs of rewards
% and effort frequencies for the right and left hands 

priors.muX0 = [];
priors.SigmaX0 = []; % NO



% No state noise for deterministic update rules
priors.a_alpha = Inf;
priors.b_alpha = 0;

options_e.priors = priors;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MODEL INVERSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[posterior,out] = VBA_NLStateSpaceModel(Y,U,f_fname_e,g_fname_e,dim_e,options_e);
