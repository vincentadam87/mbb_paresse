function [posterior,out] = invert_paresse_MATHYS(i_subject,SESS,COND,DECIDING,LEARNING,VOLATILITY,BIAS,DISPLAY)
% This script performs the inversion of the task 'paresse' for given
% - subject
% - sessions & conditions
% - model specifications
% TYPE : - the model type we are using: QL, WSLS, IDEAL BAYESIAN, CONDITION BIAS or SIMPLE BIAS.
% DECIDING :- hyperbolic or linear integration in the softmax
% UPDATE - shared or different learning rates - for the softmax only.
% BIAS - whether there is bias in learning for the hand; default is bias on,

% INPUT
% - i_subject : index of the subject
% - SESS : indices of the sessions considered (1 to 3)
% - COND : indices of the conditions considered (1 to 4 for A to D)
% - TYPE : model type
% - DECIDING : linear or hyperbolic discount
% - UPDATE : 'SHARED' or 'DIFFERENT', whether effort and reward share the
% same learning rate

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
    LEARNING = 'DIFFERENT' ;
end
%-------------------------------------------


[Y,U,IsYout,cond] = load_data_paresse(i_subject,SESS,COND);

Ntrials = size(Y,2);


% Check arguments
try LEARNING;  catch e; LEARNING = [];end  % for QL models only
try TYPE;    catch e; TYPE = [];end
try DECIDING;catch e; DECIDING = [];end
try BIAS;    catch e; BIAS =[];end
try VOLATILITY;  catch e; VOLATILITY = [];end  % for QL models only


learn_volatility = isequal(VOLATILITY,'LEARN VOLATILITY');

%----------------------------------------------------------
% Specifying models

% ------------------------------------
% -- Declare initial model (for single session)

dim_output = 1; % choices
dim_data = 3; %
% - choice (binary)
% - reward (binary)
% - effort (binary)

dim=struct('n',4*5,... % 4*5 (4 probabilities)*(5 hidden states per probability to learn)
    'p',dim_output,... % total output dimension
    'n_theta',3,... % evolution parameters(kappa/omega/theta)
    'n_phi',3,... % observation parameters beta/discount/bias
    'u',dim_data,... % data dimension per session
    'n_t',Ntrials);



options.DisplayWin = DISPLAY; % display inversion
options.GnFigs = 0;
options.binomial = 1; % Dealing with binary data
options.isYout = zeros(1,Ntrials); % Excluding data points
options.dim = dim;
options.inG.DECIDING =DECIDING; % optional information used in the observation and evomution functions
options.inG.TYPE = TYPE;
options.inF.LEARNING = LEARNING;
options.inF.TYPE = TYPE;

options.inF.lev2 = learn_volatility; % remove 3rd level (volatility learning)
options.inF.kaub = 1;%1.4; % upper bound for Kappa
options.inF.thub = 1; % upper bound for Theta
options.inF.rf = -1;

% Feedback necessity -- optional to calculate the PE after
options.skipf = zeros(1,Ntrials);
options.skipf(1) = 1;
U = [zeros(size(U,1),1),U(:,1:end-1)]; % shift data U one to the right [u1,u2,...,uN]-> [~,u1,u2,...,uN-1]

% ------------------------------------
% -- Extend to multiple sessions
Nsessions = length(SESS)*length(COND); % 3 x 4 = 12 sessions per subject

% -- Dimensions of the extended model to account for multiple sessions
dim_e = struct('n_theta',4,... %
               'n_phi',3); %

% Information for sessions
I= ones(Nsessions,1);

in_sessions = struct();
in_sessions.n_sess = Nsessions; % number of sessions
in_sessions.ind.theta = I*[1:dim_e.n_phi]; % specify parameter use for each session
in_sessions.ind.phi = I*[1:dim_e.n_phi]; % specify parameter use for each session

% Condition effect on learning : shared or different omega (= parameter 2)

if isequal(LEARNING,'SHARED')
in_sessions.ind.theta = I*[1,2,2,4]; % specify parameter use for each session
elseif isequal(LEARNING,'DIFFERENT')
in_sessions.ind.theta = I*[1,2,3,4]; % specify parameter use for each session
end

in_sessions.f_fname = @f_paresse_MATHYS_last; % handle of the shared evolution function
in_sessions.g_fname = @g_paresse_MATHYS_last; % handle of the shared observation function
in_sessions.dim_e = dim_e; % specify extended model's parameter space
in_sessions.binomial = 1;
in_sessions.DisplayWin = DISPLAY;

% ---- optional information for evolution and observation functions

in_sessions.inF = options.inF;
in_sessions.inG = options.inG;

% ----- specify parameter use for each session

[ f_fname_e,g_fname_e,dim_e,options_e ] = makeExtendedModel(dim,options,in_sessions);

%----- Defining Priors on parameters (mean and Covariance matrix)
% ----- Evolution Params
priors.muPhi = zeros(dim_e.n_phi,1);
priors.SigmaPhi = 1*eye(dim_e.n_phi);

priors.muTheta = zeros(dim_e.n_theta,1); % no learning rate
priors.SigmaTheta =1*eye(dim_e.n_theta); % no variance

% FIX KAPPA
priors.SigmaTheta(1,1) =0; % fix kappa

if isequal(LEARNING,'SHARED')
priors.SigmaTheta(3,3) =0; % do not infer 2nd omega
end

if ~isequal(VOLATILITY,'LEARN VOLATILITY')
priors.SigmaTheta(4,4) =0; % do not infer theta (the lastt parameter)
end
% ---- Observation params

priors.muPhi(1) = 0; % beta
priors.SigmaPhi(1,1) = 1;

priors.muPhi(2) = 0; % discount
priors.SigmaPhi(2,2) = 1; %




if isequal(BIAS,'BIAS OFF') % default is BIAS ON
    priors.muPhi(3) = 0; %
    priors.SigmaPhi(3,3) = 0; %
else
    priors.muPhi(3) = 0; %
    priors.SigmaPhi(3,3) = 1; %
end

% Priors on initial hidden states (values) or initial beleifs of rewards

x0 = repmat([0.5;0;1;1;1],4,1);
%x0 = [0.5;0;1;1;1];

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
