% Script to perform inversion for several models on the data-set of the
% task 'paresse' (mbb, 2011-2012, Vasilisa)


%%
%---------------------------------------------------------
% Performing the inversion for all subjects and all models
%---------------------------------------------------------

SESS = [1,2,3]; % 3 sessions per subject
COND = [1,2,3,4]; % A,B,C,D
SUBJECTS =[1:20]; % 20 subjects

Nsubjects = length(SUBJECTS);
DISPLAY = 0;
%
% COND = [1,2,3,4]; % A,B,C,D
% SESS = [3]; % 3 sessions per subject
% i_s =[9]; % 20 subjects

for i_s = SUBJECTS
    
    disp('--------------------------------------------------------')
    disp(['------Subject ',num2str(i_s),'/',num2str(length(SUBJECTS))])
    disp('--------------------------------------------------------')

    INVERSION{i_s}.i_subject = i_s;
 
    disp(['----- RL MODELS -----'])

    % RL MODELS
   [p,o] = invert_paresse_RL(i_s,SESS,COND,'QL','LINEAR','SHARED','BIAS OFF',DISPLAY);
   INVERSION{i_s}.QL.LIN.SHARED.NO_BIAS = struct('p',p,'o',o);
   [p,o] = invert_paresse_RL(i_s,SESS,COND,'QL','LINEAR','DIFFERENT','BIAS OFF',DISPLAY);  
   INVERSION{i_s}.QL.LIN.DIFFERENT.NO_BIAS = struct('p',p,'o',o);
   [p,o] = invert_paresse_RL(i_s,SESS,COND,'QL','HYPERBOLIC','SHARED','BIAS OFF',DISPLAY);
   INVERSION{i_s}.QL.HYP.SHARED.NO_BIAS = struct('p',p,'o',o);
   [p,o] = invert_paresse_RL(i_s,SESS,COND,'QL','HYPERBOLIC','DIFFERENT','BIAS OFF',DISPLAY);  
   INVERSION{i_s}.QL.HYP.DIFFERENT.NO_BIAS = struct('p',p,'o',o);
   
   [p,o] = invert_paresse_RL(i_s,SESS,COND,'QL','LINEAR','SHARED','BIAS ON',DISPLAY);
   INVERSION{i_s}.QL.LIN.SHARED.BIAS = struct('p',p,'o',o);
   [p,o] = invert_paresse_RL(i_s,SESS,COND,'QL','LINEAR','DIFFERENT','BIAS ON',DISPLAY);  
   INVERSION{i_s}.QL.LIN.DIFFERENT.BIAS = struct('p',p,'o',o);
   [p,o] = invert_paresse_RL(i_s,SESS,COND,'QL','HYPERBOLIC','SHARED','BIAS ON',DISPLAY);  
   INVERSION{i_s}.QL.HYP.SHARED.BIAS = struct('p',p,'o',o);
   [p,o] = invert_paresse_RL(i_s,SESS,COND,'QL','HYPERBOLIC','DIFFERENT','BIAS ON',DISPLAY); 
   INVERSION{i_s}.QL.HYP.DIFFERENT.BIAS = struct('p',p,'o',o);
   
   disp(['----- GREEDY MODELS -----'])
       
   % GREEDY MODEL (NO MEMORY) (always shared, both learning rates set to 1)
   [p,o] = invert_paresse_RL(i_s,SESS,COND,'WSLS','LINEAR','SHARED','BIAS OFF',DISPLAY);
   INVERSION{i_s}.GREEDY.LIN.SHARED.NO_BIAS = struct('p',p,'o',o);
   [p,o] = invert_paresse_RL(i_s,SESS,COND,'WSLS','HYPERBOLIC','SHARED','BIAS OFF',DISPLAY);
   INVERSION{i_s}.GREEDY.HYP.SHARED.NO_BIAS = struct('p',p,'o',o);
   [p,o] = invert_paresse_RL(i_s,SESS,COND,'WSLS','LINEAR','SHARED','BIAS ON',DISPLAY);
   INVERSION{i_s}.GREEDY.LIN.SHARED.BIAS = struct('p',p,'o',o);
   [p,o] = invert_paresse_RL(i_s,SESS,COND,'WSLS','HYPERBOLIC','SHARED','BIAS ON',DISPLAY);
   INVERSION{i_s}.GREEDY.HYP.SHARED.BIAS = struct('p',p,'o',o);

   
   disp(['----- STABLE BAYES MODELS -----'])

   % STABLE BAYES (STABILITY ASSUMPTION)
   invert_paresse_STABLE_BAYES(i_s,SESS,COND,'LINEAR','BIAS OFF',DISPLAY)
   INVERSION{i_s}.STABLE_BAYES.LIN.NO_BIAS = struct('p',p,'o',o);
   invert_paresse_STABLE_BAYES(i_s,SESS,COND,'LINEAR','BIAS ON',DISPLAY)
   INVERSION{i_s}.STABLE_BAYES.LIN.BIAS = struct('p',p,'o',o);
   invert_paresse_STABLE_BAYES(i_s,SESS,COND,'HYPERBOLIC','BIAS OFF',DISPLAY)
   INVERSION{i_s}.STABLE_BAYES.HYP.NO_BIAS = struct('p',p,'o',o);
   invert_paresse_STABLE_BAYES(i_s,SESS,COND,'HYPERBOLIC','BIAS ON',DISPLAY)   
   INVERSION{i_s}.STABLE_BAYES.HYP.BIAS = struct('p',p,'o',o);
   
   disp(['----- RANDOM MODELS -----'])

   % RANDOM MODEL
   [p,o] = invert_paresse_RAND(i_s,SESS,COND,'SHARED',DISPLAY)
   INVERSION{i_s}.RAND.SHARED = struct('p',p,'o',o);
   [p,o] = invert_paresse_RAND(i_s,SESS,COND,'COND',DISPLAY)
   INVERSION{i_s}.RAND.COND = struct('p',p,'o',o);
   [p,o] = invert_paresse_RAND(i_s,SESS,COND,'SESS',DISPLAY)
   INVERSION{i_s}.RAND.SESS = struct('p',p,'o',o);
   [p,o] = invert_paresse_RAND(i_s,SESS,COND,{'SESS','COND'})
   INVERSION{i_s}.RAND.SESS_COND = struct('p',p,'o',o);
   

   disp(['----- MATHYS MODELS -----'])

   % BAYESIAN + VOLATILITY
  [p,o] = invert_paresse_MATHYS(i_s,SESS,COND,'LINEAR','SHARED','LEARN VOLATILITY','BIAS OFF',DISPLAY) 
   INVERSION{i_s}.MATHYS.LIN.SHARED.LEARN_VOL.NO_BIAS = struct('p',p,'o',o);
   [p,o] = invert_paresse_MATHYS(i_s,SESS,COND,'LINEAR','SHARED','FIXED VOLATILITY','BIAS OFF',DISPLAY) 
   INVERSION{i_s}.MATHYS.LIN.SHARED.FIXED_VOL.NO_BIAS = struct('p',p,'o',o);
   [p,o] = invert_paresse_MATHYS(i_s,SESS,COND,'LINEAR','DIFFERENT','LEARN VOLATILITY','BIAS OFF',DISPLAY) 
   INVERSION{i_s}.MATHYS.LIN.DIFFERENT.LEARN_VOL.NO_BIAS = struct('p',p,'o',o);
   [p,o] = invert_paresse_MATHYS(i_s,SESS,COND,'LINEAR','DIFFERENT','FIXED VOLATILITY','BIAS OFF',DISPLAY) 
   INVERSION{i_s}.MATHYS.LIN.DIFFERENT.FIXED_VOL.NO_BIAS = struct('p',p,'o',o);
   [p,o] = invert_paresse_MATHYS(i_s,SESS,COND,'HYPERBOLIC','SHARED','LEARN VOLATILITY','BIAS OFF',DISPLAY) 
   INVERSION{i_s}.MATHYS.HYP.SHARED.LEARN_VOL.NO_BIAS= struct('p',p,'o',o);
   [p,o] = invert_paresse_MATHYS(i_s,SESS,COND,'HYPERBOLIC','SHARED','FIXED VOLATILITY','BIAS OFF',DISPLAY) 
   INVERSION{i_s}.MATHYS.HYP.SHARED.FIXED_VOL.NO_BIAS = struct('p',p,'o',o);   
  [p,o] = invert_paresse_MATHYS(i_s,SESS,COND,'HYPERBOLIC','DIFFERENT','LEARN VOLATILITY','BIAS OFF',DISPLAY) 
   INVERSION{i_s}.MATHYS.HYP.DIFFERENT.LEARN_VOL.NO_BIAS = struct('p',p,'o',o);
   [p,o] = invert_paresse_MATHYS(i_s,SESS,COND,'HYPERBOLIC','DIFFERENT','FIXED VOLATILITY','BIAS OFF',DISPLAY) 
   INVERSION{i_s}.MATHYS.HYP.DIFFERENT.FIXED_VOL.NO_BIAS = struct('p',p,'o',o);  
  [p,o] = invert_paresse_MATHYS(i_s,SESS,COND,'LINEAR','SHARED','LEARN VOLATILITY','BIAS ON',DISPLAY) 
   INVERSION{i_s}.MATHYS.LIN.SHARED.LEARN_VOL.BIAS = struct('p',p,'o',o);
   [p,o] = invert_paresse_MATHYS(i_s,SESS,COND,'LINEAR','SHARED','FIXED VOLATILITY','BIAS ON',DISPLAY) 
   INVERSION{i_s}.MATHYS.LIN.SHARED.FIXED_VOL.BIAS = struct('p',p,'o',o);
  [p,o] = invert_paresse_MATHYS(i_s,SESS,COND,'LINEAR','DIFFERENT','LEARN VOLATILITY','BIAS ON',DISPLAY)
  INVERSION{i_s}.MATHYS.LIN.DIFFERENT.LEARN_VOL.BIAS = struct('p',p,'o',o);
  [p,o] = invert_paresse_MATHYS(i_s,SESS,COND,'LINEAR','DIFFERENT','FIXED VOLATILITY','BIAS ON',DISPLAY)
  INVERSION{i_s}.MATHYS.LIN.DIFFERENT.FIXED_VOL.BIAS = struct('p',p,'o',o);
  [p,o] = invert_paresse_MATHYS(i_s,SESS,COND,'HYPERBOLIC','SHARED','LEARN VOLATILITY','BIAS ON',DISPLAY)
  INVERSION{i_s}.MATHYS.HYP.SHARED.LEARN_VOL.BIAS = struct('p',p,'o',o);
  [p,o] = invert_paresse_MATHYS(i_s,SESS,COND,'HYPERBOLIC','SHARED','FIXED VOLATILITY','BIAS ON',DISPLAY)
  INVERSION{i_s}.MATHYS.HYP.SHARED.FIXED_VOL.BIAS = struct('p',p,'o',o);
  [p,o] = invert_paresse_MATHYS(i_s,SESS,COND,'HYPERBOLIC','DIFFERENT','LEARN VOLATILITY','BIAS ON',DISPLAY) 
  INVERSION{i_s}.MATHYS.HYP.DIFFERENT.LEARN_VOL.BIAS = struct('p',p,'o',o);
  [p,o] = invert_paresse_MATHYS(i_s,SESS,COND,'HYPERBOLIC','DIFFERENT','FIXED VOLATILITY','BIAS ON',DISPLAY) 
  INVERSION{i_s}.MATHYS.HYP.DIFFERENT.FIXED_VOL.BIAS = struct('p',p,'o',o);  

   
   
end

save INVERSION_paresse INVERSION


