
% Script to perform inversion for mixed modesl - models that have different learning rules for rewards and effort outcomes several models on the data-set of the
% task 'paresse' (mbb, 2011-2012, Vincent, Vasilisa). 


%%
%---------------------------------------------------------
% Performing the inversion for all subjects and all models
%---------------------------------------------------------

SESS = [1,2,3]; % 3 sessions per subject
COND = [1,2,3,4]; % A,B,C,D
SUBJECTS =[1:20]; % 20 subjects

Nsubjects = length(SUBJECTS);

%

for i_s = SUBJECTS
    
    disp('--------------------------------------------------------')
    disp(['------Subject ',num2str(i_s),'/',num2str(length(SUBJECTS))])
    
    INVERSION{i_s}.i_subject = i_s;
 
   
    % ------ QL FOR REWARDS and WSLS for EFFORT LINEAR   
    [posterior,out] = invert_paresse_mixed(i_s,SESS,COND,'QL','WSLS','LINEAR');
    INVERSION{i_s}.QLrew_WSeff.LIN.posterior = posterior;
    INVERSION{i_s}.QLrew_WSeff.LIN.out = out;
    
    % ------  QL FOR REWARDS and WSLS for EFFORT HYPERBOLIC   
    [posterior,out] = invert_paresse_mixed(i_s,SESS,COND,'QL','WSLS','HYPERBOLIC');
    INVERSION{i_s}.QLrew_WSLSeff.HYP.posterior = posterior;
    INVERSION{i_s}.QLrew_WSLSeff.HYP.out = out;
   
%     
    % ----QL FOR EFFORTS and WSLS for REWARD LINEAR    
    [posterior,out] = invert_paresse_mixed(i_s,SESS,COND,'WSLS','QL','LINEAR');
    INVERSION{i_s}.QLeff_WSrew.LIN.posterior = posterior;
    INVERSION{i_s}.QLeff_WSrew.LIN.out = out;
%    
    % ----QL FOR EFFORTS and WSLS for REWARD HYPERBOLIC   
    [posterior,out] = invert_paresse_mixed(i_s,SESS,COND,'WSLS','QL','HYPERBOLIC');
    INVERSION{i_s}.QLeff_WSrew.HYP.posterior = posterior;
    INVERSION{i_s}.QLeff_WSrew.HYP.out = out;
    
    % ---- QL FOR REWARDS and IDEAL OBSERVER FOR EFFORTS LINEAR    
    [posterior,out] = invert_paresse_mixed(i_s,SESS,COND,'QL','IDEAL','LINEAR');
    INVERSION{i_s}.QLrew_IDEALeff.LIN.posterior = posterior;
    INVERSION{i_s}.QLrew_IDEALeff.LIN.out = out;
    
    % ---- QL FOR REWARDS and IDEAL OBSERVER FOR EFFORTS HYPERBOLIC   
    [posterior,out] = invert_paresse_mixed(i_s,SESS,COND,'QL','IDEAL','HYPERBOLIC');
    INVERSION{i_s}.QLrew_IDEALeff.HYP.posterior = posterior;
    INVERSION{i_s}.QLrew_IDEALeff.HYP.out = out;  
 
     
  % ---- QL FOR EFFORTS and IDEAL OBSERVER FOR REWARDS LINEAR     
    [posterior,out] = invert_paresse_mixed(i_s,SESS,COND,'IDEAL','QL','LINEAR');
    INVERSION{i_s}.QLeff_IDEALrew.LIN.posterior = posterior;
    INVERSION{i_s}.QLeff_IDEALrew.LIN.out = out;

  % ---- QL FOR EFFORTS and IDEAL OBSERVER FOR REWARDS HYPERBOLIC      
    [posterior,out] = invert_paresse_mixed(i_s,SESS,COND,'IDEAL','QL','HYPERBOLIC');
    INVERSION{i_s}.QLeff_IDEALrew.HYP.posterior = posterior;
    INVERSION{i_s}.QLeff_IDEALrew.HYP.out = out;
    
    
    % ---- WSLS FOR REWARDS and IDEAL OBSERVER FOR EFFORTS LINEAR       
    [posterior,out] = invert_paresse_mixed(i_s,SESS,COND,'WSLS','IDEAL','LINEAR');
    INVERSION{i_s}.WSrew_IDEALeff.LIN.posterior = posterior;
    INVERSION{i_s}.WSrew_IDEALeff.LIN.out = out;
    
    
    % ---- WSLS FOR REWARDS and IDEAL OBSERVER FOR EFFORTS HYPERBOLIC       
    [posterior,out] = invert_paresse_mixed(i_s,SESS,COND,'WSLS','IDEAL','HYPERBOLIC');
    INVERSION{i_s}.WSrew_IDEALeff.HYP.posterior = posterior;
    INVERSION{i_s}.WSrew_IDEALeff.LIN.out = out;
    
    
    % ---- WSLS FOR EFFORTS and IDEAL OBSERVER FOR REWARDS LINEAR        
    [posterior,out] = invert_paresse_mixed(i_s,SESS,COND,'IDEAL','WSLS','LINEAR');
    INVERSION{i_s}.WSeff_IDEALrew.LIN.posterior = posterior;
    INVERSION{i_s}.WSeff_IDEALrew.LIN.out = out;
    
    
    % ---- WSLS FOR EFFORTS and IDEAL OBSERVER FOR REWARDS HYPERBOLIC        
    [posterior,out] = invert_paresse_mixed(i_s,SESS,COND,'IDEAL','WSLS','HYPERBOLIC');
    INVERSION{i_s}.WSeff_IDEALrew.HYP.posterior = posterior;
    INVERSION{i_s}.WSeff_IDEALrew.HYP.out = out;
    
    % ---- WSLS FOR EFFORTS and WSLS  FOR REWARDS LINEAR        
    [posterior,out] = invert_paresse_mixed(i_s,SESS,COND,'WSLS','WSLS','LINEAR');
    INVERSION{i_s}.WSeff_WSrew.LIN.posterior = posterior;
    INVERSION{i_s}.WSeff_WSrew.LIN.out = out;
    
     % ---- WSLS FOR EFFORTS and WSLS  FOR REWARDS HYPERBOLIC        
    [posterior,out] = invert_paresse_mixed(i_s,SESS,COND,'WSLS','WSLS','HYPERBOLIC');
    INVERSION{i_s}.WSeff_WSrew.HYP.posterior = posterior;
    INVERSION{i_s}.WSeff_WSrew.HYP.out = out;
    
     % ---- WSLS FOR EFFORTS and WSLS  FOR REWARDS LINEAR 
      % ----- for the full QL model test the hypothesis about different or
    % shared Learning rates 
    
    [posterior,out] = invert_paresse_mixed(i_s,SESS,COND,'QL','QL','LINEAR','SHARED');
    INVERSION{i_s}.QLeff_QLrew.LIN.SHARED.posterior = posterior;
    INVERSION{i_s}.QLeff_QLrew.LIN.SHARED.out = out;
    
    [posterior,out] = invert_paresse_mixed(i_s,SESS,COND,'QL','QL','LINEAR','DIFFERENT');
    INVERSION{i_s}.QLeff_QLrew.LIN.DIFF.posterior = posterior;
    INVERSION{i_s}.QLeff_QLrew.LIN.DIFF.out = out;
    
    
     % ---- WSLS FOR EFFORTS and WSLS  FOR REWARDS HYPERBOLIC        
    [posterior,out] = invert_paresse_mixed(i_s,SESS,COND,'QL','QL','HYPERBOLIC','SHARED');
    INVERSION{i_s}.QLeff_QLrew.HYP.SHARED.posterior = posterior;
    INVERSION{i_s}.QLeff_QLrew.HYP.SHARED.out = out;
    
    [posterior,out] = invert_paresse_mixed(i_s,SESS,COND,'QL','QL','HYPERBOLIC','DIFFERENT');
    INVERSION{i_s}.QLeff_QLrew.HYP.DIFF.posterior = posterior;
    INVERSION{i_s}.QLeff_QLrew.HYP.DIFF.out = out;
    
     % ---- WSLS FOR EFFORTS and WSLS  FOR REWARDS LINEAR        
    [posterior,out] = invert_paresse_mixed(i_s,SESS,COND,'IDEAL','IDEAL','LINEAR');
    INVERSION{i_s}.IDEALeff_IDEALrew.LIN.posterior = posterior;
    INVERSION{i_s}.IDEALeff_IDEALrew.LIN.out = out;
    
     % ---- WSLS FOR EFFORTS and WSLS  FOR REWARDS HYPERBOLIC        
    [posterior,out] = invert_paresse_mixed(i_s,SESS,COND,'IDEAL','IDEAL','HYPERBOLIC');
    INVERSION{i_s}.IDEALeff_IDEALrew.HYP.posterior = posterior;
    INVERSION{i_s}.IDEALeff_IDEALrew.HYP.out = out;
    
    
%     % ------- BIAS PER HAND
%     
%     [posterior,out] = invert_paresse_mixed(i_s,SESS,COND,'BIAS_HAND');
%     INVERSION{i_s}.BIAS.HAND.posterior = posterior;
%     INVERSION{i_s}.BIAS.HAND.out = out;
% %     
%     
%     % ------- SIMPLE BIAS PER CONDITION 
%     
%     [posterior,out] = invert_paresse_mixed(i_s,SESS,COND,'BIAS_COND');
%     INVERSION{i_s}.BIAS.COND.posterior = posterior;
%     INVERSION{i_s}.BIAS.COND.out = out;
%     
%    
end
%

 save PARESSE_INVERSION_mixed INVERSION


% %% TO FIND the best parameters of the group - concatenate all the subjects and specify the model for which you want to do a group inversion 
% SUBJECTS = [1:20]; 
% SESS=[1,2,3];
% COND=[1,2,3,4]; 
% INVERSION_fix = struct(); 
% [posterior,out] = invert_fixed_effect_last(SUBJECTS,SESS,COND,'QL','LINEAR','DIFFERENT','BIAS OFF'); 
% INVERSION_fix.posterior = posterior;
% INVERSION_fix.out = out;
% 
% save PARESSE_INVERSION_fix_QL_LIN_DIFF INVERSION_fix 
% %[1,2] = subjects



%%
% reorganize inversions (from previously erroneous launch)
%
%
% for i_s = SUBJECTS
%     INVERSION{i_s}.SINGLE.REWARD =  INVERSION{i_s}.SHARED.REWARD;
%     INVERSION{i_s}.SINGLE.EFFORT =  INVERSION{i_s}.SHARED.EFFORT;
%     INVERSION{i_s}.SHARED = rmfield(INVERSION{i_s}.SHARED,{'REWARD','EFFORT'});
%     INVERSION{i_s}.DIFFERENT=   rmfield(INVERSION{i_s}.DIFFERENT,{'REWARD','EFFORT'});
%     INVERSION{i_s}= rmfield(INVERSION{i_s},'HYP');
% end

%%
%---------------------------------------------------------
% Model comparison at the group level
%load ....
% SESS = [1,2,3]; % 3 sessions per subject
% COND = [1,2,3,4]; % A,B,C,D
SUBJECTS =[1:20]; % 20 subjects

Nsubjects =length(SUBJECTS);

model_names = {'QLrew_WSeff_LIN',...
    'QLrew_WSeff_HYP',...
    'QLeff_WSrew_LIN',...
    'QLeff_WSrew_HYP',...
    'WSrew_IDEALeff_LIN',...
    'WSrew_IDEALeff_HYP',...
    'WSeff_IDEALrew_LIN',...
    'WSeff_IDEALrew_HYP',...
   'QLrew_IDEALeff_LIN',...
   'QLrew_IDEALeff_HYP', ...
   'QLeff_QLrew_LIN_SHARED',...
   'QLeff_QLrew_LIN_DIFF'
   'QLeff_QLrew_HYP_SHARED',
   'QLeff_QLrew_HYP_DIFF',
   'WSeff_WSrew_LIN',
   'WSeff_WSrew_HYP',
   'IDEALeff_IDEALrew_LIN',
   'IDEALeff_IDEALrew_HYP'};


Nmodels = size(model_names,2);
%LogEv = zeros(Nsubjects,Nmodels); % Matrix of Log evidences for each subject and model
LogEv = zeros(Nsubjects,Nmodels);

for i_s = 1:Nsubjects
    %-------------- 
   % LOG EVIDENCE pour les modeles 
    LogEv(i_s,1) = INVERSION{i_s}.QLrew_WSeff_LIN.out.F;
    LogEv(i_s,2) = INVERSION{i_s}.QLrew_WSeff_HYP.out.F;
    
    
    LogEv(i_s,3) = INVERSION{i_s}.QLeff_WSrew_LIN.out.F;
    LogEv(i_s,4) = INVERSION{i_s}.QLeff_WSrew_HYP.out.F;
    
    LogEv(i_s,5) = INVERSION_sing{i_s}.WSrew_IDEALeff_LIN.out.F;
    LogEv(i_s,6) = INVERSION_sing{i_s}.WSrew_IDEALeff_HYP.out.F;
    
    LogEv(i_s,7) = INVERSION{i_s}.WSeff_IDEALrew_LIN.out.F;
    LogEv(i_s,8) = INVERSION{i_s}.WSeff_IDEALrew_HYP.out.F;
    
    LogEv(i_s,9) = INVERSION{i_s}.QLrew_IDEALeff_LIN.out.F;
    LogEv(i_s,10) =INVERSION{i_s}.QLrew_IDEALeff_HYP.out.F;
    
    LogEv(i_s,11) = INVERSION{i_s}.QLeff_IDEALrew_LIN.out.F;
    LogEv(i_s,12) =INVERSION{i_s}.QLeff_IDEALrew_HYP.out.F;
    
    LogEv(i_s,13) = INVERSION{i_s}.QLrew_QLeff_LIN_SHARED.out.F;
    LogEv(i_s,14) =INVERSION{i_s}.QLrew_QLeff_HYP_SHARED.out.F;
   
    LogEv(i_s,15) = INVERSION{i_s}.QLrew_QLeff_LIN_DIFF.out.F;
    LogEv(i_s,16) =INVERSION{i_s}.QLrew_QLeff_HYP_DIFF.out.F;
   
    LogEv(i_s,17) = INVERSION{i_s}.WSrew_WSeff_LIN.out.F;
    LogEv(i_s,18) =INVERSION{i_s}.WSrew_WSeff_HYP.out.F;
    
     LogEv(i_s,19) =INVERSION{i_s}.IDEALrew_IDEALeff_LIN.out.F;
    LogEv(i_s,20) =INVERSION{i_s}.IDEALrew_IDEALeff_HYP.out.F;
end

%%
%-- Comparing all models individually
out_allmodels=group_level_analysis(LogEv,'RFX',model_names);
pause(1); % to update display

%-- Comparing by groups

%--- HYPERBOLIC vs LINEAIRE 
partition_LIN_vs_HYP = {[1,3,7,9,11,13,15,17,19],[2,4,8,10,12,14,16,18,20]};
out_LIN_vs_HYP=group_level_analysis(LogEv,'RFX',{'LINEAR','HYPERBOLIC'},partition_LIN_vs_HYP);
pause(1); % to update display

%--- QL reward or QL effort vs QL both outcomes 
partition_QLfull_vs_QLpartial = {[1,2,9,10],[3,4,11,12],[13,14,15,16]};
out_QLfull_vs_QLpartial=group_level_analysis(LogEv,'RFX',{'QL_reward','QL_effort','QL_rew&eff'},partition_QLfull_vs_QLpartial)
pause(1)

partition_SHARED_vs_DIFF = {[13,14],[15,16]};
out_SHARED_vs_DIFF=group_level_analysis(LogEv,'RFX',{'QL_reward','QL_effort','QL_rew&eff'},partition_SHARED_vs_DIFF)
pause(1)

%% 
% muPhi=zeros(numel(SUBJECTS),2);
% muPhi2=zeros(numel(SUBJECTS),2);
% alphasLIN=zeros(numel(SUBJECTS),2);
% alphasHYP=zeros(numel(SUBJECTS),2);
% 
% for i_s=SUBJECTS
%     muPhi(i_s,1)=INVERSION{1,i_s}.QL.DIFFERENT.LIN.posterior.muPhi(1,1); 
%     muPhi(i_s,2)=INVERSION{1,i_s}.QL.DIFFERENT.LIN.posterior.muPhi(2,1);
%     
%     
%     
%     muPhi2(i_s,1)=INVERSION{1,i_s}.QL.DIFFERENT.HYP.posterior.muPhi(1,1); 
%     muPhi2(i_s,2)=INVERSION{1,i_s}.QL.DIFFERENT.HYP.posterior.muPhi(2,1);
%     % convert to the range from 0 to 1 
%     
%     alphasLIN(:,1) = sigm(muPhi(:,1));
%     alphasLIN(:,2) = sigm(muPhi(:,2));
%     
%     alphasHYP(:,1) = sigm(muPhi2(:,1));
%     alphasHYP(:,2) = sigm(muPhi2(:,2));
%     
% end
%     
%     
    
%% EXTRACT THE CORRELATION between parameters 

% correlation_theta = zeros(1,numel(SUBJECTS)); 
% 
% for i_s=1:length(SUBJECTS)
%     
%     posterior=INVERSION{1,i_s}.QL.DIFFERENT.LIN.posterior; 
%     out=INVERSION{1,i_s}.QL.DIFFERENT.LIN.out;    
%     [Ctheta,Cphi]=extract_paramscorr(posterior, out); 
%    
%     correlation_theta(1,i_s)=Ctheta(1,2); % correlation between alphas 
%     
%     VBA_ReDisplay(posterior,out)
% end
    
    
    
    
    
    
