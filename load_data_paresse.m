function [ Y,U,IsYout,cond,sess,Counts] = load_data_paresse(i_subject,SESS,COND)
%loads data for RL inversion


% Inversion of data from the PARESSE tasks

% Different conditions COND
% A : high effort, left better
% B : low effort , right better
% C : right harder, high reward
% D : left harder, low reward

% Session
% 96 trials : 24 per condition

% Modeling
% - Value of effort : 80% or 20%
% - Reward : 0.4 or 0.2
% - Probabilities of high reward/effort : .75 or .25

%--- Formating data

% Experimenter data per session : u
% - 1 : choice
% - 2 : reward
% - 3 : effort
% Separate data for the bayesian model under stability assumption : counts
% - 1 : total # of right choices
% - 2: total # of left choices
% - 3: # of big R (wins) | right choices
% - 4: # of big R (wins) | left choices
% - 5: # of small E (wins) | right choices
% - 6: # of small E (wins) | left choices

% --- u(4:9) is used for Baysian Ideal Learner

% Data : y
% - 1 : choice

%--- Task contingencies
%fileparts(mfilename('fullpath'))
 %root = '/home/vasilisa/Desktop/PARESSE_Vincent/data/paresse_alldata/';
%root = '/mnt/data/IRMDATA/U610_PARESSE/behaviour/COMP_MODEL/PARESSE_Vincent/data/paresse_alldata/';
%root = 'C:\Users\vincent.adam\Desktop\Vincent ADAM\Matlab\Demo_toolbox\Data\PARESSE\data\paresse_alldata\';

if isunix || ismac
root = [fileparts(mfilename('fullpath')),'/data/'];
elseif ispc
root = [fileparts(mfilename('fullpath')),'\data\'];
end

N =24;
%loading one session


Y =[];
U = [];
Counts = [];
cond = []; % vector of conditions
sess =[];
% data columns
% - 1 : index of subject
% - 2 : index of trial
% - 3 : index of condition
% - 4 : --- (time)
% - 5 : ---
% - 6 : ---
% - 7 : ---
% - 8 : ---
% - 9 : ---
% - 10 : ---
% - 11 : choice : left or right (binary)
% - 12 : ---
% - 13 : ? (format : 0.5222)
% - 14 : Reward (0.2 -0.4 )
% - 15 : effort (0.2 -0.8 )
% - 16 : ? (format : 8.9797e3)
% - 17 : ? (format : 1.4895e4)

for i_sess = SESS
    
    filename = ['BMCNtestSub',num2str(i_subject),'ses',num2str(i_sess)]; %A
    load([root,filename])
    
    for i_cond = COND
        
        I = find(data(:,3)==i_cond);
        u = zeros(3,24); % u = zeros(9,24); 
        
        % binaries rewards and efforts:  rewards 0.4 = 1 and 0.2 = 0 and efforts 0.8 = 1 and 0.2 = 0
        R=data(I,14);
        I_highR = find(R==0.4);
        R(:) =0; R(I_highR)=1;
        
        E=data(I,15);
        I_highE = find(E==0.8);
        E(:)=0;E(I_highE)=1;
                
        y = (data(I,11)-1)'; % choice to binary
        
        u(1,:) = y; % choice
        u(2,:)=R;
        u(3,:)=E;
        
        
        counts = zeros(6,24); % count how many big rewards and small efforts (wins) where per each hand
        counts(1,:)=cumsum(y); % total right choices
        counts(2,:)=cumsum(~y); % total left choices
        counts(3,:)=cumsum(y&R');
        counts(4,:)=cumsum(~y&R');
        counts(5,:)=cumsum(y&(~E')); % small efforts given right choice
        counts(6,:)=cumsum(~y&(~E'));
  %      u(4:9,:)=counts;
        
        
        Y = [Y;y];
        U = [U;u];
        Counts = [Counts;counts];
        sess = [sess; i_sess];

    end
    
    cond = [cond;COND'];
    
end

IsYout = 0*Y;


end

