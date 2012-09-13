function [ fx ] = f_paresse_RL(x,P,u,in)
% INPUT
% - x : Qvalues
%     1 : QRleft
%     2 : QRright
%     3 : QEleft
%     4 : QEright
% - P : parameters
%     1 : learning rate for reward (sigmoid mapped)
%     2 : learning rate for effort (sigmoid mapped)
% - u : input
%     1 : choice (left = 0, right = 1)
%     2 : reward
%     3 : effort

fx = x; % copy all
R = u(2); % reward obtained
E = u(3); % effort done

if isequal(in.TYPE,'QL')
alphaR = sigm(P(1),struct('INV',0)); % learning rate for reward
alphaE = sigm(P(2),struct('INV',0)); % learning rate for effort
elseif isequal(in.TYPE,'WSLS')
alphaR = 1; % stick to last feedback
alphaE = 1; % stick to last feedback
end
 
% update
a = u(1); % the chosen side (0 = left, 1 = right)
iu =[1,3]+a; % [1,3] or [2,4] which side to update
fx(iu(1)) = x(iu(1)) + alphaR*(R-x(iu(1))); % update reward
fx(iu(2)) = x(iu(2)) + alphaE*(E-x(iu(2))); % update effort

end

