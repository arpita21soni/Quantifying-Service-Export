clc;
clear;
close all;

disp('Running factual...');
run('factual.m');

% Line 930 Counter = productivity of H is same as in 1995
% u_h_t = repmat(u_h_t(:,1), 1, size(u_h_t,2));
% disp('Running counterfactual...');
% run('counter1.m');

% Line 1010 Counter = investment efficiency is same as in 1995
% Ax_t = repmat(Ax_t(:,1), 1, size(Ax_t,2));
% disp('Running counterfactual...');
% run('counter2.m');

% Line 1010 Counter = trade cost is same as in 1995
% tao_h_t = repmat(tao_h_t(:,:,1), 1, size(tao_h_t,3));
% disp('Running counterfactual...');
% run('counter2.m');

disp('All done.');
