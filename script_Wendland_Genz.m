%% script_S716_Genz
%
% Description: 
% Script to perform a stability and error analysis of Wendland's compactly supported
% RBFs for Genz' test functions 
%
% Author: Jan Glaubitz 
% Date: June 22, 2021 

clc, clear

%% Free parameters 
N = 20^2; % number of data points 
CC = 100; % number of tests for Genz 
noise_level = 0; % amount of uniform noise (0 means no noise, while a>0 mean 10^(-a))
points = 'equid'; % type of data points (equid, Halton, random) 
kernel = 'Wendland'; % kernel (G, MQ, IQ, Wendland, TPS, cubic) 
order = 1; % order (for Wendland function)
d = 0; % polynomial degree 

%% prepare script 

% dimension and precision 
a = 0; b = 1; dim = 2; % domain is [0,1]^2
precision = 32; % use usual double precision 

% domain, data points, and RBb 
X = generate_points( dim, a, b, N, points); % generate data points 
DM = DistanceMatrix( X, X ); % matrix containing the distances between points
h = min( DM + 42*eye(N), [], 'all' ); % measure of distance between data points
rbf = initialize_RBF( kernel, dim, order ); % initialize RBF 

% set up shape parameters 
mid = ceil((1/h)^(1/10)); % exponent corresponding to h
exp_start = mid-6; exp_end = mid+1; n = 100;
Ep = logspace(exp_start,exp_end,n)'; % generates a vector of n logarithmically spaced points

% values of interest 
cond_nr = []; % condition number 
opt = []; % optimal values 
s = []; % stability values 
err = []; % errors 

for i=1:n
    
    %% Update shape parameter 
    ep = Ep(i); % vector of shape parameters  
    [i, n]
    
    %% Compute the condition number 
    cond_number = Cond( a, b, rbf, ep, X, d ); 

    %% Compute moments of the RBFs 
    m_RBF = RBF_moments( a, b, kernel, rbf, ep, X );
    
    %% Compute RBF-CF weights (in a different manner)
    w = compute_weights( a, b, rbf, ep, X, m_RBF, d, precision );
    
    %% Compute the stability measure s_N and otimal value C_N[1] 
    if d < 0 
        opt_value = sum(w); % optimal value 
    else 
        opt_value = (b-a)^2; % optimal value
    end
    stab_measure = sum(abs(w)); % stability measure 
    
    %% Genz test functions 
    error = Genz( a, b, X, w, CC, noise_level ); % average errors 
    
    %% Store values 
    cond_nr = [cond_nr; cond_number]; % condition number 
    opt = [opt; opt_value]; % optimal values 
    s = [s; stab_measure]; % stability values 
    err = [err; error]; % error of CF 
    
end 

%format shortE
%[min_error,index] = min( err(:,1) );
%[p_count, order, d, 1, min_error, Ep(index), s(index)]
%[min_error,index] = min( err(:,4) );
%[p_count, order, d, 4, min_error, Ep(index), s(index)]

%% Illustrate results - Genz 1
figure(1) 
p = plot( Ep,s,'b-.', Ep,opt,'k:', Ep,err(:,1),'r-'); 
line = xline(1/h,'-','1/h', 'LineWidth',3.5); 
set(line, 'LineWidth',3.5, 'FontSize',26 )
set(p, 'LineWidth',3.5)
set(gca, 'FontSize', 26)  % Increasing ticks fontsize
xlim([ Ep(1); Ep(end) ]) 
xlabel('$\varepsilon$','Interpreter','latex') 
xticks(10.^[exp_start+1 0.5*(exp_start+1+exp_end) exp_end]);
set(gca, 'XScale', 'log') 
set(gca, 'YScale', 'log') 
if d >= 0
    lgnd = legend('$\|C_N\|_{\infty}$','$\|I\|_{\infty}$','$| C_N[g_1] - I[g_1] |$','Interpreter','latex','Location','best');
else
    lgnd = legend('$\|C_N\|_{\infty}$','$C_N[1]$','$| C_N[g_1] - I[g_1] |$','Interpreter','latex','Location','best');
end
set(lgnd, 'Interpreter','latex', 'FontSize',26, 'color','none')
grid on

%% Illustrate results - Genz 2
figure(2) 
p = plot( Ep,s,'b-.', Ep,opt,'k:', Ep,err(:,2),'r-'); 
line = xline(1/h,'-','1/h', 'LineWidth',3.5); 
set(line, 'LineWidth',3.5, 'FontSize',26 )
set(p, 'LineWidth',3.5)
set(gca, 'FontSize', 26)  % Increasing ticks fontsize
xlim([ Ep(1); Ep(end) ]) 
xlabel('$\varepsilon$','Interpreter','latex') 
xticks(10.^[exp_start+1 0.5*(exp_start+exp_end+1) exp_end]);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log') 
if d >= 0
    lgnd = legend('$\|C_N\|_{\infty}$','$\|I\|_{\infty}$','$| C_N[g_2] - I[g_2] |$','Interpreter','latex','Location','best');
else
    lgnd = legend('$\|C_N\|_{\infty}$','$C_N[1]$','$| C_N[g_2] - I[g_2] |$','Interpreter','latex','Location','best');
end
set(lgnd, 'Interpreter','latex', 'FontSize',26, 'color','none')
grid on

%% Illustrate results - Genz 3
figure(3) 
p = plot( Ep,s,'b-.', Ep,opt,'k:', Ep,err(:,3),'r-'); 
line = xline(1/h,'-','1/h', 'LineWidth',3.5); 
set(line, 'LineWidth',3.5, 'FontSize',26 )
set(p, 'LineWidth',3.5)
set(gca, 'FontSize', 26)  % Increasing ticks fontsize
xlim([ Ep(1); Ep(end) ]) 
xlabel('$\varepsilon$','Interpreter','latex') 
xticks(10.^[exp_start+1 0.5*(exp_start+exp_end+1) exp_end]);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log') 
if d >= 0
    lgnd = legend('$\|C_N\|_{\infty}$','$\|I\|_{\infty}$','$| C_N[g_3] - I[g_3] |$','Interpreter','latex','Location','best');
else
    lgnd = legend('$\|C_N\|_{\infty}$','$C_N[1]$','$| C_N[g_3] - I[g_3] |$','Interpreter','latex','Location','best');
end
set(lgnd, 'Interpreter','latex', 'FontSize',26, 'color','none')
grid on

%% Illustrate results - Genz 4
figure(4) 
p = plot( Ep,s,'b-.', Ep,opt,'k:', Ep,err(:,4),'r-'); 
line = xline(1/h,'-','1/h', 'LineWidth',3.5); 
set(line, 'LineWidth',3.5, 'FontSize',26 )
set(p, 'LineWidth',3.5)
set(gca, 'FontSize', 26)  % Increasing ticks fontsize
xlim([ Ep(1); Ep(end) ]) 
xlabel('$\varepsilon$','Interpreter','latex') 
xticks(10.^[exp_start+1 0.5*(exp_start+exp_end+1) exp_end]);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log') 
if d >= 0
    lgnd = legend('$\|C_N\|_{\infty}$','$\|I\|_{\infty}$','$| C_N[g_4] - I[g_4] |$','Interpreter','latex','Location','best');
else
    lgnd = legend('$\|C_N\|_{\infty}$','$C_N[1]$','$| C_N[g_4] - I[g_4] |$','Interpreter','latex','Location','best');
end
set(lgnd, 'Interpreter','latex', 'FontSize',26, 'color','none')
grid on