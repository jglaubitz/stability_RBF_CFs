%% script_PHS_Genz
%
% Description: 
% Script to perform a stability and error analysis of PHS for Genz' test functions 
%
% Author: Jan Glaubitz 
% Date: June 22, 2021 

clc, clear

%% Free parameters 
CC = 100; % number of tests for Genz 
noise_level = 0; % amount of uniform noise (0 means no noise, while a>0 mean 10^(-a))
points = 'equid'; % type of data points (equid, Halton, random) 
kernel = 'TPS'; % kernel (G, MQ, IQ, Wendland, TPS, cubic, quintic)  
ep = 1; % shape parameter 
order = 2; % order (for Wendland function)
d = 1; % polynomial degree 

%% prepare script 

% dimension and precision 
a = 0; b = 1; dim = 2; % domain is [0,1]^2
precision = 32; % use usual double precision 
NN = (2:10:62).^2;

% values of interest 
cond_nr = []; % condition number 
opt = []; % optimal values 
s = []; % stability values 
err = []; % errors 

for i=1:length(NN)
    
    %% Update data points 
    N = NN(i); % number of data points 
    [i, length(NN)]
    
    % data points, and RBb 
    X = generate_points( dim, a, b, N, points); % generate data points 
    rbf = initialize_RBF( kernel, dim, order ); % initialize RBF 
    
    %% Compute the condition number 
    cond_number = Cond( a, b, rbf, ep, X, d ); 

    %% Compute moments of the RBFs 
    m_RBF = RBF_moments( a, b, kernel, rbf, ep, X );
    %m_RBF2 = RBF_moments( a, b, 'numint', rbf, ep, X );
    
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

%% Illustrate results - Genz 1
figure(1) 
p = plot( NN,s,'b-.', NN,opt,'k:', NN,err(:,1),'r-'); 
set(p, 'LineWidth',3.5)
set(gca, 'FontSize', 26)  % Increasing ticks fontsize
xlim([ NN(1); NN(end) ]) 
xlabel('$N$','Interpreter','latex') 
xticks(10.^[1 2 3]);
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
p = plot( NN,s,'b-.', NN,opt,'k:', NN,err(:,2),'r-'); 
set(p, 'LineWidth',3.5)
set(gca, 'FontSize', 26)  % Increasing ticks fontsize
xlim([ NN(1); NN(end) ]) 
xlabel('$N$','Interpreter','latex') 
xticks(10.^[1 2 3 4]);
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
p = plot( NN,s,'b-.', NN,opt,'k:', NN,err(:,3),'r-'); 
set(p, 'LineWidth',3.5)
set(gca, 'FontSize', 26)  % Increasing ticks fontsize
xlim([ NN(1); NN(end) ]) 
xlabel('$N$','Interpreter','latex') 
xticks(10.^[1 2 3 4]);
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
p = plot( NN,s,'b-.', NN,opt,'k:', NN,err(:,4),'r-'); 
set(p, 'LineWidth',3.5)
set(gca, 'FontSize', 26)  % Increasing ticks fontsize
xlim([ NN(1); NN(end) ]) 
xlabel('$N$','Interpreter','latex') 
xticks(10.^[1 2 3 4]);
set(gca, 'XScale', 'log') 
set(gca, 'YScale', 'log') 
if d >= 0
    lgnd = legend('$\|C_N\|_{\infty}$','$\|I\|_{\infty}$','$| C_N[g_4] - I[g_4] |$','Interpreter','latex','Location','best');
else
    lgnd = legend('$\|C_N\|_{\infty}$','$C_N[1]$','$| C_N[g_4] - I[g_4] |$','Interpreter','latex','Location','best');
end
set(lgnd, 'Interpreter','latex', 'FontSize',26, 'color','none')
grid on

%% Illustrate results - Genz 1 & 4 combined 
figure(5) 
p = plot( NN,s,'b-.', NN,err(:,1), 'k--', NN,err(:,4),'r-'); 
set(p, 'LineWidth',3.5) 
% yline for optimal stability value
aux_line = yline(opt(1),':','$\|I\|_{\infty}$','Interpreter','latex'); 
set(aux_line, 'LineWidth',3.5, 'FontSize',26 );
aux_line.LabelVerticalAlignment = 'bottom'; 
lgnd = legend('$\|C_N\|_{\infty}$','$| C_N[g_1] - I[g_1] |$','$| C_N[g_4] - I[g_4] |$','Interpreter','latex','Location','southwest');
% axes 
set(gca, 'FontSize', 26)  % Increasing ticks fontsize
xlim([ NN(1); NN(end) ]) 
xlabel('$N$','Interpreter','latex') 
xticks(10.^[1 2 3 4]); 
set(gca, 'XScale', 'log') 
set(gca, 'YScale', 'log') 
if d >= 0
    lgnd = legend('$\|C_N\|_{\infty}$','$| C_N[g_1] - I[g_1] |$','$| C_N[g_4] - I[g_4] |$','Interpreter','latex','Location','southwest');
else
    lgnd = legend('$\|C_N\|_{\infty}$','$| C_N[g_1] - I[g_1] |$','$| C_N[g_4] - I[g_4] |$','Interpreter','latex','Location','southwest');
end
set(lgnd, 'Interpreter','latex', 'FontSize',26, 'color','none')
grid on