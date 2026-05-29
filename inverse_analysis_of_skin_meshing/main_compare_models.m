clc; clear; close all;

%% Run both models
% res1 = run_GOH_1fiber;
res2 = run_GOH_2fiber;
res3 = run_Ogden;

fprintf('\nModel comparison (SSE)\n')
fprintf('-----------------------\n')
% fprintf('1F GOH   SSE = %.6e\n', res1.SSE)
fprintf('2F GOH   SSE = %.6e\n', res2.SSE)
fprintf('Ogden    SSE = %.6e\n', res3.SSE)

fprintf('\nStandard deviation of residuals\n')
fprintf('--------------------------------\n')
% fprintf('1F GOH   std = %.6e\n', res1.stdError)
fprintf('2F GOH   std = %.6e\n', res2.stdError)
fprintf('Ogden    std = %.6e\n\n', res3.stdError)

%% Load experimental data
% Get root folder from function location
stretch_exp_perp = res2.stretch_exp_perp;
stress_exp_perp  = res2.stress_exp_perp;

stretch_exp_para = res2.stretch_exp_para;
stress_exp_para  = res2.stress_exp_para;


%% Plot settings
fontSize   = 20;
lineWidth1 = 3;
markerSize = 8;

figure; hold on; box on; grid on;

%% --- Experimental data ---
h_exp1 = plot(stretch_exp_perp, stress_exp_perp, ...
    'ro', 'MarkerSize', markerSize, 'LineWidth', 1.5);

h_exp2 = plot(stretch_exp_para, stress_exp_para, ...
    'go', 'MarkerSize', markerSize, 'LineWidth', 1.5);
h_mean = plot(res3.stretch_exp_mean, res3.stress_exp_mean,'k--','LineWidth',2);
%% --- 1-Fiber (solid lines) ---
% h1 = plot(res1.stretch_perp, res1.stress_perp, ...
%     'r-', 'LineWidth', lineWidth1);
% 
% h2 = plot(res1.stretch_para, res1.stress_para, ...
%     'g-', 'LineWidth', lineWidth1);

%% --- 2-Fiber (triangle markers) ---
h3 = plot(res2.stretch_perp, res2.stress_perp, ...
    'r-^', 'LineWidth', lineWidth1, ...
    'MarkerSize', 6, 'MarkerIndices', 1:5:length(res2.stretch_perp));

h4 = plot(res2.stretch_para, res2.stress_para, ...
    'g-^', 'LineWidth', lineWidth1, ...
    'MarkerSize', 6, 'MarkerIndices', 1:5:length(res2.stretch_para));

%% --- Ogden (blue dashed squares) ---
h5 = plot(res3.stretch_perp, res3.stress_perp, ...
    'b--s', 'LineWidth', lineWidth1, ...
    'MarkerSize', 6, ...
    'MarkerIndices', 1:5:length(res3.stretch_perp));

%% Formatting
xlabel('\lambda Stretch [.]','FontSize',fontSize)
ylabel('\sigma Cauchy stress [MPa]','FontSize',fontSize)
title('Comparison of GOH vs Ogden Model','FontSize',fontSize)

% legend([h_exp1 h_exp2 h_mean h1 h2 h3 h4 h5], ...
%     {'Exp perp','Exp para','Exp mean', ...
%  '1F perp','1F para', ...
%  '2F perp','2F para', ...
%  'Ogden'})
legend([h_exp1 h_exp2 h_mean h3 h4 h5], ...
{'Exp perp','Exp para','Exp mean', ...
 'GOH (2F) perp','GOH (2F) para', ...
 'Ogden'})

set(gca,'FontSize',fontSize)
axis square
axis tight
Add main model comparison script
