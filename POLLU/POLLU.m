clear
format longE

rng('default') % set random seed


% =======================================
% Input parameters
t0 = 0; % initial time
tf = 60; % final time
t_tot = (t0:0.001:tf)';
y0 = [0, 0.2, 0, 0.04, 0, 0, 0.1, 0.3, 0.01, 0, 0, 0, 0, 0, 0, 0, 0.007, 0, 0, 0];  % Initial Values
L = 10;  % number of neurons
LB = -1; % Lower boundary for weight and bias samplings
UB = 1; % Upper boundary for weight and bias samplings
weight = unifrnd(LB, UB, L, 1);
bias = unifrnd(LB, UB, L, 1);
type_act = 2; % type activation function
x = linspace(0, 1, 10)';
% Iterative least-square parameters
IterMax = 100;
IterTol = 1e-9;


% =======================================
% x-tfc
tStart = tic;  % start stopwatch timer
[y, ydot, training_err_xtfc_paper] = xtfc(weight, bias, type_act, x, t_tot, y0, @pollu_function, @pollu_j_function, IterMax, IterTol);
xtfc_elapsedtime = toc(tStart);
fprintf('The elapsed time for x-tfc is: %g \n', xtfc_elapsedtime);

ygrad = gradient(y')'./gradient(t_tot);

ydot_ = pollu_function(y);
Loss_dot = ydot - ydot_;
Loss_grad = ygrad - ydot_;

training_err_xtfc_dot = sum(sqrt(mean(Loss_dot.^2)));
training_err_xtfc_grad = sum(sqrt(mean(Loss_grad.^2)));


%% =======================================
% MATLAB ode15s solver

options = odeset('RelTol', IterTol);
% test the automatic detection of a DAE.

tStart = tic;
[t_15s, y_15s] = ode15s(@pollu_ode15s_function, t_tot, y0', options);
ode15s_elapsedtime = toc(tStart);

fprintf('The elapsed time for ode15s is: %g \n', ode15s_elapsedtime );

ydot = gradient(y_15s')'./gradient(t_tot);

ydot_ = pollu_function(y_15s);
Loss = ydot - ydot_;

training_err_ode15s = sum(sqrt(mean(Loss.^2)));


% =======================================
% Errors
err_ode15s = abs(y_15s - y);

fprintf('\n')
fprintf('The average training error for X-TFC (paper method) is: %g \n', training_err_xtfc_paper)
fprintf('The average training error for X-TFC is: %g \n', training_err_xtfc_dot)
fprintf('The average training error for X-TFC (gradient(y)) is: %g \n', training_err_xtfc_grad)
fprintf('The average training error for ode15s is: %g \n', training_err_ode15s)


% =======================================
% PLOTS
figure(1)
subplot(4, 5, 1)
hold on
grid on
plot(t_15s, y_15s(:, 1), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 1), '--', 'LineWidth', 3, 'Color', [17 50 50]/100)
ylabel('y1')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 2)
hold on
grid on
plot(t_15s, y_15s(:, 2), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 2), '--', 'LineWidth', 3, 'Color', [17 50 50]/100)
ylabel('y2')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 3)
hold on
grid on
plot(t_15s, y_15s(:, 3), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 3), '--', 'LineWidth', 3, 'Color', [17 50 50]/100)
ylabel('y3')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 4)
hold on
grid on
plot(t_15s, y_15s(:, 4), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 4), '--', 'LineWidth', 3, 'Color', [17 50 50]/100)
ylabel('y4')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 5)
hold on
grid on
plot(t_15s, y_15s(:, 5), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 5), '--', 'LineWidth', 3, 'Color', [17 50 50]/100)
ylabel('y5')
hYLabel = get(gca, 'YLabel');
legend('X-TFC','ode15s')
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 6)
hold on
grid on
plot(t_15s, y_15s(:, 6), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 6), '--', 'LineWidth', 3, 'Color', [17 50 50]/100)
ylabel('y6')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 7)
hold on
grid on
plot(t_15s, y_15s(:, 7), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 6), '--', 'LineWidth', 3, 'Color', [17 50 50]/100)
ylabel('y7')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 8)
hold on
grid on
plot(t_15s, y_15s(:, 8), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 8), '--', 'LineWidth', 3, 'Color', [17 50 50]/100)
ylabel('y8')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 9)
hold on
grid on
plot(t_15s, y_15s(:, 9), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 9), '--', 'LineWidth', 3, 'Color', [17 50 50]/100)
ylabel('y9')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 10)
hold on
grid on
plot(t_15s, y_15s(:, 10), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 10), '--', 'LineWidth', 3, 'Color', [17 50 50]/100)
ylabel('y10')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 11)
hold on
grid on
plot(t_15s, y_15s(:, 11), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 11), '--', 'LineWidth', 3, 'Color', [17 50 50]/100)
ylabel('y11')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 12)
hold on
grid on
plot(t_15s, y_15s(:, 12), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 12), '--', 'LineWidth', 3, 'Color', [17 50 50]/100)
ylabel('y12')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 13)
hold on
grid on
plot(t_15s, y_15s(:, 13), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 13), '--', 'LineWidth', 3, 'Color', [17 50 50]/100)
ylabel('y13')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 14)
hold on
grid on
plot(t_15s, y_15s(:, 14), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 14), '--', 'LineWidth', 3, 'Color', [17 50 50]/100)
ylabel('y14')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 15)
hold on
grid on
plot(t_15s, y_15s(:, 15), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 15), '--', 'LineWidth', 3, 'Color', [17 50 50]/100)
ylabel('y15')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 16)
hold on
grid on
plot(t_15s, y_15s(:, 16), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 16), '--', 'LineWidth', 3, 'Color', [17 50 50]/100)
ylabel('y16')
xlabel('time (min)')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 17)
hold on
grid on
plot(t_15s, y_15s(:, 17), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 17), '--', 'LineWidth', 3, 'Color', [17 50 50]/100)
ylabel('y17')
xlabel('time (min)')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 18)
hold on
grid on
plot(t_15s, y_15s(:, 18), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 18), '--', 'LineWidth', 3, 'Color', [17 50 50]/100)
ylabel('y18')
xlabel('time (min)')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 19)
hold on
grid on
plot(t_15s, y_15s(:, 19), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 19), '--', 'LineWidth', 3, 'Color', [17 50 50]/100)
xlabel('time (min)')
ylabel('y19')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 20)
hold on
grid on
plot(t_15s, y_15s(:, 20), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 20), '--', 'LineWidth', 3, 'Color', [17 50 50]/100)
xlabel('time (min)')
ylabel('y20')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on


figure(2)

subplot(4, 5, 1)
hold on
grid on
plot(t_15s, err_ode15s(:, 1), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
xlim([t0 tf])
ylabel('y1')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
set(gca, 'YScale', 'log')
box on

subplot(4, 5, 2)
hold on
grid on
plot(t_15s, err_ode15s(:, 2), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
xlim([t0 tf])
ylabel('y2')
set(gca, 'YScale', 'log')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 3)
hold on
grid on
plot(t_15s, err_ode15s(:, 3), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
xlim([t0 tf])
ylabel('y3')
set(gca, 'YScale', 'log')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 4)
hold on
grid on
plot(t_15s, err_ode15s(:, 4), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
xlim([t0 tf])
ylabel('y4')
set(gca, 'YScale', 'log')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 5)
hold on
grid on
plot(t_15s, err_ode15s(:, 5), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
xlim([t0 tf])
ylabel('y5')
set(gca, 'YScale', 'log')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 6)
hold on
grid on
plot(t_15s, err_ode15s(:, 6), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
xlim([t0 tf])
ylabel('y6')
set(gca, 'YScale', 'log')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 7)
hold on
grid on
plot(t_15s, err_ode15s(:, 7), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
xlim([t0 tf])
ylabel('y7')
set(gca, 'YScale', 'log')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 8)
hold on
grid on
plot(t_15s, err_ode15s(:, 8), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
xlim([t0 tf])
ylabel('y8')
set(gca, 'YScale', 'log')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 9)
hold on
grid on
plot(t_15s, err_ode15s(:, 9), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
xlim([t0 tf])
ylabel('y9')
set(gca, 'YScale', 'log')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 10)
hold on
grid on
plot(t_15s, err_ode15s(:, 10), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
xlim([t0 tf])
ylabel('y10')
set(gca, 'YScale', 'log')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 11)
hold on
grid on
plot(t_15s, err_ode15s(:, 11), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
xlim([t0 tf])
ylabel('y11')
set(gca, 'YScale', 'log')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 12)
hold on
grid on
plot(t_15s, err_ode15s(:, 12), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
xlim([t0 tf])
ylabel('y12')
set(gca, 'YScale', 'log')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on


subplot(4, 5, 13)
hold on
grid on
plot(t_15s, err_ode15s(:, 13), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
xlim([t0 tf])
ylabel('y13')
set(gca, 'YScale', 'log')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 14)
hold on
grid on
plot(t_15s, err_ode15s(:, 14), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
xlim([t0 tf])
ylabel('y14')
set(gca, 'YScale', 'log')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 15)
hold on
grid on
plot(t_15s, err_ode15s(:, 15), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
xlim([t0 tf])
ylabel('y15')
set(gca, 'YScale', 'log')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 16)
hold on
grid on
plot(t_15s, err_ode15s(:, 16), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
xlim([t0 tf])
xlabel('time (min)')
ylabel('y16')
set(gca, 'YScale', 'log')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 17)
hold on
grid on
plot(t_15s, err_ode15s(:, 17), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
xlim([t0 tf])
xlabel('time (min)')
ylabel('y17')
set(gca, 'YScale', 'log')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 18)
hold on
grid on
plot(t_15s, err_ode15s(:, 18), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
xlim([t0 tf])
xlabel('time (min)')
ylabel('y18')
set(gca, 'YScale', 'log')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 19)
hold on
grid on
plot(t_15s, err_ode15s(:, 19), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
xlim([t0 tf])
xlabel('time (min)')
ylabel('y19')
set(gca, 'YScale', 'log')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on

subplot(4, 5, 20)
hold on
grid on
plot(t_15s, err_ode15s(:, 20), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
xlim([t0 tf])
xlabel('time (min)')
set(gca, 'YScale', 'log')
hYLabel = get(gca, 'YLabel');
set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
box on
waitfor(ylabel('y20'))

