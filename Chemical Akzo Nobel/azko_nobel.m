clear
format longE

rng('default')  % set random seed


% =======================================
% Input parameters
t0 = 0;  % initial time
tf = 180;  % final time
t_tot = (t0:0.01:tf)';
% Initial Values
y0 = [0.437, 0.00123, 0, 0, 0, 0.367];
L = 20;  % number of neurons
LB = -1;  % Lower boundary for weight and bias samplings
UB = 1;  % Upper boundary for weight and bias samplings
% Activation function definition
weight = unifrnd(LB, UB, L, 1);
bias = unifrnd(LB, UB, L, 1);
type_act = 2;  % type activation function
x = linspace(0, 1, 50)';
% iterative least-square parameters
IterMax = 100;
IterTol = 1e-12;


% =======================================
% x-tfc
tStart = tic;  % start stopwatch timer
[y, ydot, training_err_xtfc_paper] = xtfc(weight, bias, type_act, x, t_tot, y0, @akzonobel_function, @akzonobel_j_function, IterMax, IterTol);
xtfc_elapsedtime = toc(tStart);
fprintf('The elapsed time for x-tfc is: %g \n', xtfc_elapsedtime);

ygrad = gradient(y')'./gradient(t_tot);

ydot_ = akzonobel_function(y);
Loss_dot = ydot - ydot_;
Loss_grad = ygrad - ydot_;

training_err_xtfc_dot = sum(sqrt(mean(Loss_dot.^2)));
training_err_xtfc_grad = sum(sqrt(mean(Loss_grad.^2)));


%% =======================================
% MATLAB ode15s solver

options = odeset('RelTol',IterTol);
% test the automatic detection of a DAE.

tStart = tic;
[t_15s, y_15s] = ode15s(@akzonobel_ode15s_function, t_tot', y0', options);
ode15s_elapsedtime = toc(tStart);

fprintf('The elapsed time for ode15s is: %g \n', ode15s_elapsedtime );

ydot = gradient(y_15s')'./gradient(t_tot);

ydot_ = akzonobel_function(y_15s);
Loss = ydot - ydot_;

training_err_ode15s = sum(sqrt(mean(Loss.^2)));


% =======================================
%% errors
err_ode15s = abs(y_15s - y);

fprintf('\n')
fprintf('The average training error for X-TFC (paper method) is: %g \n', training_err_xtfc_paper)
fprintf('The average training error for X-TFC is: %g \n', training_err_xtfc_dot)
fprintf('The average training error for X-TFC (gradient(y)) is: %g \n', training_err_xtfc_grad)
fprintf('The average training error for ode15s is: %g \n', training_err_ode15s)


% =======================================
%% PLOTS
subplot(4, 3, 1)
set(gca,'Fontsize', 12)
hold on
grid on
plot(t_15s, y_15s(:,1), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 1), '--','LineWidth', 3 , 'Color', [17 50 50]/100)
xlim([t0 tf])
ylabel('concentration')
box on
title('y1', 'FontWeight', 'Normal')

subplot(4, 3, 2)
set(gca, 'Fontsize', 12)
hold on
grid on
plot(t_15s, y_15s(:,2), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 2), '--', 'LineWidth', 3, 'Color', [17 50 50]/100)
xlim([t0 tf])
box on
title('y2', 'FontWeight', 'Normal')

subplot(4, 3, 3)
set(gca, 'Fontsize', 12)
hold on
grid on
plot(t_15s, y_15s(:,3), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 3),'--','LineWidth', 3, 'Color', [17 50 50]/100)
xlim([t0 tf])
legend('ode15s', 'X-TFC')
box on
title('y3', 'FontWeight', 'Normal')


subplot(4, 3, 4)
set(gca, 'Fontsize', 12)
hold on
grid on
plot(t_15s, err_ode15s(:, 1), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
xlim([t0 tf])
ylabel('abs(error)')
set(gca, 'YScale', 'log')
box on

subplot(4, 3, 5)
set(gca,'Fontsize', 12)
hold on
grid on
plot(t_15s, err_ode15s(:, 2), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
xlim([t0 tf])
set(gca, 'YScale', 'log')
box on

subplot(4, 3, 6)
set(gca, 'Fontsize', 12)
hold on
grid on
plot(t_15s,err_ode15s(:, 3), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
xlim([t0 tf])
set(gca, 'YScale', 'log')
box on


subplot(4, 3, 7)
set(gca,'Fontsize', 12)
hold on
grid on
plot(t_15s,y_15s(:, 4), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot,y(:, 4), '--', 'LineWidth', 3 , 'Color', [17 50 50]/100)
xlim([t0 tf])
ylabel('concentration')
box on
title('y4', 'FontWeight', 'Normal')

subplot(4, 3, 8)
set(gca, 'Fontsize', 12)
hold on
grid on
plot(t_15s, y_15s(:,5), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 5), '--', 'LineWidth', 3, 'Color', [17 50 50]/100)
xlim([t0 tf])
box on
title('y5', 'FontWeight', 'Normal')

subplot(4, 3, 9)
set(gca, 'Fontsize', 12)
hold on
grid on
plot(t_15s, y_15s(:,6), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 6),'--', 'LineWidth', 3, 'Color', [17 50 50]/100)
xlim([t0 tf])
box on
title('y6', 'FontWeight', 'Normal')


subplot(4, 3, 10)
set(gca, 'Fontsize', 12)
hold on
grid on
plot(t_15s, err_ode15s(:, 4), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
xlim([t0 tf])
xlabel('time (s)')
ylabel('abs(error)')
set(gca, 'YScale', 'log')
box on

subplot(4, 3, 11)
set(gca, 'Fontsize', 12)
hold on
grid on
plot(t_15s, err_ode15s(:, 5), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
xlim([t0 tf])
xlabel('time (s)')
set(gca, 'YScale', 'log')
box on

subplot(4, 3, 12)
set(gca, 'Fontsize', 12)
hold on
grid on
plot(t_15s, err_ode15s(:, 6), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
xlim([t0 tf])
xlabel('time (s)')
set(gca, 'YScale', 'log')
box on


figure(2)
set(gca, 'Fontsize', 12)
hold on
grid on
plot(t_tot, y(:, 1), t_tot, y(:, 2), t_tot, y(:, 3), t_tot, y(:, 4), t_tot, y(:, 5), t_tot, y(:, 6), 'LineWidth', 1.5)
xlim([t0 tf])
xlabel('time(s)')
legend('y1', 'y2', 'y3', 'y4', 'y5', 'y6')
waitfor(title('All concentrations', 'FontWeight', 'Normal'))

