clear
format longE

rng('default') % set random seed


% =======================================
% Input parameters
t_0 = 1e-5;  % initial time
t_f = 1e5;  % final time
t_step = 2000;
t_tot = (t_0:t_step:t_f)';
n_t = length(t_tot);
t_tot = zeros(n_t + 1, 1);
t_tot(2:end) = logspace(log10(t_0), log10(t_f), n_t)';
y0 = [1, 0, 0];  % Initial Values
L = 10;  % number of neurons
LB = - 1;  % Lower boundary for weight and bias samplings
UB = 1;  % Upper boundary for weight and bias samplings
% Activation function definition
weight = unifrnd(LB, UB, L, 1);
bias = unifrnd(LB, UB, L, 1);
type_act = 2;  % type activation function
x = linspace(0, 1, L)';
% Iterative least-square parameters
IterMax = 100;
IterTol = 1e-12;


% =======================================
% x-tfc
tStart = tic;  % start stopwatch timer
[y, ydot, training_err_xtfc_paper] = xtfc(weight, bias, type_act, x, t_tot, y0, @rober_function, @rober_j_function, IterMax, IterTol);
xtfc_elapsedtime = toc(tStart);
fprintf('The elapsed time for x-tfc is: %g \n', xtfc_elapsedtime);

ygrad = gradient(y')'./gradient(t_tot);

ydot_ = rober_function(y);
Loss_dot = ydot - ydot_;
Loss_grad = ygrad - ydot_;

training_err_xtfc_dot = sum(sqrt(mean(Loss_dot.^2)));
training_err_xtfc_grad = sum(sqrt(mean(Loss_grad.^2)));


% =======================================
% MATLAB ode15s solver

options = odeset('RelTol', IterTol);
% test the automatic detection of a DAE.

tStart = tic;
[t_15s, y_15s] = ode15s(@rober_ode15s_function, t_tot, y0', options);
ode15s_elapsedtime = toc(tStart);

fprintf('The elapsed time for ode15s is: %g \n', ode15s_elapsedtime);

ydot = gradient(y_15s')'./gradient(t_tot);

ydot_ = rober_function(y_15s);
Loss = ydot - ydot_;

training_err_ode15s = sum(sqrt(mean(Loss.^2)));


% =======================================
% MATLAB ode15i solver

y0_ode15i = [1; 0; 1e-3];
yp0 = [0; 0; 0];

tStart = tic;
[y0, yp0] = decic(@rober_ode15i_function, 0, y0_ode15i, [1 1 0], yp0, [], options); 
[t_15i, y_15i] = ode15i(@rober_ode15i_function, t_tot, y0_ode15i, yp0, options);
ode15i_elapsedtime = toc(tStart);

fprintf('The elapsed time for ode15i is: %g \n', ode15i_elapsedtime);

ydot = gradient(y_15i')'./gradient(t_tot);

ydot_ = rober_function(y_15s);
Loss = ydot - ydot_;

training_err_ode15i = sum(sqrt(mean(Loss.^2)));


% =======================================
% Errors
err_ode15s = abs(y_15s - y);

fprintf('\n')
fprintf('The average training error for X-TFC (paper method) is: %g \n', training_err_xtfc_paper)
fprintf('The average training error for X-TFC is: %g \n', training_err_xtfc_dot)
fprintf('The average training error for X-TFC (gradient(y)) is: %g \n', training_err_xtfc_grad)
fprintf('The average training error for ode15s is: %g \n', training_err_ode15s)
fprintf('The average training error for ode15i is: %g \n', training_err_ode15i)


% =======================================
% PLOTS
subplot(3, 2, 1)
set(gca, 'Fontsize', 12)
hold on  % add one more curve or subplot to this same figure
grid on
plot(t_15s, y_15s(:, 1), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 1), '--', 'LineWidth', 3 , 'Color', [17 50 50]/100)
ylabel('y_1')
set(gca, 'XScale', 'log')
box on
title('concentration', 'FontWeight', 'Normal')
legend('ode15s', 'X-TFC')
xticklabels([])

subplot(3, 2, 2)
set(gca, 'Fontsize', 12)
hold on
grid on
plot(t_15s, err_ode15s(:, 1), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
ylabel('abs(error)')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
box on
title('abs(error)', 'FontWeight', 'Normal')
xticklabels([])

subplot(3, 2, 3)
set(gca, 'Fontsize', 12)
hold on
grid on
plot(t_15s, y_15s(:,2), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 2), '--', 'LineWidth', 3, 'Color', [17 50 50]/100)
set(gca, 'XScale', 'log')
ylabel('y_2')
box on
xticklabels([])

subplot(3, 2, 4)
set(gca, 'Fontsize', 12)
hold on
grid on
plot(t_15s, err_ode15s(:, 2), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
box on
xticklabels([])


subplot(3, 2, 5)
set(gca, 'Fontsize', 12)
hold on
grid on
plot(t_15s, y_15s(:,3), 'LineWidth', 3, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 3), '--', 'LineWidth', 3, 'Color', [17 50 50]/100)
xlabel('time (s)')
set(gca, 'XScale', 'log')
ylabel('y_3')
box on

subplot(3, 2, 6)
set(gca, 'Fontsize', 12)
hold on
grid on
plot(t_15s, err_ode15s(:, 3), '*', 'LineWidth', 1.3, 'Color', [17 50 50]/100)
xlabel('time (s)')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
box on



figure(2)
set(gca, 'Fontsize', 12)
hold on
grid on
plot(t_tot, y(:, 1), '-o')
xlabel('time(s)')
set(gca, 'XScale','log')
xlim([1e-4, t_tot(end)])
title('y1', 'FontWeight', 'Normal')

figure(3)
set(gca, 'Fontsize', 12)
hold on
grid on
plot(t_tot, y(:, 2), '-o')
xlabel('time(s)')
set(gca, 'XScale', 'log')
xlim([1e-4, t_tot(end)])
title('y2', 'FontWeight', 'Normal')

figure(4)
set(gca, 'Fontsize', 12)
hold on
grid on
plot(t_tot, y(:, 3), '-o')
xlabel('time(s)')
set(gca, 'XScale', 'log')
xlim([1e-4, t_tot(end)])
waitfor(title('y3', 'FontWeight', 'Normal'))

