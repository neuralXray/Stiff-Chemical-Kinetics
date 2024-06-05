clear
format longE

rng('default') % set random seed


% =======================================
% Input parameters

t0 = 0;  % initial time
tf = 40;  % 40, 120, 5000: final time
t_step = 0.01;
t_step_ode15 = 0.0001;
t_tot = (t0:t_step:tf)';
y0 = [0.066, 0, 0, 0, 0.066, 0.002, 0];
L = 20;  % number of neurons
LB = - 1;  % Lower boundary for weight and bias samplings
UB = 1; %  Upper boundary for weight and bias samplings
% Activation function definition
weight = unifrnd(LB, UB, L, 1);
bias = unifrnd(LB, UB, L, 1);
type_act = 2;  % activation functions
x = linspace(0, 1, 20)';
% iterative least-square parameters
IterMax = 100;
IterTol = 1e-9;


% =======================================
% x-tfc
tStart = tic;  % start stopwatch timer
[y, ydot, training_err_xtfc_paper] = xtfc(weight, bias, type_act, x, t_tot, y0, @belousov_zhabotinsky_function, @belousov_zhabotinsky_j_function, IterMax, IterTol);
xtfc_elapsedtime = toc(tStart);
fprintf('The elapsed time for x-tfc is: %g \n', xtfc_elapsedtime);

ygrad = gradient(y')'./gradient(t_tot);

ydot_ = belousov_zhabotinsky_function(y);
Loss_dot = ydot - ydot_;
Loss_grad = ygrad - ydot_;

training_err_xtfc_dot = sum(sqrt(mean(Loss_dot.^2)));
training_err_xtfc_grad = sum(sqrt(mean(Loss_grad.^2)));


%% =======================================
% MATLAB ode15s solver

options = odeset('RelTol', IterTol);
% test the automatic detection of a DAE.

tStart = tic;
[t_15s, y_15s] = ode15s(@belousov_zhabotinsky_ode15s_function, (t0:t_step_ode15:tf)', y0', options);
ode15s_elapsedtime = toc(tStart);

fprintf('The elapsed time for ode15s is: %g \n', ode15s_elapsedtime );

ydot = gradient(y_15s')'./gradient(t_15s);

ydot_ = belousov_zhabotinsky_function(y_15s);
Loss = ydot - ydot_;

training_err_ode15s = sum(sqrt(mean(Loss.^2)));


%% =======================================
% MATLAB ode15i solver

yp0 = [0; 0; 0; 0; 0; 0; 0];

tStart = tic;
[y0, yp0] = decic(@belousov_zhabotinsky_ode15i_function, 0, y0, [0 1 0 0 0 0 0], yp0, [], options); 
[t_15i, y_15i] = ode15i(@belousov_zhabotinsky_ode15i_function, (t0:t_step_ode15:tf)', y0', yp0, options);
ode15i_elapsedtime = toc(tStart);

fprintf('The elapsed time for ode15i is: %g \n', ode15i_elapsedtime);

ydot = gradient(y_15i')'./gradient(t_15i);

ydot_ = belousov_zhabotinsky_function(y_15s);
Loss = ydot - ydot_;

training_err_ode15i = sum(sqrt(mean(Loss.^2)));



% =======================================
% Errors

fprintf('\n')
fprintf('The average training error for X-TFC (paper method) is: %g \n', training_err_xtfc_paper)
fprintf('The average training error for X-TFC is: %g \n', training_err_xtfc_dot)
fprintf('The average training error for X-TFC (gradient(y)) is: %g \n', training_err_xtfc_grad)
fprintf('The average training error for ode15s is: %g \n', training_err_ode15s)
fprintf('The average training error for ode15i is: %g \n', training_err_ode15i)


% =======================================
% PLOTS

subplot(3, 3, 1)
set(gca, 'Fontsize', 12)
hold on
grid on
plot(t_tot, y(:, 1), 'LineWidth', 3 , 'Color', [17 50 50]/100)
plot(t_15i, y_15i(:,1), '--', 'LineWidth', 3, 'Color', [90 61.5 25]/100)
plot(t_15s, y_15s(:,1), ':', 'LineWidth', 3, 'Color', [80 32.5 9]/100)
xlim([t0 tf])
ylabel('concentration')
box on
title('y1', 'FontWeight', 'Normal')

subplot(3, 3, 2)
set(gca, 'Fontsize', 12)
hold on
grid on
plot(t_tot, y(:, 2), 'LineWidth', 3, 'Color', [17 50 50]/100)
plot(t_15i, y_15i(:,2), '--', 'LineWidth', 3, 'Color', [90 61.5 25]/100)
plot(t_15s, y_15s(:,2), ':', 'LineWidth', 3, 'Color', [80 32.5 9]/100)
xlim([t0 tf])
box on
title('y2', 'FontWeight', 'Normal')

subplot(3, 3, 3)
set(gca, 'Fontsize', 12)
hold on
grid on
plot(t_tot, y(:, 3), 'LineWidth', 3, 'Color', [17 50 50]/100)
plot(t_15i, y_15i(:,3), '--', 'LineWidth', 3, 'Color', [90 61.5 25]/100)
plot(t_15s, y_15s(:,3), ':', 'LineWidth', 3, 'Color', [80 32.5 9]/100)
xlim([t0 tf])
legend('X-TFC', 'ode15i', 'ode15s')
box on
title('y3', 'FontWeight', 'Normal')

subplot(3, 3, 4)
set(gca, 'Fontsize', 12)
hold on
grid on
plot(t_tot, y(:, 4), 'LineWidth', 3 , 'Color', [17 50 50]/100)
plot(t_15i, y_15i(:,4), '--', 'LineWidth', 3, 'Color', [90 61.5 25]/100)
plot(t_15s, y_15s(:,4), ':', 'LineWidth', 3, 'Color', [80 32.5 9]/100)
xlim([t0 tf])
ylabel('concentration')
box on
title('y4', 'FontWeight', 'Normal')

subplot(3, 3, 5)
set(gca, 'Fontsize', 12)
hold on
grid on
plot(t_tot, y(:, 5), 'LineWidth', 3, 'Color', [17 50 50]/100)
plot(t_15i, y_15i(:, 5), '--', 'LineWidth', 3, 'Color', [90 61.5 25]/100)
plot(t_15s, y_15s(:, 5), ':', 'LineWidth', 3, 'Color', [80 32.5 9]/100)
xlim([t0 tf])
xlabel('time (s)')
box on
title('y5', 'FontWeight', 'Normal')

subplot(3, 3, 6)
set(gca, 'Fontsize', 12)
hold on
grid on
plot(t_tot, y(:, 6), 'LineWidth', 3, 'Color', [17 50 50]/100)
plot(t_15i, y_15i(:, 6), '--', 'LineWidth', 3, 'Color', [90 61.5 25]/100)
plot(t_15s, y_15s(:, 6), ':', 'LineWidth', 3, 'Color', [80 32.5 9]/100)
xlim([t0 tf])
xlabel('time (s)')
box on
title('y6', 'FontWeight', 'Normal')

subplot(3 ,3, 7)
set(gca, 'Fontsize', 12)
hold on
grid on
plot(t_tot, y(:, 7), 'LineWidth', 3 , 'Color', [17 50 50]/100)
plot(t_15i, y_15i(:, 7), '--', 'LineWidth', 3, 'Color', [90 61.5 25]/100)
plot(t_15s, y_15s(:, 7), ':', 'LineWidth', 3, 'Color', [80 32.5 9]/100)
xlim([t0 tf])
xlabel('time (s)')
ylabel('concentration')
box on
title('y7', 'FontWeight', 'Normal')


figure(2)
set(gca, 'Fontsize', 12)
hold on
grid on
title('X-TFC solutions')
plot(t_tot, y(:, 2), 'LineWidth', 2, 'Color', [17 50 50]/100)
plot(t_tot, y(:, 4), 'LineWidth', 2, 'Color', [80 32.5 9]/100)
plot(t_tot, y(:, 6), 'LineWidth', 2, 'Color', [90 61.5 25]/100)
plot(t_tot, y(:, 7), 'LineWidth', 2)
xlabel('time(s)')
waitfor(legend('Y', 'P', 'Z', 'Q'))

