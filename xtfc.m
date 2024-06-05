function [y, ydot, training_err] = xtfc(weight, bias, type_act, x, t_tot, y0, dydt_problem, J_problem, IterMax, IterTol)
% type_act: activation function
%{
1= Logistic;
2= TanH;
3= Sine;
4= Cosine;
5= Gaussian; the best so far w/ m=11
6= ArcTan;
7= Hyperbolic Sine;
8= SoftPlus
9= Bent Identity;
10= Inverse Hyperbolic Sine
11= Softsign
%}

L = length(weight);
n_x = length(x);
n_y = length(y0);
n_t = length(t_tot);

y = zeros(n_t, n_y);
ydot = zeros(n_t, n_y);

y(1, :) = y0;
ydot(1, :) = zeros(1, n_y);

[h, hd] = act(x, weight, bias, type_act);
h_0 = h - h(1, :);

training_err = 0;
for i = 1:(n_t - 1)
    ci = (x(end) - x(1)) / (t_tot(i + 1) - t_tot(i));
    %[t_tot(i), y0]

    betai = zeros(L, n_y);

    % Build Constrained Expressions
    yi = h_0*betai + y0;
    ydoti = ci*hd*betai;

    ydoti_ = dydt_problem(yi);
    % Build the Losses
    Loss = ydoti - ydoti_;

    % X-TFC ILS loop
    l2 = [2 1];
    iter = 0;
    while abs(l2(2)) > IterTol && abs(l2(1) - l2(2)) > IterTol &&  iter < IterMax
        l2(1) = l2(2);
        J = J_problem(h_0, hd, ci, yi);

        % xi variation
        dbetai = lsqminnorm(J, reshape(Loss, [n_x*n_y, 1]));

        % update xi
        betai = betai - reshape(dbetai, [L, n_y]);

        % Re-Build Constrained Expressions
        yi = h_0*betai + y0;
        ydoti = ci*hd*betai;

        ydoti_ = dydt_problem(yi);
        % Re-Build the Losses
        Loss = ydoti - ydoti_;

        l2(2) = norm(Loss);
        iter = iter + 1;
    end

    y0 = yi(end, :);

    y(i + 1, :) = y0;
    ydot(i + 1, :) = ydoti(end, :);
    training_err = training_err + sqrt(mean(Loss(:, 1).^2)) + sqrt(mean(Loss(:, 2).^2)) + sqrt(mean(Loss(:, 3).^2));
end

training_err = training_err/(length(t_tot) - 1);

end

