function [act, actd] = act(x, w, b, type_Activation)

switch type_Activation
    case 1 % logistic
        act = 1./(exp(- b' - x*w') + 1);
        actd = (w'.*exp(- b' - x*w'))./(exp(- b' - x*w') + 1).^2;

    case 2 % tanH
        act = (exp(b' + x*w') - exp(- b' - x*w'))./(exp(b' + x*w') + exp(- b' - x*w'));
        actd = (w'.*exp(b' + x*w') + w'.*exp(- b' - x*w'))./(exp(b' + x*w') + exp(- b' - x*w')) - ((exp(b' + x*w') - exp(- b' - x*w')).*(w'.*exp(b' + x*w') - w'.*exp(- b' - x*w')))./(exp(b' + x*w') + exp(- b' - x*w')).^2;

    case 3 % sin
        act = sin(b' + x*w');
        actd = w'.*cos(b' + x*w');

    case 4 % cos
        act = cos(b' + x*w');
        actd = -w'.*sin(b' + x*w');

    case 5 % Gaussian
        act = exp(-(b' + x*w').^2);
        actd = -2*w'.*exp(-(b' + x*w').^2).*(b' + x*w');

    case 6 % ArcTan
        act = atan(b' + x*w');
        actd = w'./((b' + x*w').^2 + 1);

    case 7 % Hyp Sine
        act = sinh(b' + x*w');
        actd = w'.*cosh(b' + x*w');

    case 8 % SoftPlus
        act = log(exp(b' + x*w') + 1);
        actd = (w'.*exp(b' + x*w'))./(exp(b' + x*w') + 1);

    case 9 % Bent Identity
        act = b' + ((b' + x*w').^2 + 1).^(1/2)/2 + x*w' - 1/2;
        actd = w + (w*(b' + x*w'))/(2*((b' + x*w').^2 + 1).^(1/2));

    case 10 % Inverse Hyp Sine
        act = log(b' + ((b' + x*w').^2 + 1).^(1/2) + x*w');
        actd = (w' + (w*(b' + x*w'))./((b' + x*w').^2 + 1).^(1/2))./(b' + ((b' + x*w').^2 + 1).^(1/2) + x*w');

    case 11 % Softsign
        act = (b' + x*w')/(abs(b' + x*w') + 1);
        actd = w'./(abs(b' + x*w') + 1) - (w'.*sign(b' + x*w').*(b' + x*w'))./(abs(b' + x*w') + 1).^2;

end

end

