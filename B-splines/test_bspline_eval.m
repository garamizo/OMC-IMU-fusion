%% Test range
clear

bord = 6;  % basis order

tt = linspace(0, 5, 1000)';
dt = 0.11;
kk = (tt - tt(1)) / dt;

% minimum number of knots needed
nknot = ceil((tt(end) - tt(1)) / dt + 1e-10) + bord - 1;

figure, set(gcf, 'position', [339   476   913   391])

for i = 1 : nknot

    c = zeros(nknot, 1);
    c(i) = 1;

    y = bspline_eval(c, kk, bord);

    plot(tt, y)
    hold on
end
ylim([0, 1])
grid on
title(sprintf('%d Order', bord))

tlast = (nknot(end) - bord + 1) * dt;  % last possible interpolation
kk = (tlast - tt(1)) / dt;
ypass = bspline_eval(c, kk - 1e-10, bord);
% yfail = bspline_eval(c, kk + 1e-10, bord);

xline(tlast, '--')

%% Test composition
clear

bord = 6;  % basis order

tt = linspace(0, 5.1, 1000)';
dt = 0.1;
kk = (tt - tt(1)) / dt;

% minimum number of knots needed
nknot = ceil((tt(end) - tt(1)) / dt + 1e-10) + bord - 1;
c = randn(nknot, 1);

y = bspline_eval(c, kk, bord);

figure, set(gcf, 'position', [339   476   913   391])
plot(tt, y, 'k', 'linewidth', 2)
ylim([-2, 2])
grid on
title(sprintf('%d Order', bord))

hold on
for i = 1 : nknot
    ceq = zeros(nknot, 1);
    ceq(i) = c(i);
    y = bspline_eval(ceq, kk, bord);
    
    plot(tt, y)
end

%% Test multiple dimensions
clear

bord = 6;  % basis order

tt = linspace(0, 5.1, 1000)';
dt = 0.1;
kk = (tt - tt(1)) / dt;

% minimum number of knots needed
nknot = ceil((tt(end) - tt(1)) / dt + 1e-10) + bord - 1;
c = randn(nknot, 3);

y = bspline_eval(c, kk, bord);

figure, set(gcf, 'position', [339   476   913   391])
plot(tt, y)
ylim([-2, 2])
grid on
title(sprintf('%d Order', bord))

%% Test derivative
clear

bord = 6;  % basis order

tt = linspace(0, 5.1, 1000)';
Fs = 1/(tt(2)-tt(1));
dt = 0.1;
kk = (tt - tt(1)) / dt;

% minimum number of knots needed
nknot = ceil((tt(end) - tt(1)) / dt + 1e-10) + bord - 1;
c = polyval([1, -1, -1, 1], (1 : nknot)'/nknot);
c = [c, 1-2*c];  % add 2nd channel

[y, yd] = bspline_eval(c, kk, bord);
yd = yd / dt;
% yd2 = deriv_sgolay(detrend(y,0), Fs, [2, 7]);
yd2 = [y(2,:)-y(1,:); (y(3:end,:)-y(1:end-2,:))/2; y(end,:)-y(end-1,:)] * Fs;

figure, set(gcf, 'position', [339   476   913   391])
h1 = subplot(211);
plot(tt, y), grid on

h2 = subplot(212);
plot(tt, yd)
hold on, plot(tt, yd2, '--')
legend('analytical', 'numerical (ref)')
grid on
title(sprintf('%d Order', bord))

%% Test 2nd derivative
clear

bord = 6;  % basis order

tt = linspace(0, 5.1, 1000)';
Fs = 1/(tt(2)-tt(1));
dt = 0.1;
kk = (tt - tt(1)) / dt;

% minimum number of knots needed
nknot = ceil((tt(end) - tt(1)) / dt + 1e-10) + bord - 1;
c = polyval([1, -1, -1, 1], (1 : nknot)'/nknot);
c = [c, 1-2*c];  % add 2nd channel

[y, yd, ydd] = bspline_eval(c, kk, bord);
yd = yd / dt;
ydd = ydd / dt^2;
[yd2, ~, ydd2] = deriv_central(y, Fs);

figure, set(gcf, 'position', [339   476   913   391])
h1 = subplot(211);
plot(tt, y), grid on

h2 = subplot(212);
plot(tt, ydd)
hold on, plot(tt, ydd2, '--')
legend('analytical', 'numerical (ref)')
grid on
title(sprintf('%d Order', bord))


