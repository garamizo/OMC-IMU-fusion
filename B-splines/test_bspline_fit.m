%% Test fit
clear

% profile on
Fs = 800;
t = (0 : 1/Fs : 40)';
y = 10*polyval([1, -2, 1, 1], t / 200.1) + sin(2*pi*5 * t) + randn(length(t), 1) * 0.1;
y = [y, 25-y*2];

bord = 5;
dt = 5e-3;

tic
c = bspline_fit(bord, dt, t, y);
toc
% profile viewer

%% Test fit with derivative
clear

% profile on

bord = 6;
dt = 5e-3;

Fs = 800;
t = (0 : 1/Fs : 60)';

y = 10 * polyval([1, -2, 1, 1], t / t(end)) + sin(2*pi*20 * t);
y = [y, 10-2*y, detrend(y).^2];

t2 = t;
[yd, ~, ydd] = deriv_central(y, Fs);

t_ds = t(1:5:end);
y_ds = y(1:5:end,:);

y_ds = y_ds + randn(size(y_ds)) .* std(y_ds) * 5e-1;
yd = yd + randn(size(yd)) .* std(yd) * 1e-1 + [1, 2, 3];

% hold on, plot(t, y, '.-')
% figure, plot(t2, yd)

tic
% c = bspline_fit(bord, dt, t_ds, y_ds);
c = bspline_fit(bord, dt, t_ds, y_ds, t2, yd);
toc

hold on, plot(t, y, 'k--')
% profile viewer

%% Test fit with 2nd derivative
clear

bord = 6;
dt = 5e-3;

Fs = 800;
t = (0 : 1/Fs : 60)';

y = 10 * polyval([1, -2, 1, 1], t / t(end)) + sin(2*pi*20 * t);
y = [y, 10-2*y, detrend(y).^2];

t2 = t;
[yd, ~, ydd] = deriv_central(y, Fs);

t_ds = t(1:5:end);
y_ds = y(1:5:end,:);

y_ds = y_ds + randn(size(y_ds)) .* std(y_ds) * 5e-1;
ydd = ydd + randn(size(ydd)) .* std(ydd) * 0.3e-1;

% hold on, plot(t, y, '.-')
% figure, plot(t2, yd)

tic
% c = bspline_fit(bord, dt, t_ds, y_ds);
c = bspline_fit(bord, dt, t_ds, y_ds, t2, ydd, 'second');
toc

hold on, plot(t, y, 'k--')
