%% Extrinsic IMU-OMC calibration
%% Load dataset
clear, clc
% load('C:\Users\garamizo\Documents\GitHub\mocap-utils\vibrating_platform\data\plate_extrinsic_calibration.mat', 'trial')
load('data\cable_driven_prosthesis\cable_driven_prosthesis1.mat', 'trial')
tshift_init = [-49.2260, 0];

Fs_omc = 183;

i = 1;
t = trial(i);
% trange = [70, 122];  % trial 1
% trange = [90, 100];  % trial 1
trange = [3, 122];  % trial 1

% shank ===========================
% mtime = t.mtime;
% quat = t.squat;
% trans = t.strans;
% mtrans = t.mstrans;
% time_imu = t.stime + tshift_init(i);
% w_imu = t.w2;
% a_imu = t.a2;

% foot ============================
mtime = t.mtime;
quat = t.fquat;
trans = t.ftrans;
mtrans = t.mftrans;
time_imu = t.stime + tshift_init(i);
w_imu = t.w1;
a_imu = t.a1;

% viz ========================================================
[wd, w] = angular_rates(quat, Fs_omc, [3, 7]);
[v, ~, a] = deriv_sgolay(trans, Fs_omc, [3, 7]);            
figure, 
% h1 = subplot(221); plot(t.mtime, w)
h1 = subplot(221); plot(t.mtime, quat)
ylabel('Body orientation [deg]'), grid on, xline(trange(1)), xline(trange(2))
% h4 = subplot(222); plot(t.mtime, a)
h4 = subplot(222); plot(t.mtime, trans)
ylabel('Body velocity [m/s]'), grid on, xline(trange(1)), xline(trange(2))
h2 = subplot(223); plot(time_imu, w_imu)
ylabel('\omega [rad/s]'), grid on, xline(trange(1)), xline(trange(2))
h3 = subplot(224); plot(time_imu, a_imu)
ylabel('a [m/s^2]'), grid on, xline(trange(1)), xline(trange(2))
linkaxes([h1, h2, h3, h4], 'x')

%% Select time range

rows = find(mtime > trange(1) & mtime < trange(2));
irows = find(time_imu > trange(1) & time_imu < trange(2));
time0 = min(mtime(rows(1)), time_imu(irows(1)));
mtime = mtime(rows,:) - time0;
time_imu = time_imu(irows,:) - time0;

quat = quat(rows,:);
trans = trans(rows,:);
mtrans = mtrans(rows,:,:);

w_imu = w_imu(irows,:);
a_imu = a_imu(irows,:);

if any(isnan(trans(:))) || any(isnan(quat(:)))
    warning('Interpolating gaps')
    trans = fillgaps(trans);
    quat = quatnormalize(fillgaps(quat));
end

t0 = tic;
[cq, cs, cwbias, cabias, Tw, Ta, r, g13, tshift, ri, x, params] = calibrate_OMC_IMU( ...
    mtime, quat, trans, mtrans, time_imu, w_imu, a_imu);
toc(t0)
