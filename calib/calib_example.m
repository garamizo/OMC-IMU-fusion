%% Extrinsic IMU-OMC calibration
%% Load dataset
clear, %clc
% load('C:\Users\garamizo\Documents\GitHub\mocap-utils\vibrating_platform\data\plate_extrinsic_calibration.mat', 'trial')
% load('data\cable_driven_prosthesis1.mat', 'trial')
% tshift_init = [-49.2260, 0];

% load('data\cable_driven_prosthesis2.mat', 'trial')
% tshift_init = [11.1765-16.7456, 7.52077-5.08493, 2.17377-5.07619, 6.8-1.02965+0, 1.85355-3.68149];

% load('data\cable_driven_prosthesis3.mat', 'trial')
% tshift_init = [11.1765-16.7456, 7.52077-5.08493, 2.17377-5.07619, 6.8-1.02965+0, 1.85355-3.68149];

% load('data\cable_driven_prosthesis4.mat', 'trial')
load('data\trial_plate_imu_extrinsics.mat', 'trial')



Fs_omc = 183;

i = 1;
t = trial(i);
% trange = [70, 122];  % trial 1
% trange = [90, 100];  % trial 1
% trange = [3, 122];  % trial 1
% trange = [14, 76];  % trial 1, shank and foot
trange = [11, 64];  % plate

% shank ===========================
% mtime = t.mtime;
% quat = t.squat;
% trans = t.strans;
% mtrans = t.mstrans;
% time_imu = t.stime2;
% w_imu = t.w2;
% a_imu = t.a2;

% foot ============================
% mtime = t.mtime;
% quat = t.fquat;
% trans = t.ftrans;
% mtrans = t.mftrans;
% time_imu = t.stime1;
% w_imu = t.w1;
% a_imu = t.a1;

% plate ============================
mtime = t.mtime;
quat = t.pquat;
trans = t.ptrans;
mtrans = t.mptrans;
time_imu = t.time_imu;
w_imu = t.w;
a_imu = t.a;

% viz ========================================================
[wd, w] = angular_rates(quat, Fs_omc, [3, 7]);
[v, ~, a] = deriv_sgolay(trans, Fs_omc, [3, 7]);            
figure, 
% h1 = subplot(221); plot(t.mtime, w)
h1 = subplot(221); plot(t.mtime, quat)
ylabel('Body orientation [deg]'), grid on, xline(trange(1)), xline(trange(2))
% h4 = subplot(222); plot(t.mtime, a)
h4 = subplot(222); plot(t.mtime, a)
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

uisave({'Tw', 'Ta', 'r', 'g13', 'tshift', 'ri', 'x', 'params', 'cq', 'cs', 'cwbias', 'cabias'}, 'calib_plate.mat')

%%
calib.Tw = Tw;

[q2, q3, q1] = quat2angle(quat, 'YZX'); 
qXZY = ([q1, q2, q3]);

[q2, q3, q1] = quat2angle(quat, 'YZX'); 
qYZX = unwrap([q1, q2, q3]);

w = QTMParser.body_rates(quat, trans, 183, [5, 15]);

figure, 
h1 = subplot(321); plot(mtime, qXZY, '.-'), grid on
h2 = subplot(323); plot(mtime, qYZX, '.-'), grid on
h3 = subplot(325); plot(mtime, quat(:,1:4), '.-'), grid on
% h3 = subplot(313); plot(mtime, log_quat(quat), '.-'), grid on
h4 = subplot(322); plot(time_imu, (calib.Tw \ (w_imu'))', '.-'), grid on
h5 = subplot(324); plot(mtime, w, '.-'), grid on

linkaxes([h1, h2, h3, h4, h5], 'x')
% xlim([56, 58])