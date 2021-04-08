%% Fusing IMU-OMC
%% Load dataset
clear, %clc

load('data/cable_driven_prosthesis4.mat', 'trial')
% 1) IMU extrinsic
% 2) IMU extrinsic (bugged)
% 3) Wrench (bugged)
% 4) Wrench

% load('data/cable_driven_prosthesis2.mat', 'trial', 'omc_files')
% 1) IMU extrinsic
% 2) Wrench
% 3) Active + Ext
% 4) Active
% 5) Walk

i = 1;
t = trial(i);

% trange = [70, 122];  % trial 1
% trange = [90, 100];  % trial 1
% trange = [5, 120];  % trial 4
trange = [14, 76];  % trial 1

% trange = [110, 120];  % trial 5
tout = t.ptime;

% shank ===========================
% mtime1 = t.mtime;
% quat1 = t.squat;
% trans1 = t.strans;
% mtrans1 = t.mstrans;
% time_imu1 = t.stime2 + tshift;
% w_imu1 = t.w2;
% a_imu1 = t.a2;
% calib = load('data/shank_imu2_calib_full.mat');

% foot ============================
mtime1 = t.mtime;
quat1 = t.fquat;
trans1 = t.ftrans;
mtrans1 = t.mftrans;
time_imu1 = t.stime1;
w_imu1 = t.w1;
a_imu1 = t.a1;
% calib = load('data/foot_imu1_calib_full.mat');
calib = load('data/calib_foot_day4.mat');

% viz ========================================================
Fs_omc = 1 / mean(diff(mtime1));

figure, 
% h1 = subplot(221); plot(t.mtime, w)
h1 = subplot(221); plot(mtime1, quat1, '.-')
ylabel('Quat [1]'), grid on, xline(trange(1)), xline(trange(2))
% h4 = subplot(222); plot(t.mtime, a)
h4 = subplot(222); plot(mtime1, trans1, '.-')
ylabel('Position [m]'), grid on, xline(trange(1)), xline(trange(2))
h2 = subplot(223); plot(time_imu1, w_imu1, '.-')
ylabel('\omega [rad/s]'), grid on, xline(trange(1)), xline(trange(2))
h3 = subplot(224); plot(time_imu1, a_imu1, '.-')
ylabel('a [m/s^2]'), grid on, xline(trange(1)), xline(trange(2))
linkaxes([h1, h2, h3, h4], 'x')

%% Select time range

% rows = find(mtime1 > trange(1) & mtime1 < trange(2));
% irows = find(time_imu1 > trange(1) & time_imu1 < trange(2));
% time0 = min(mtime1(rows(1)), time_imu1(irows(1)));
% 
% mtime = mtime1(rows,:) - time0;
% time_imu = time_imu1(irows,:) - time0;
% 
% quat = quat1(rows,:);
% trans = trans1(rows,:);
% mtrans = mtrans1(rows,:,:);
% 
% w_imu = w_imu1(irows,:);
% a_imu = a_imu1(irows,:);
% 
% if any(isnan(trans(:))) || any(isnan(quat(:)))
%     warning('Interpolating gaps')
%     trans = fillgaps(trans);
%     quat = quatnormalize(fillgaps(quat));
% end

mtime = mtime1;
quat = quat1;
trans = trans1;
mtrans = mtrans1;
time_imu = time_imu1;
w_imu = w_imu1;
a_imu = a_imu1;

start_opt_time = tic;
[cq, cs, cwbias, cabias, ushift, out_calib, quatf, transf] = fuse_OMC_IMU( ...
    mtime, quat, trans, mtrans, time_imu, w_imu, a_imu, calib, tout, trange);
toc(start_opt_time)

%%

[q1, q3, q2] = quat2angle(quat1, 'YZY'); 
qXZY = ([q1, q2, q3]);

[q2, q3, q1] = quat2angle(quat1, 'YZX'); 
qYZX = unwrap([q1, q2, q3]);

w = QTMParser.body_rates(quat1, trans1, 183, [5, 15]);

figure, 
h1 = subplot(321); plot(mtime1, qXZY, '.-'), grid on
h2 = subplot(323); plot(mtime1, qYZX, '.-'), grid on
h3 = subplot(325); plot(mtime1, quat1(:,1:4), '.-'), grid on
% h3 = subplot(313); plot(mtime, log_quat(quat), '.-'), grid on
h4 = subplot(322); plot(time_imu1, (calib.Tw \ (w_imu1'))', '.-'), grid on
h5 = subplot(324); plot(mtime1, w, '.-'), grid on

linkaxes([h1, h2, h3, h4, h5], 'x')
% xlim([56, 58])

%% Upsample 
dT = 30e-3;
tshift = -dT * ushift;

Fs = 1 / mean(diff(tout));
[w, wd, v, a] = QTMParser.body_rates(quat1, trans1, Fs_omc, [4, 13]);
[wf, wdf, vf, af] = QTMParser.body_rates(quatf, transf, Fs, [4, 5]);

figure, 
h1 = subplot(211); plot(mtime1, trans1, '.-')
hold on, plot(tout + tshift, transf, '-'), grid on
% xlim(trange)
h2 = subplot(212); plot(mtime1, quat1, '.-')
hold on, plot(tout + tshift, quatf, '-'), grid on
% xlim(trange)
linkaxes([h1, h2], 'x')

% figure, plot(mtime, w, '.-')
% hold on, plot(time, wf, '-'), grid on
% xlim([1, mtime(end)-1])
% figure, plot(mtime, a, '.-')
% hold on, plot(time, af, '-'), grid on
% xlim([1, mtime(end)-1])

%%
% trans2 = trans;
% quat2 = quat;
tt = mtime;
trans2 = interp1(tout, transf, mtime);
quat2 = interp1(tout, quatf, mtime);


% ri = calib.params.ri;
ri = out_calib.ri;
mtransh = zeros(size(mtrans));
for i = 1 : size(mtrans, 3)
    mtransh(:,:,i) = quatrotate(quatinv(quat2), ri(:,i)') + trans2;
end

% figure
% h1 = subplot(311); plot(mtime, squeeze(mtrans(:,1,:)))
% hold on, set(gca, 'ColorOrderIndex', 1), plot(mtime, squeeze(mtransh(:,1,:)), '--')
% h2 = subplot(312); plot(mtime, squeeze(mtrans(:,2,:)))
% hold on, set(gca, 'ColorOrderIndex', 1), plot(mtime, squeeze(mtransh(:,2,:)), '--')
% h3 = subplot(313); plot(mtime, squeeze(mtrans(:,3,:)))
% hold on, set(gca, 'ColorOrderIndex', 1), plot(mtime, squeeze(mtransh(:,3,:)), '--')
% linkaxes([h1, h2, h3], 'x')

figure
h1 = subplot(311); plot(tt, squeeze(mtrans(:,1,:) - mtransh(:,1,:)), '.-'), grid on
h2 = subplot(312); plot(tt, squeeze(mtrans(:,2,:) - mtransh(:,2,:)), '.-'), grid on
h3 = subplot(313); plot(tt, squeeze(mtrans(:,3,:) - mtransh(:,3,:)), '.-'), grid on
linkaxes([h1, h2, h3], 'x')
xlim(trange)
