%% Fusing IMU-OMC
%% Load dataset
clear, clc
load('data/cable_driven_prosthesis2.mat', 'trial', 'omc_files')
% 1) IMU extrinsic
% 2) Wrench
% 3) Active + Ext
% 4) Active
% 5) Walk
tshifts = [11.1765-16.7456, 7.52077-5.08493, 2.17377-5.07619, 6.8-1.02965+0, 1.85355-3.68149];

i = 1;
t = trial(i);
tshift = tshifts(i);

% trange = [70, 122];  % trial 1
% trange = [90, 100];  % trial 1
trange = [50, 90];  % trial 4
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
time_imu1 = t.stime1 + tshift;
w_imu1 = t.w1;
a_imu1 = t.a1;
calib = load('data/foot_imu1_calib_full.mat');

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
[cq, cs, cwbias, cabias, tshift, out_calib, quatf, transf] = fuse_OMC_IMU( ...
    mtime, quat, trans, mtrans, time_imu, w_imu, a_imu, calib, tout, trange);
toc(start_opt_time)

%% Upsample 

time_norm = (tout - time0) / out_calib.params.dT + tshift;
time_norm_bias = (tout - time0) / out_calib.params.dT_bias;

% Euler YZX
time_norm_shift = time_norm;
time_norm_shift(time_norm_shift < 0) = 0;
time_norm_shift(time_norm_shift > length(cq)-out_calib.params.bord) = length(cq)-out_calib.params.bord;

time_norm_shift_bias = time_norm_bias;
time_norm_shift_bias(time_norm_shift_bias < 0) = 0;
time_norm_shift_bias(time_norm_shift_bias > length(cwbias)-out_calib.params.bord_bias) = length(cwbias)-out_calib.params.bord_bias;

[qf, qdf, qddf] = bspline_eval(cq, time_norm_shift, out_calib.params.bord);
[sf, sdf, sddf] = bspline_eval(cs, time_norm_shift, out_calib.params.bord);
wbias = bspline_eval(cwbias, time_norm_shift_bias, out_calib.params.bord_bias);

quatf = angle2quat(qf(:,2), qf(:,3), qf(:,1), 'YZX');
transf = sf - quatrotate(quatinv(quatf), calib.r');

Fs = 1 / mean(diff(tout));
[w, wd, v, a] = QTMParser.body_rates(quat1, trans1, Fs_omc, [4, 13]);
[wf, wdf, vf, af] = QTMParser.body_rates(quatf, transf, Fs, [4, 5]);

figure, plot(mtime1, trans1, '.-')
hold on, plot(tout, transf, '-'), grid on
xlim(trange)
figure, plot(mtime1, quat1, '.-')
hold on, plot(tout, -quatf, '-'), grid on
xlim(trange)


% figure, plot(mtime, w, '.-')
% hold on, plot(time, wf, '-'), grid on
% xlim([1, mtime(end)-1])
% figure, plot(mtime, a, '.-')
% hold on, plot(time, af, '-'), grid on
% xlim([1, mtime(end)-1])


