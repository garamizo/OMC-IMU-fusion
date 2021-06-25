%% Extrinsic IMU-OMC calibration
clear, %clc
addpath(genpath("C:\Users\garamizo\Documents\GitHub\OMC_IMU_fusion"))

% load('C:\Users\garamizo\Documents\GitHub\cdprosthesis-desktop\MATLAB\data\cable_driven_prosthesis_calib_0510.mat'); i = 3;
% load('C:\Users\garamizo\Documents\GitHub\cdprosthesis-desktop\MATLAB\data\cable_driven_prosthesis_calib_0602.mat'); i = 1;
load('C:\Users\garamizo\Documents\GitHub\cdprosthesis-desktop\MATLAB\data\cable_driven_prosthesis4.mat'); i = 1;

t = trial(i);

%%

cal = Calibration_OMC_IMU(t.mtime, t.squat, t.strans, t.mstrans, t.packet.t_sensor_sync, ...
                                                                 t.packet.W1, ...
                                                                 t.packet.A1);
% cal = Calibration_OMC_IMU(t.mtime, t.fquat, t.ftrans, t.mftrans, t.packet.t_sensor_sync, ...
%                                                                  t.packet.W2, ...
%                                                                  t.packet.A2);
cal.trange_ = [7.5, 160];

% figure, cal.plot_data()

tic
[cq1, cs1, cwbias1, cabias1, Tw1, Ta1, r1, g131, tshift1, ri1, x] = cal.calibrate_LM();
[Tw1, r1]    
toc

toffset = tshift1*cal.dT_;

return

%% Debug each stage

[Tw, Ta, wbias, abias, r, g] = calibrate_lsq(cal);  
cal.fit_bspline_orientation(Tw);
cal.fit_bspline_translation(Ta, r, abias);  % TODO Fix 


%% Check if OMC has wrong labeled marker

% trans = t.ftrans;
% quat = t.fquat;
% mtrans = t.mftrans;
% 
% dlen = size(trans, 1);
% nmarkers = size(mtrans, 3);
% 
% mres = zeros(dlen, nmarkers);
% rms = zeros(3, nmarkers);
% for i = 1 : nmarkers
%     rms(:,i) = nanmedian(quatrotate(quat, mtrans(:,:,i) - trans));
%     
%     res = mtrans(:,:,i) - trans - quatrotate(quatinv(quat), rms(:,i)');
%     mres(:,i) = sqrt(sum(res.^2, 2));
% end
% rms
% 
% [bft, aft] = fir1(50, 0.5/(183/2));
% mres = filtfilt(bft, aft, fillmissing(mres, 'spline', 'EndValues', 'nearest'));
% 
% figure, 
% h1 = subplot(211); plot(t.mtime, quat)
% h2 = subplot(212); plot(t.mtime, mres)
% linkaxes([h1, h2], 'x')
