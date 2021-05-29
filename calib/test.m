clear, clc
addpath(genpath("C:\Users\garamizo\Documents\GitHub\OMC_IMU_fusion"))

load('C:\Users\garamizo\Documents\GitHub\cdprosthesis-desktop\MATLAB\data\cable_driven_prosthesis_calib_0510.mat')

%%

t = trial(3);

cal = Calibration_OMC_IMU(t.mtime, t.fquat, t.ftrans, t.mftrans, t.stime2, t.w2, t.a2);
% cal = Calibration_OMC_IMU(t.mtime, t.fquat, t.ftrans, t.mftrans, t.stime1, t.w1, t.a1);

cal.trange_ = [14, 155];

% figure
% cal.plot_data()

%%

[cq1, cs1, cwbias1, cabias1, Tw1, Ta1, r1, g131, tshift1, ri1, x] = cal.calibrate_LM();
[Tw1, Ta1, r1]

%%

cal = Calibration_OMC_IMU(t.mtime, t.fquat, t.ftrans, t.mftrans, t.stime2 - tshift1*cal.dT_, t.w2, t.a2);
cal.trange_ = [14, 155];
[cq1, cs1, cwbias1, cabias1, Tw1, Ta1, r1, g131, tshift1, ri1, x] = cal.calibrate_LM();
[Tw1, Ta1, r1]


%%

[Tw, Ta, wbias, abias, r, g] = calibrate_lsq(cal)

%%

cal.fit_bspline_orientation(Tw);

%%

cal.fit_bspline_translation(Tw, r, abias);

%%

%%
mtime = t.mtime;
quat = t.squat;

itime = t.stime2;
w_imu = t.w2;


[rotm_fcn, w_fcn] = Calibration_OMC_IMU.generate_EulerYZX_transformation();

[r2, r3, r1] = quat2angle(quat, 'YZX');
p = unwrap([r1, r2, r3]);
pd = deriv_sgolay(p, 183, [3, 7]);

w = interp1(itime, (w_imu - wbias') * inv(Tw)', mtime);
pdh = calib_pd(quat, w);
ph = unwrap(calib_p(quat));

figure, 
h1 = subplot(121); plot(p, '.')
set(gca, 'ColorOrderIndex', 1)
hold on, plot(ph)
h2 = subplot(122); plot(pd, '.-')
set(gca, 'ColorOrderIndex', 1)
hold on, plot(pdh)
linkaxes([h1, h2], 'x')

%%

dlen = length(mtime);
Rh = zeros(3, 3, dlen);
wh = zeros(dlen, 3);
for i = 1 : dlen
    Rh(:,:,i) = rotm_fcn(p(i,:)');
    wh(i,:) = w_fcn(p(i,:)', pd(i,:)');
end

figure, plot(mtime, wh, '.-')
set(gca, 'ColorOrderIndex', 1)
hold on
plot(t.mtime, w)

R = quat2rotm(quat);


%%
p = [0.1, 0.2, 0.3];
R = Rz(p(3)) * Rx(p(1)) * Ry(p(2))

Rh = angle2dcm(p(3), p(1), p(2), 'ZXY')'

%%

p = sym('p', [3, 1]);
pd = sym('pd', [3, 1]);

R1 = eye(3);
R2 = Rx(p(1)).' * Rz(p(3)).';
R3 = Ry(p(2)).';
S = [R1(1,:); R2(2,:); R3(3,:)];

w = S * pd

%%




R = quat2rotm_sym(quat);

%%

pt = str2sym('[pt1(t); pt2(t); pt3(t)]');
quatt = str2sym('[qt0(t), qt1(t), qt2(t), qt3(t)]');

eqs = diff(quat2rotm_sym(quatt), t) == diff(R, t);

p = sym('p', [3, 1]);
pd = sym('pd', [3, 1]);
quat

%%


% w = sym('w', [3, 1]);

% p = sym('p', [3, 1]);
% Rp = Ry(p(2)) * Rz(p(3)) * Rx(p(1))




matlabFunction(pd.', 'Vars', {quat, w.'}, 'File', 'calib_to_body_rates')

%%

syms x y z

Ry(y) * Rz(z) * Rx(x)

%%

trans = t.ftrans;
quat = t.fquat;
mtrans = t.mftrans;

dlen = size(trans, 1);
nmarkers = size(mtrans, 3);

mres = zeros(dlen, nmarkers);
rms = zeros(3, nmarkers);
for i = 1 : nmarkers
    rms(:,i) = nanmedian(quatrotate(quat, mtrans(:,:,i) - trans));
    
    res = mtrans(:,:,i) - trans - quatrotate(quatinv(quat), rms(:,i)');
    mres(:,i) = sqrt(sum(res.^2, 2));
end
rms

[bft, aft] = fir1(50, 0.5/(183/2));
mres = filtfilt(bft, aft, fillmissing(mres, 'spline', 'EndValues', 'nearest'));

figure, 
h1 = subplot(211); plot(t.mtime, quat)
h2 = subplot(212); plot(t.mtime, mres)
linkaxes([h1, h2], 'x')
