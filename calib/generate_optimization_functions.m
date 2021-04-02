%% Calibrate
% Assume diagonal cov noise
clear, 

poly_order = 4;
poly_order_bias = 2;

% orientation representation ============================
% Rodrigues params
skew_func = @(p) [0, -p(3), p(2)
                  p(3), 0, -p(1)
                  -p(2), p(1), 0];
% C_func = @(p) eye(3) * (1 - 2*p.'*p) + 2*p*p.' + 2*sqrt(1 - p.'*p)*skew_func(p);
% S_func = @(p) 2*((eye(3) + skew_func(p)*skew_func(p))/sqrt(1 - p.'*p) - skew_func(p));
% S_func = @(p) eye(3) + 2 * (skew_func(p)*skew_func(p) - skew_func(p)) / (1 + p.'*p);  % Rodrigues, from book
% S_func = @(p) eye(3) + 2 * skew_func(p)^2 - 2*sqrt(1 - p.'*p)*skew_func(p);  % Euler params, from book
% to_body_rates = @(p, pd) S_func(p) * pd;

% % Euler angles YZX
C_func = @(p) Ry(p(2)) * Rz(p(3)) * Rx(p(1));
to_body_rates = @(p, pd) [pd(1); 0; 0] + Rx(p(1)).' * [0; 0; pd(3)] + (Rz(p(3)) * Rx(p(1))).' * [0; pd(2); 0];

% B-spline basis ============================================
syms t dt
% plate pose, translated to the IMU center
cq = sym('cq', [poly_order+1, 3]);
cs = sym('cs', [poly_order+1, 3]);
cwbias = sym('cwbias', [poly_order_bias+1, 3]);
cabias = sym('cabias', [poly_order_bias+1, 3]);

M = 1;
for k = 1 : poly_order
    M = ([M; zeros(1,k)] * ([diag(1:(k)), zeros(k,1)] + [zeros(k,1), diag((k-1):-1:0)]) + ...
         [zeros(1,k); M] * ([-eye(k), zeros(k,1)] + [zeros(k,1), eye(k)]) ) / k;
end
M = sym(M);

M_bias = 1;
for k = 1 : poly_order_bias
    M_bias = ([M_bias; zeros(1,k)] * ([diag(1:(k)), zeros(k,1)] + [zeros(k,1), diag((k-1):-1:0)]) + ...
         [zeros(1,k); M_bias] * ([-eye(k), zeros(k,1)] + [zeros(k,1), eye(k)]) ) / k;
end
M_bias = sym(M_bias);

% states ===========================================
% sensor pose
qp = ((t.^(0:poly_order)) * M * cq).';
qdp = diff(qp, t) / dt;
wp = to_body_rates(qp, qdp);

ss = ((t.^(0:poly_order)) * M * cs).';
vs = diff(ss, t) / dt;
as = diff(vs, t) / dt;

% biases
syms t_bias dt_bias
wbias = ((t_bias.^(0:poly_order_bias)) * M_bias * cwbias).';
wbiasd = diff(wbias, t_bias) / dt_bias;
abias = ((t_bias.^(0:poly_order_bias)) * M_bias * cabias).';
abiasd = diff(abias, t_bias) / dt_bias;

% extrinsic parameters ==============================
r = sym('r', [3, 1]);  % plate to sensor
Ta = sym('Ta', [3, 3]);
Tw = sym('Tw', [3, 3]);
g13 = sym('g', [2, 1]);
g = [g13(1); -sqrt(9.80136^2 + g13.'*g13); g13(2)];
ri = sym('ri', [3, 1]);  % plate to marker, plate frame
syms tshift

% intrinsic parameters  =============================
Nwb_inv = sym('Nwb_inv', [3, 1]);
Nab_inv = sym('Nab_inv', [3, 1]);
Nw_inv = sym('Nw', [3, 1]);
Na_inv = sym('Na', [3, 1]);

% measurements  =====================================
w_meas = sym('w_meas', [3, 1]);
a_meas = sym('a_meas', [3, 1]);

% Model ==============================================
% IMU
R = C_func(qp);  % plate orientation

w_imu = Tw * wp + wbias;
a_imu = Ta * R.' * (as - g) + abias;

sbody = ss - R * r;   % plate
rj = sbody + R * ri;  % marker
rj_shift = subs(rj, t, t+tshift);

% vars_all = [cq(:); cs(:); cwbias(:); cabias(:); Tw(:); Ta(:); r; g13; ri];
vars_imu = [cq(:); cs(:); cwbias(:); cabias(:); Tw(:); Ta(:); r; g13];
vars_bias = [cwbias(:); cabias(:)];
% vars_omc = [cq(:); cs(:); r];
vars_omc = [cq(:); cs(:); r; tshift; ri];
% vars_time = [cq(:); cs(:); cwbias(:); cabias(:); Tw(:); Ta(:); g13];
% vars_time = [cq(:); cs(:); cwbias(:); cabias(:)];

% cost = (w_imu - w_meas).' * diag(Nw_inv) * (w_imu - w_meas) + ...
%        (a_imu - a_meas).' * diag(Na_inv) * (a_imu - a_meas);
% grad = jacobian(cost, vars_time);
% hess = jacobian(grad, vars_time);
% hess_flat = reshape(hess, 1, []);

w_imu_jac = jacobian(w_imu, vars_imu);
a_imu_jac = jacobian(a_imu, vars_imu);
rj_jac = jacobian(rj_shift, vars_omc);
wbiasd_jac = jacobian(wbiasd, vars_bias);
abiasd_jac = jacobian(abiasd, vars_bias);


tic
matlabFunction(w_imu, a_imu, w_imu_jac, a_imu_jac, 'Vars', ...
    {vars_imu.', t, dt, t_bias, dt_bias}, ...
    'File', 'calib_cost_imu', 'Optimize', true)
matlabFunction(wbiasd, abiasd, wbiasd_jac, abiasd_jac, 'Vars', ...
    {vars_bias.', t_bias, dt_bias}, ...
    'File', 'calib_cost_imubias', 'Optimize', true)
matlabFunction(rj_shift, rj_jac, 'Vars', ...
    {vars_omc.', t, dt}, ...
    'File', 'calib_cost_omc', 'Optimize', true)

% matlabFunction(cost, grad, hess_flat, 'Vars', ...
%     {vars_time.', w_meas.', a_meas.', Nw_inv, Na_inv, t, dt, t_bias, dt_bias, Tw, Ta, g13}, ...
%     'File', 'fuse4_cost_imu', 'Optimize', true)
% matlabFunction(w_imu.', a_imu.', 'Vars', ...
%     {vars_time.', t, dt, t_bias, dt_bias, Tw, Ta, g13}, ...
%     'File', 'fuse4_eval_imu', 'Optimize', true)
toc


