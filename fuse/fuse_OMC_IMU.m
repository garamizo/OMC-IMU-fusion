function [cq1, cs1, cwbias1, cabias1, tshift1, out_calib, quatf, transf] = fuse_OMC_IMU( ...
    time, quat, trans, mtrans, time_imu, w_imu, a_imu, calib, tout, trange, dT)
% [cq(:); cs(:); cwbias(:); cabias(:); Tw(:); Ta(:); r; g13; tshift; ri(:)]
% The smoothing of high frequency data needs a filtering formulation of the
% continuous-time batch estimation, not the global optimization that we use
% for calibration. The difference is that motions during calibration are of
% low frequency, thus, a large time step can be used. 
% 
% addpath('C:\Users\garamizo\Documents\MATLAB\B-splines')

% assert(~any(isnan(trans(:))), 'Cannot have gap')
if any(isnan(trans(:)))
    warning('Interpolating gaps')
    quat = quatnormalize(fillgaps(quat, 3, 1));
    trans = fillgaps(trans, 3, 1);
end

%% Fusing parameters

bord = calib.params.bord;
bord_bias = calib.params.bord_bias;
% dT = calib.params.dT;
% dT = 5e-3;
dT_bias = calib.params.dT_bias;


% Intrinsic params --------------------------------
Nr = calib.params.Nr;  % from QTM calibration, at 183 Hz
w_white = calib.params.w_white;
a_white = calib.params.a_white;
w_walk = calib.params.w_walk / 1;
a_walk = calib.params.a_walk / 1;

Fs_omc = 1 / median(diff(time));
Fs_imu = 1 / median(diff(time_imu));
assert(abs(Fs_omc - calib.params.Fs_omc) < 1.0, 'OMC sampling rate is too different from calib')
% assert(abs(Fs_imu - calib.params.Fs_imu) < 1.0, 'IMU sampling rate is too different from calib')

Nw = (w_white * sqrt(Fs_imu)).^2;
Na = (a_white * sqrt(Fs_imu)).^2;
Nwb = w_walk.^2;
Nab = a_walk.^2;

%% Initialize extrinsic params

Tw = calib.Tw;
Ta = calib.Ta;
wbias = median(calib.cwbias)';
abias = median(calib.cabias)';
r = calib.r;
g13 = calib.g13;
ri = calib.ri;
tshift = 0;

nmarkers = size(mtrans, 3);

%% Initialize motion params
t0 = max([time(1), tout(1), trange(1)]);
tend = min([time(end), tout(end), trange(end)]);

% crop data ==========================
rows = time >= t0 & time <= tend;
quat = quat(rows,:);
trans = trans(rows,:);
mtrans = mtrans(rows,:,:);
time = time(rows);

rows = time_imu >= time(1) & time_imu < time(end);
a_imu = a_imu(rows,:);
w_imu = w_imu(rows,:);
time_imu = time_imu(rows);
% =====================================

time_norm = (time - t0) / dT;
time2_norm = (time_imu - t0) / dT;
time_norm_bias = (time_imu - t0) / dT_bias;

nknot = ceil((tend - t0) / dT + 1e-10) + bord - 1;
nknot_bias = ceil((tend - t0) / dT_bias + 1e-10) + bord_bias - 1;

% Init orientation ================================
[cq, ~] = bspline_fit_orientation(time_norm, quat, time2_norm, w_imu, Tw, bord, dT, nknot);
% cq = fillgaps(cq, 3, 1);

% Init translation ================================
[cs, ~] = bspline_fit_translation(time_norm, quat, trans, time2_norm, a_imu, Ta, r, abias, bord, dT, nknot);
% cs = fillgaps(cs, 3, 1);

% Init biases =====================================
cwbias = repmat(wbias', [nknot_bias, 1]);
cabias = repmat(abias', [nknot_bias, 1]);

%% Optimize

x0 = [cq(:); cs(:); cwbias(:); cabias(:); tshift];
% x0 = [cq(:); cs(:); cwbias(:); cabias(:); Tw(:); Ta(:); r; g13];

iters = 30;
X = zeros(iters, length(x0));
Fval = zeros(iters, 1);
x = x0;
for i = 1 : iters
    tic
    [dx, fval] = solve_gauss_newton_sparse_fast(x, time_norm, mtrans, time2_norm, w_imu, a_imu, time_norm_bias, ...
        bord, bord_bias, dT, dT_bias, nknot, nknot_bias, Nw, Na, Nr, Nwb, Nab, Tw, Ta, r, g13, ri);
    titer = toc;

    if any(isnan(dx))
        warning('Went NAN')
        break
    end
%     dx(nknot*6 + nknot_bias*6 + [1:23, 24+(1:nmarkers*3)]) = 0;
    x = x + dx;
    
    X(i,:) = x;
    Fval(i) = fval;
%     fprintf('Iter %2d)\tCost: %.7e\tDuration: %.3f\n', i, fval, titer)
    
    if i > 1 && (Fval(i-1) - fval)/Fval(i-1) < 1/10000
%         fprintf('Converged!\n')
        break
    end
end

cq1 = reshape(x(1:3*nknot), [], 3);
cs1 = reshape(x(3*nknot + (1:3*nknot)), [], 3);
cwbias1 = reshape(x(6*nknot + (1:3*nknot_bias)), [], 3);
cabias1 = reshape(x(6*nknot + 3*nknot_bias + (1:3*nknot_bias)), [], 3);
% Tw1 = reshape(x(6*nknot + 6*nknot_bias + (1:9)), 3, 3);
% Ta1 = reshape(x(6*nknot + 6*nknot_bias + (10:18)), 3, 3);
% r1 = reshape(x(6*nknot + 6*nknot_bias + (19:21)), 3, 1);
% g131 = reshape(x(6*nknot + 6*nknot_bias + (22:23)), 2, 1);
tshift1 = reshape(x(6*nknot + 6*nknot_bias + 1), 1, 1);
% ri1 = reshape(x(6*nknot + 6*nknot_bias + 24 + (1:nmarkers*3)), 3, []);

% this should not change much from calib
out_calib = [];%struct('Tw', Tw1, 'Ta', Ta1, 'r', r1, 'g13', g131, 'tshift', tshift1, 'ri', ri1, ...
%     'params', struct('bord', bord, 'bord_bias', bord_bias, 'dT', dT, 'dT_bias', dT_bias, ...
%     'Fs_omc', Fs_omc, 'Fs_imu', Fs_imu));

% solve_gauss_newton_sparse_fast(x, time_norm, mtrans, time2_norm, w_imu, a_imu, time_norm_bias, ...
%         bord, bord_bias, dT, dT_bias, nknot, nknot_bias, Nw, Na, Nr, Nwb, Nab, Tw, Ta, r, g13, ri);
    
% eval on desired time ============================================
time_norm = (tout - t0) / dT + tshift;

% Euler YZX
time_norm_shift = time_norm;
time_norm_shift(time_norm_shift < 0) = 0;
time_norm_shift(time_norm_shift > length(cq1)-bord) = length(cq1)-bord;

qf = bspline_eval(cq, time_norm_shift, calib.params.bord);
sf = bspline_eval(cs, time_norm_shift, calib.params.bord);

quatf = angle2quat(qf(:,2), qf(:,3), qf(:,1), 'YZX');
transf = sf - quatrotate(quatinv(quatf), calib.r');

%%

% figure, 
% h1 = subplot(211); plot(time, trans, '.-')
% hold on, plot(tout, transf, '-'), grid on
% % xlim(trange)
% h2 = subplot(212); plot(time, quat, '.-')
% hold on, plot(tout, quatf, '-'), grid on
% % xlim(trange)
% linkaxes([h1, h2], 'x')

end

function A = generate_predictor_matrix(bord, time_norm, nknot, order)
    % time_norm := (time - t0) / dT
    % A := length(time_norm) x nknot
    
    % minimum number of knots needed
    nsamples = length(time_norm);
    
    % sparse representation of predictor matrix A
    nonzero_len = ceil(bord / mean(diff(time_norm)));
    aa = zeros(nonzero_len, nknot);
    jj = zeros(nonzero_len, nknot);
    ii = zeros(nonzero_len, nknot);
    
    % single B-spline basis
    ci = zeros(nknot, 1);
    ci(1) = 1;
    
    for j = 1 : nknot

        k0 = find(time_norm >= j, 1, 'first');
        
        if isempty(k0)  % border cases
            k0 = nsamples + 1;
        elseif k0 - nonzero_len < 1
            k0 = nonzero_len + 1;
        end
        
        nonzero_rows = k0 + (-nonzero_len : -1);
        ii(:,j) = nonzero_rows;
        jj(:,j) = j;
        if order == 1
            aa(:,j) = bspline_eval(ci, time_norm(nonzero_rows), bord);
        elseif order == 2
            [~, aa(:,j)] = bspline_eval(ci, time_norm(nonzero_rows), bord);
        elseif order == 3
            [~, ~, aa(:,j)] = bspline_eval(ci, time_norm(nonzero_rows), bord);
        end
        
        ci = circshift(ci, 1);
    end
    
    valid = ii(:) > 0 & jj(:) > 0;
    A = sparse(ii(valid), jj(valid), aa(valid), nsamples, nknot);
end


function [cq, wbias] = bspline_fit_orientation(...
    kk, quat, kk2, wimu, Tw, bord, dt, nknot)

    [q2, q3, q1] = quat2angle(quat, 'YZX');
    q = unwrap([q1, q2, q3]);
    w = wimu * inv(Tw)';
    
    % to angular velocity
    S = @(q) [1, cos(q(1))*sin(q(3)), 0
              0, cos(q(1))*cos(q(3)), sin(q(1))
              0, -sin(q(1)), cos(q(1))];
    qd = zeros(size(w));
    q_interp = interp1(kk, q, kk2);
    for i = 1 : size(w, 1)
        qd(i,:) = (S(q_interp(i,:)) \ w(i,:)')';
    end
    
    % prior =================================================
%     A0 = speye(nknot);
%     b0 = zeros(nknot, 3);

    % position ==============================================
    A1 = generate_predictor_matrix(bord, kk, nknot, 1);
    b1 = q;

    % velocity ================================================
    A2 = generate_predictor_matrix(bord, kk2, nknot, 2) / dt;
    b2 = qd;

    % combine ==================================================
    A = [A1, zeros(size(A1,1), 1)
          A2,  ones(size(A2,1), 1)];
    b = [b1; b2];
    W = blkdiag(speye(length(kk))/2.5e-5, speye(length(kk2))/0.1e-3);  % weight
    cq_bias = (A' * W * A) \ (A' * W * b);
%     cq_bias = [A1 \ b1; zeros(1, 3)];
    
    cq = cq_bias(1:end-1,:);
    wbias = cq_bias(end,:) * Tw';
    
    if nargout == 0
        bh = A * cq_bias;
        b1h = bh(1:size(b1,1),:);
        b2h = bh(size(b1,1)+1:end,:);

        fit = 1 - [goodnessOfFit(b1h, b1, 'NMSE'), goodnessOfFit(b2h, b2, 'NMSE')];
        fprintf('Orientation fitness (NMSE): [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f]\n', ...
            fit(1), fit(2), fit(3), fit(4), fit(5), fit(6))
    
        figure
        h1 = subplot(211); plot(kk, b1, kk, b1h)
        legend([compose("meas %d", 1:size(b1,2)), compose("fit %d", 1:size(b1,2))])
        h2 = subplot(212); plot(kk2, b2, kk2, b2h)
        legend([compose("meas %d", 1:size(b1,2)), compose("fit %d", 1:size(b1,2))])
        linkaxes([h1, h2], 'x')
    end
end

function [cs, grav] = bspline_fit_translation(...
    kk, quat, trans, kk2, aimu, Ta, r, abias, bord, dt, nknot)
    % fit cs as the IMU position

    trans_imu = trans + quatrotate(quatinv(quat), r');
    quat_interp = quatnormalize(interp1(kk, quat, kk2, 'pchip', 'extrap'));
    as = quatrotate(quatinv(quat_interp), (aimu - abias') * inv(Ta)');

    % position ==============================================
    A1 = generate_predictor_matrix(bord, kk, nknot, 1);
    b1 = trans_imu;

    % acceleration ================================================
    A2 = generate_predictor_matrix(bord, kk2, nknot, 3) / dt^2;
    b2 = as;

    % combine ==================================================
    A = [A1, zeros(size(A1,1), 1)
          A2,  ones(size(A2,1), 1)];
    b = [b1; b2];
    W = blkdiag(speye(length(kk)), 0.1*speye(length(kk2)));  % weight
    cs_grav = (A' * W * A) \ (A' * W * b);
    
%     cs_grav = [A1 \ b1; 0, -9.81, 0];
    
    cs = cs_grav(1:end-1,:);
    grav = cs_grav(end,:);
    
    if nargout == 0
        bh = A * cs_grav;
        b1h = bh(1:size(b1,1),:);
        b2h = bh(size(b1,1)+1:end,:);

        fit = 1 - [goodnessOfFit(b1h, b1, 'NMSE'), goodnessOfFit(b2h, b2, 'NMSE')];
        fprintf('Orientation fitness (NMSE): [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f]\n', ...
            fit(1), fit(2), fit(3), fit(4), fit(5), fit(6))
    
        figure
        h1 = subplot(211); plot(kk, b1, kk, b1h)
        legend([compose("meas %d", 1:size(b1,2)), compose("fit %d", 1:size(b1,2))])
        h2 = subplot(212); plot(kk2, b2, kk2, b2h)
        legend([compose("meas %d", 1:size(b1,2)), compose("fit %d", 1:size(b1,2))])
        linkaxes([h1, h2], 'x')
    end
end


function [dx, cost] = solve_gauss_newton_sparse_fast(x, time_norm, mtrans, time2_norm, w_imu, a_imu, time_norm_bias, ...
    bord, bord_bias, dT, dT_bias, nknot, nknot_bias, Nw, Na, Nr, Nwb, Nab, Tw, Ta, r, g13, ri)

    integ_npoints = 5;
    
    nmarkers = size(mtrans, 3);
    wh = zeros(length(time2_norm), 3);
    ah = zeros(length(time2_norm), 3);
    mtransh = zeros(length(time_norm), 3, nmarkers);

    % [cq(:); cs(:); cwbias(:); cabias(:); tshift]
    index = reshape((1 : bord)' + (0 : 2)*nknot, [], 1);
    index_bias = reshape((1 : bord_bias)' + (0 : 2)*nknot_bias, [], 1);
    index_omcparam = 1;  % tshift
    
    time_norm = time_norm + x(nknot*6 + nknot_bias*6 + 1);
    x(nknot*6 + nknot_bias*6 + 1) = 0;

    Nw_inv = diag(1./Nw);
    Na_inv = diag(1./Na);
    Nr_inv = diag(1./Nr);
    Nwb_inv = diag(1./Nwb);
    Nab_inv = diag(1./Nab);
    
    A_sparse = zeros(length(x) + ...
        (length(index)*2 + length(index_bias)*2)^2 * (nknot + nknot_bias) + ...
        (length(index)*2 + length(index_omcparam))^2 * nknot * nmarkers + ...
        (length(index_bias)*2)^2 * (nknot_bias-bord_bias) * integ_npoints, 1);
    i_sparse = ones(length(x) + ...
        (length(index)*2 + length(index_bias)*2)^2 * (nknot + nknot_bias) + ...
        (length(index)*2 + length(index_omcparam))^2 * nknot * nmarkers + ...
        (length(index_bias)*2)^2 * (nknot_bias-bord_bias) * integ_npoints, 1);
    j_sparse = i_sparse;

    b = zeros(length(x), 1);
    cost = 0;
    
    % Prior ============================================================
    % [cq(:); cs(:); cwbias(:); cabias(:); tshift]

    % make prior (regularization)
    c_confidence = 1e-5;
    buf = 1 : length(x);
    A_sparse(buf,:) = [c_confidence * ones(nknot*3,1)/(0.2*pi/180)^2
                    c_confidence * ones(nknot*3,1)/(1e-3)^2
                    c_confidence * ones(nknot_bias*3,1)/(1e-2)^2
                    c_confidence * ones(nknot_bias*3,1)/(1e-5)^2
                    c_confidence * ones(1,1)/(10e-3 / dT)^2];
    i_sparse(buf) = buf;
    j_sparse(buf) = buf;
    buf = buf + length(x);
    

    % IMU ==============================================================
    buflen = (length(index)*2 + length(index_bias)*2)^2;
    buf = buf(1) - 1 + (1 : buflen);
    jac = zeros(length(index)*2 + length(index_bias)*2);
    k0_last = floor(time2_norm(1)); 
    k0_bias_last = floor(time_norm_bias(1));
    for k = 1 : length(time2_norm)

        k0 = floor(time2_norm(k));
        k0_bias = floor(time_norm_bias(k));
        
        if k0_last ~= k0 || k0_bias_last ~= k0_bias
            % save progress
            i_sparse(buf) = reshape(repmat(rows, [1, length(rows)]), [], 1);
            j_sparse(buf) = reshape(repmat(rows', [length(rows), 1]), [], 1);
            A_sparse(buf) = reshape(jac, [], 1);
            buf = buf + buflen;
            jac = jac * 0;
        end

        rows = [k0 + index
                k0 + index + nknot*3
                k0_bias + index_bias + nknot*6
                k0_bias + index_bias + nknot*6 + nknot_bias*3];
        vars = x(rows);

        [wh(k,:), ah(k,:), wjac, ajac] = fuse_cost_imu(...
            vars.', mod(time2_norm(k), 1), dT, mod(time_norm_bias(k), 1), dT_bias, Tw, Ta, g13);

        jac = jac + wjac' * Nw_inv * wjac ...
                  + ajac' * Na_inv * ajac;
        b(rows) = b(rows) - wjac' * Nw_inv * (w_imu(k,:) - wh(k,:))' ...
                          - ajac' * Na_inv * (a_imu(k,:) - ah(k,:))';
        cost = cost + (w_imu(k,:) - wh(k,:)) * Nw_inv * (w_imu(k,:) - wh(k,:))' / 2 + ...
                      (a_imu(k,:) - ah(k,:)) * Na_inv * (a_imu(k,:) - ah(k,:))' / 2;
        
        k0_last = k0;
        k0_bias_last = k0_bias;
    end
    i_sparse(buf) = reshape(repmat(rows, [1, length(rows)]), [], 1);
    j_sparse(buf) = reshape(repmat(rows', [length(rows), 1]), [], 1);
    A_sparse(buf) = reshape(jac, [], 1);
    buf = buf + buflen;

   % OMC ===============================================================
    buflen = (length(index)*2 + length(index_omcparam))^2;
    buf = buf(1) - 1 + (1 : buflen);
    jac = zeros(length(index)*2 + length(index_omcparam));
    k0_last = floor(time_norm(1)); 
    rows = [];
    for i = 1 : nmarkers
        for k = 1 : length(time_norm)

            if any(isnan(mtrans(k,:,i)))
                continue
            end
            
            k0 = floor(time_norm(k));
            if k0 < 0 || k0 >= nknot-bord
                continue
            end
            
            if k0_last ~= k0 && ~isempty(rows)
                i_sparse(buf) = reshape(repmat(rows, [1, length(rows)]), [], 1);
                j_sparse(buf) = reshape(repmat(rows', [length(rows), 1]), [], 1);
                A_sparse(buf) = reshape(jac, [], 1);
                buf = buf + buflen;
                jac = jac * 0;
            end
            
            rows = [k0 + index
                    k0 + index + nknot*3
                    nknot*6 + nknot_bias*6 + index_omcparam];
            vars = x(rows);

            [mtransh(k,:,i), rjac] = fuse_cost_omc(...
                vars.', mod(time_norm(k), 1), dT, r, ri(:,i));
            jac = jac + rjac' * Nr_inv * rjac;
            b(rows) = b(rows) - rjac' * Nr_inv * (mtrans(k,:,i) - mtransh(k,:,i))';
            cost = cost + (mtrans(k,:,i) - mtransh(k,:,i)) * Nr_inv * (mtrans(k,:,i) - mtransh(k,:,i))' / 2;
            
            k0_last = k0;
        end
    end
    i_sparse(buf) = reshape(repmat(rows, [1, length(rows)]), [], 1);
    j_sparse(buf) = reshape(repmat(rows', [length(rows), 1]), [], 1);
    A_sparse(buf) = reshape(jac, [], 1);
    buf = buf + buflen;

    % Bias =============================================================
    % Approximate integral with Euler

    buflen = (length(index_bias)*2)^2;
    buf = buf(1) - 1 + (1 : buflen);
    for k = 1 : nknot_bias-bord_bias

        k0_bias = k - 1;
        rows = [k0_bias + index_bias + nknot*6
                k0_bias + index_bias + nknot*6 + nknot_bias*3];
        vars = x(rows);
        for kk = 1 : integ_npoints
            t = (kk - 1) / integ_npoints;
            [wbiasd, abiasd, wbiasd_jac, abiasd_jac] = fuse_cost_imubias( ...
                vars', t, dT_bias);

            i_sparse(buf) = reshape(repmat(rows, [1, length(rows)]), [], 1);
            j_sparse(buf) = reshape(repmat(rows', [length(rows), 1]), [], 1);
            A_sparse(buf) = reshape(wbiasd_jac' * Nwb_inv * wbiasd_jac * dT_bias / integ_npoints + ...
                abiasd_jac' * Nab_inv * abiasd_jac * dT_bias / integ_npoints, [], 1);
            buf = buf + buflen;
            
            b(rows) = b(rows) + wbiasd_jac' * Nwb_inv * wbiasd * dT_bias / integ_npoints ...
                              + abiasd_jac' * Nab_inv * abiasd * dT_bias / integ_npoints;
            cost = cost + wbiasd' * Nwb_inv * wbiasd * dT_bias / integ_npoints + ...
                          abiasd' * Nab_inv * abiasd * dT_bias / integ_npoints;
        end
    end
           
    if nargout > 0
        A = sparse(i_sparse, j_sparse, A_sparse, length(x), length(x));
        dx = -A \ b;  % fails with singular matrices
%         dx = -pinv(A) * b;  % does not support sparse
%         dx = lsqminnorm(A, -b);  % too slow
    
    else
        figure, 
        h1 = subplot(221); plot(time2_norm, w_imu, time2_norm, wh, '--')
        h2 = subplot(223); plot(time2_norm, a_imu, time2_norm, ah, '--')
        
        h3 = subplot(322); plot(time_norm, squeeze(mtrans(:,1,:)))
        hold on, set(gca, 'ColorOrderIndex', 1), plot(time_norm, squeeze(mtransh(:,1,:)), '--')
        h4 = subplot(324); plot(time_norm, squeeze(mtrans(:,2,:)))
        hold on, set(gca, 'ColorOrderIndex', 1), plot(time_norm, squeeze(mtransh(:,2,:)), '--')
        h5 = subplot(326); plot(time_norm, squeeze(mtrans(:,3,:)))
        hold on, set(gca, 'ColorOrderIndex', 1), plot(time_norm, squeeze(mtransh(:,3,:)), '--')
        linkaxes([h1, h2, h3, h4, h5], 'x')
    end
end

