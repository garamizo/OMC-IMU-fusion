classdef Calibration_OMC_IMU < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mtime_
        quat_
        trans_
        mtrans_
        itime_
        w_
        a_
        
        trange_
        
        % params =========================================
        bord_ = 4 + 1
        bord_bias_ = 2 + 1
        dT_ = 15e-3
        dT_bias_ = 100e-3
        integ_npoints_ = 5
        
        sgolay_params_ = [4, 15]  % for LSQ method
        
        % Intrinsic params --------------------------------
        Nr_ = (1*ones(1, 3) * 0.4e-3).^2;  % from QTM calibration, at 183 Hz
        w_white = [0.3573, 0.2735, 0.2739] * 1e-3;
        a_white = [0.6213, 0.7069, 0.9750] * 1e-3;
        w_walk = [1.2448, 2.2278, 0.1894] * 1e-5;
        a_walk = [0.9647, 1.7339, 2.4046] * 1e-6;
    end
    
    methods
        function obj = Calibration_OMC_IMU(mtime, quat, trans, mtrans, itime, w, a)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            
            itime_fix = obj.sync_quat_w(mtime, quat, itime, w);
            
            
            obj.set_data(mtime, quat, trans, mtrans, itime_fix, w, a);
            obj.trange_ = [max(mtime(1), itime(1)), min(mtime(end), itime(end))];
        end
        
        function plot_data(obj)     
            [mtime, quat, trans, mtrans, time_imu, w_imu, a_imu] = obj.get_raw_data();
            trange = obj.trange_;
            
            h1 = subplot(231); plot(mtime, quat)
            ylabel('Body quat'), grid on, xline(trange(1)), xline(trange(end))
            h4 = subplot(234); plot(mtime, trans)
            ylabel('Body trans'), grid on, xline(trange(1)), xline(trange(end))
            
            h2 = subplot(233); plot(time_imu, w_imu)
            ylabel('\omega [rad/s]'), grid on, xline(trange(1)), xline(trange(end))
            h3 = subplot(236); plot(time_imu, a_imu)
            ylabel('a [m/s^2]'), grid on, xline(trange(1)), xline(trange(end))
            
            h7 = subplot(332); plot(mtime, squeeze(mtrans(:,1,:)))
            ylabel('Body marker x'), grid on, xline(trange(1)), xline(trange(end))
            h5 = subplot(335); plot(mtime, squeeze(mtrans(:,2,:)))
            ylabel('Body marker y'), grid on, xline(trange(1)), xline(trange(end))
            h6 = subplot(338); plot(mtime, squeeze(mtrans(:,3,:)))
            ylabel('Body marker z'), grid on, xline(trange(1)), xline(trange(end))
            
            linkaxes([h1, h2, h3, h4, h5, h6, h7], 'x')
        end

        function [cq1, cs1, cwbias1, cabias1, T1, r1, g131, tshift1, ri1, x] = calibrate_LM(obj)
            
            [mtime, quat, trans, mtrans, itime, w_imu, a_imu] = obj.get_data();
            [time_norm, time2_norm, time_norm_bias, nknot, nknot_bias, integ_npoints, ...
                Nw, Na, Nr, Nwb, Nab] = obj.get_bspline_vars();
            bord = obj.bord_;
            bord_bias = obj.bord_bias_;
            dT = obj.dT_;
            dT_bias = obj.dT_bias_;
            trange = obj.trange_;
            
            % Estimate initial approximate solution ================================
            [Tw, Ta, wbias, abias, r, g] = obj.calibrate_lsq();
            T = quat2rotm(rotm2quat(0.5*Tw + 0.5*Ta));
            [cq, wbias2] = obj.fit_bspline_orientation(T);
            [cs, g2] = obj.fit_bspline_translation(T, r, abias);
            cwbias = repmat(wbias', [nknot_bias, 1]);
            cabias = repmat(abias', [nknot_bias, 1]);
            
            qq = [atan2(-T(2,3), T(2,2))
                  atan2(-T(3,1), T(1,1))
                  asin(T(2,1))];
            
            g13 = g([1,3]);
            
            nmarkers = size(mtrans, 3);
            ri = zeros(3, nmarkers);
            for j = 1 : nmarkers
                ri(:,j) = nanmedian(quatrotate(quat, mtrans(:,:,j) - trans));
            end
            tshift = 0;
            
%             qq, r, g13
            
            x0 = [cq(:); cs(:); cwbias(:); cabias(:); qq; r; g13; tshift; ri(:)];
            
            % Setup variables ===========================================
            % [cq(:); cs(:); cwbias(:); cabias(:); qq; r; g13; tshift; ri(:)]
            index = reshape((1 : bord)' + (0 : 2)*nknot, [], 1);
            index_bias = reshape((1 : bord_bias)' + (0 : 2)*nknot_bias, [], 1);
            index_imuparam = (1:8)';  % [qq; r; g13]
            index_omcparam_fix = [4:6, 9]';  % r, tshift
            index_omcparam_var = (10:12)';  % ri (first of)

            Nw_inv = diag(1./Nw);
            Na_inv = diag(1./Na);
            Nr_inv = diag(1./Nr);
            Nwb_inv = diag(1./Nwb);
            Nab_inv = diag(1./Nab);
            
            % Iterate ==================================================
            max_iters = 10;
            X = zeros(length(x0), max_iters);
            Fval = Inf(max_iters, 1);
            Lambda = zeros(max_iters, 1);
            x = x0;
 
            lambda = 1e-8;
            
            for iter = 1 : max_iters
                [Ai, bi, fval] = LM_iteration(x);
                
                if mod(iter, 1) == 0
                    fprintf("%d) cost = %.5e,\t Lambda = %.2e", iter, fval, lambda)
                end

                if fval < min(Fval)
                    A = Ai;
                    b = bi;
                    xopt = x;
                    
                    if iter > 1 && iter - 1 == iopt
                        lambda = lambda / 3;
                    end
                    iopt = iter;
                    fprintf(" *\n");
                    
                else
                    lambda = lambda * 3 * 10;
                    fprintf("\n");

                    if iter - iopt > 5
                        break
                    end
                end
                
                X(:,iter) = x;
                Fval(iter) = fval;
                Lambda(iter) = lambda;
                
                valid_rows = any(A ~= 0, 1);
                dx = zeros(length(x), 1);
                dx(valid_rows) = -(A(valid_rows,valid_rows) + lambda*diag(diag(A(valid_rows,valid_rows)))) \ b(valid_rows);

                x = xopt + dx;
                lambda = max(lambda / 3, 1e-20);
            end
            x = xopt;
            
            cq1 = reshape(x(1:3*nknot), [], 3);
            cs1 = reshape(x(3*nknot + (1:3*nknot)), [], 3);
            cwbias1 = reshape(x(6*nknot + (1:3*nknot_bias)), [], 3);
            cabias1 = reshape(x(6*nknot + 3*nknot_bias + (1:3*nknot_bias)), [], 3);
            qq1 = reshape(x(6*nknot + 6*nknot_bias + (1:3)), 3, 1);
            T1 = Ry(qq1(2)) * Rz(qq1(3)) * Rx(qq1(1));
            r1 = reshape(x(6*nknot + 6*nknot_bias + (4:6)), 3, 1);
            g131 = reshape(x(6*nknot + 6*nknot_bias + (7:8)), 2, 1);
            tshift1 = reshape(x(6*nknot + 6*nknot_bias + 9), 1, 1);
            ri1 = reshape(x(6*nknot + 6*nknot_bias + 9 + (1:nmarkers*3)), 3, []);
            
            figure, plot_LM(x0), sgtitle("before")
            figure, plot_LM(x), sgtitle("after")

            function [A, b, cost] = LM_iteration(x)

%                 tshifti = 0;
                tshifti = x(nknot*6 + nknot_bias*6 + 9);
                x(nknot*6 + nknot_bias*6 + 9) = 0;
                
                A_sparse = zeros(length(x) + ...
                    (length(index)*2 + length(index_bias)*2 + length(index_imuparam))^2 * (nknot + nknot_bias) + ...
                    (length(index)*2 + length(index_omcparam_fix) + length(index_omcparam_var))^2 * nknot * nmarkers + ...
                    (length(index_bias)*2)^2 * (nknot_bias-bord_bias) * integ_npoints, 1);
                i_sparse = ones(length(x) + ...
                    (length(index)*2 + length(index_bias)*2 + length(index_imuparam))^2 * (nknot + nknot_bias) + ...
                    (length(index)*2 + length(index_omcparam_fix) + length(index_omcparam_var))^2 * nknot * nmarkers + ...
                    (length(index_bias)*2)^2 * (nknot_bias-bord_bias) * integ_npoints, 1);
                j_sparse = i_sparse;

                b = zeros(length(x), 1);
                cost = 0;
                buf = 1;

                % IMU ==============================================================
                buflen = (length(index)*2 + length(index_bias)*2 + length(index_imuparam))^2;
                buf = buf(1) - 1 + (1 : buflen);
                jac = zeros(length(index)*2 + length(index_bias)*2 + length(index_imuparam));
                k0_last = 0; 
                k0_bias_last = 0;
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
                            k0_bias + index_bias + nknot*6 + nknot_bias*3
                            nknot*6 + nknot_bias*6 + index_imuparam];
                    vars = x(rows);

                    [wh, ah, wjac, ajac] = Calibration_OMC_IMU_cost_imu(...
                        vars.', mod(time2_norm(k), 1), dT, mod(time_norm_bias(k), 1), dT_bias);

                    jac = jac + wjac' * Nw_inv * wjac ...
                              + ajac' * Na_inv * ajac;
                    b(rows) = b(rows) - wjac' * Nw_inv * (w_imu(k,:) - wh')' ...
                                      - ajac' * Na_inv * (a_imu(k,:) - ah')';
                    cost = cost + (w_imu(k,:) - wh') * Nw_inv * (w_imu(k,:) - wh')' / 2 + ...
                                  (a_imu(k,:) - ah') * Na_inv * (a_imu(k,:) - ah')' / 2;

                    k0_last = k0;
                    k0_bias_last = k0_bias;
                end
                i_sparse(buf) = reshape(repmat(rows, [1, length(rows)]), [], 1);
                j_sparse(buf) = reshape(repmat(rows', [length(rows), 1]), [], 1);
                A_sparse(buf) = reshape(jac, [], 1);
                buf = buf + buflen;

               % OMC ===============================================================
                buflen = (length(index)*2 + length(index_omcparam_fix) + length(index_omcparam_var))^2;
                buf = buf(1) - 1 + (1 : buflen);
                jac = zeros(length(index)*2 + length(index_omcparam_fix) + length(index_omcparam_var));
                k0_last = 0; 
                rows = [];
                for i = 1 : nmarkers
                    for k = 1 : length(time_norm)

                        if any(isnan(mtrans(k,:,i)))
                            continue
                        end

                        k0 = floor(time_norm(k) + tshifti);
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
                                nknot*6 + nknot_bias*6 + index_omcparam_fix
                                nknot*6 + nknot_bias*6 + index_omcparam_var + 3*(i-1)];
                        vars = x(rows);

                        [mtransh, rjac] = Calibration_OMC_IMU_cost_omc(...
                            vars.', mod(time_norm(k) + tshifti, 1), dT);
                        jac = jac + rjac' * Nr_inv * rjac;
                        b(rows) = b(rows) - rjac' * Nr_inv * (mtrans(k,:,i) - mtransh')';
                        cost = cost + (mtrans(k,:,i) - mtransh') * Nr_inv * (mtrans(k,:,i) - mtransh')' / 2;

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
                        [wbiasd, abiasd, wbiasd_jac, abiasd_jac] = Calibration_OMC_IMU_cost_imubias( ...
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

                A = sparse(i_sparse, j_sparse, A_sparse, length(x), length(x));
            end
            
            function plot_LM(x)
%                 tshifti = 0;
                tshifti = x(nknot*6 + nknot_bias*6 + 9);
                x(nknot*6 + nknot_bias*6 + 9) = 0;

                % IMU ==============================================================
                wh = zeros(length(time2_norm), 3);
                ah = zeros(length(time2_norm), 3);
                for k = 1 : length(time2_norm)

                    k0 = floor(time2_norm(k));
                    k0_bias = floor(time_norm_bias(k));

                    rows = [k0 + index
                            k0 + index + nknot*3
                            k0_bias + index_bias + nknot*6
                            k0_bias + index_bias + nknot*6 + nknot_bias*3
                            nknot*6 + nknot_bias*6 + index_imuparam];
                    vars = x(rows);

                    [wh(k,:), ah(k,:)] = Calibration_OMC_IMU_cost_imu(...
                        vars.', mod(time2_norm(k), 1), dT, mod(time_norm_bias(k), 1), dT_bias);
                end
                cost_w = sum((w_imu - wh).^2 .* diag(Nw_inv)', 2);
                cost_a = sum((a_imu - ah).^2 .* diag(Na_inv)', 2);

               % OMC ===============================================================
                mtransh = NaN(length(time_norm), 3, nmarkers);
                for i = 1 : nmarkers
                    for k = 1 : length(time_norm)

                        if any(isnan(mtrans(k,:,i)))
                            continue
                        end

                        k0 = floor(time_norm(k) + tshifti);
                        if k0 < 0 || k0 >= nknot-bord
                            continue
                        end

                        rows = [k0 + index
                                k0 + index + nknot*3
                                nknot*6 + nknot_bias*6 + index_omcparam_fix
                                nknot*6 + nknot_bias*6 + index_omcparam_var + 3*(i-1)];
                        vars = x(rows);

                        mtransh(k,:,i) = Calibration_OMC_IMU_cost_omc(...
                            vars.', mod(time_norm(k) + tshifti, 1), dT);
                    end
                end
                cost_omc = squeeze(sum((mtrans - mtransh).^2 .* diag(Nr_inv)', 2));

                h1 = subplot(331); plot(mtime, cost_omc)
                ylabel('OMC cost'), grid on, xline(trange(1)), xline(trange(end))
                h4 = subplot(334); plot(itime, cost_w)
                ylabel('w cost'), grid on, xline(trange(1)), xline(trange(end))
                h8 = subplot(337); plot(itime, cost_a)
                ylabel('a cost'), grid on, xline(trange(1)), xline(trange(end))

                h2 = subplot(233); plot(itime, w_imu, '.-')
                hold on, set(gca, 'ColorOrderIndex', 1)
                plot(itime, wh)
                ylabel('\omega [rad/s]'), grid on, xline(trange(1)), xline(trange(end))
                h3 = subplot(236); plot(itime, a_imu, '.-')
                hold on, set(gca, 'ColorOrderIndex', 1)
                plot(itime, ah)
                ylabel('a [m/s^2]'), grid on, xline(trange(1)), xline(trange(end))

                h7 = subplot(332); plot(mtime, squeeze(mtrans(:,1,:)), '.-')
                hold on, set(gca, 'ColorOrderIndex', 1)
                plot(mtime, squeeze(mtransh(:,1,:)))
                ylabel('Body marker x'), grid on, xline(trange(1)), xline(trange(end))
                h5 = subplot(335); plot(mtime, squeeze(mtrans(:,2,:)), '.-')
                hold on, set(gca, 'ColorOrderIndex', 1)
                plot(mtime, squeeze(mtransh(:,2,:)))
                ylabel('Body marker y'), grid on, xline(trange(1)), xline(trange(end))
                h6 = subplot(338); plot(mtime, squeeze(mtrans(:,3,:)), '.-')
                hold on, set(gca, 'ColorOrderIndex', 1)
                plot(mtime, squeeze(mtransh(:,3,:)))
                ylabel('Body marker z'), grid on, xline(trange(1)), xline(trange(end))
                sgtitle("LM OMC IMU calib")

                linkaxes([h1, h2, h3, h4, h5, h6, h7, h8], 'x')
            end
        end
             
        function [mtime_norm, itime_norm, btime_norm, nknot, nknot_bias, ...
                integ_npoints, Nw, Na, Nr, Nwb, Nab] = get_bspline_vars(obj)
            
            [time, ~, ~, ~, time_imu] = obj.get_data();
            dT = obj.dT_;
            dT_bias = obj.dT_bias_;
            bord = obj.bord_;
            bord_bias = obj.bord_bias_;
            
            t0 = time(1);
            mtime_norm = (time - t0) / dT;
            itime_norm = (time_imu - t0) / dT;
            btime_norm = (time_imu - t0) / dT_bias;

            tend = time(end);
            nknot = ceil((tend - t0) / dT + 1e-10) + bord - 1;
            nknot_bias = ceil((tend - t0) / dT_bias + 1e-10) + bord_bias - 1;
            
            Fs_imu = 1 / median(diff(time_imu));
            integ_npoints = obj.integ_npoints_;
            Nw = (obj.w_white * sqrt(Fs_imu)).^2;
            Na = (obj.a_white * sqrt(Fs_imu)).^2;
            Nwb = obj.w_walk.^2;
            Nab = obj.a_walk.^2;
            Nr = obj.Nr_;
        end
        
        function [cq, wbias] = fit_bspline_orientation(obj, Tw)

            angle_std = 0.2e-3;
            vel_std = 0.1e-3;
 
            [kk, kk2, ~, nknot] = obj.get_bspline_vars();
            [mtime, quat, ~, ~, itime, wimu] = obj.get_data();
            bord = obj.bord_;
            dt = obj.dT_;
            
            w = wimu * inv(Tw)';
            quat_i = quatnormalize(interp1(kk, quat, kk2));
            qd = Calibration_OMC_IMU_pd(quat_i, w);
            q = unwrap(Calibration_OMC_IMU_p(quat));

            % position ==============================================
            A1 = obj.generate_predictor_matrix(bord, kk, nknot, 1);
            b1 = q;

            % velocity ================================================
            A2 = obj.generate_predictor_matrix(bord, kk2, nknot, 2) / dt;
            b2 = qd;

            % combine ==================================================
            A = [zeros(size(A1,1), 1), A1
                 ones(size(A2,1), 1), A2];
            b = [b1; b2];
            W = blkdiag(speye(length(kk )) / angle_std, ...
                        speye(length(kk2)) / vel_std);  % weight
   
            rows = any(A ~= 0, 1);
            if all(rows)
                cq_bias = (A' * W * A) \ (A' * W * b);
            else
                warning("Not all B coefficients can be estimated...");
                cq_bias = (A(:,rows)' * W * A(:,rows)) \ (A(:,rows)' * W * b);
                cq_bias = [cq_bias(1,:)
                           interp1(find(rows(2:end)), cq_bias(2:end,:), 1:length(rows)-1, 'linear', 'extrap')];
            end
            
            wbias = cq_bias(1,:) * Tw';
            cq = cq_bias(2:end,:);

            if nargout == 0
                bh = A * cq_bias;
                b1h = bh(1:size(b1,1),:);
                b2h = bh(size(b1,1)+1:end,:);

                fit = 1 - [goodnessOfFit(b1h, b1, 'NMSE'), goodnessOfFit(b2h, b2, 'NMSE')];
                fprintf('Orientation fitness (NMSE): [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f]\n', ...
                    fit(1), fit(2), fit(3), fit(4), fit(5), fit(6))

                figure
                h1 = subplot(221); plot(mtime, b1, '--'), title("angle, OMC")
                hold on, set(gca, 'ColorOrderIndex', 1), plot(mtime, b1h)
                h3 = subplot(223); plot(mtime, b1 - b1h), title('error')
                h2 = subplot(222); plot(itime, b2, '--'), title("velocity, IMU")
                hold on, set(gca, 'ColorOrderIndex', 1), plot(itime, b2h)
                h4 = subplot(224); plot(itime, b2 - b2h), title('error'), sgtitle('Fit B-spline orientation')
                linkaxes([h1, h2, h3, h4], 'x')
            end
        end

        function [cs, grav] = fit_bspline_translation(obj, Ta, r, abias)

            trans_std = 1;
            acc_std = 10;
 
            [kk, kk2, ~, nknot] = obj.get_bspline_vars();
            [mtime, quat, trans, ~, itime, ~, aimu, gaps] = obj.get_data();
            bord = obj.bord_;
            dt = obj.dT_;
            
            kk(gaps) = [];
            mtime(gaps) = [];
            quat(gaps,:) = [];
            trans(gaps,:) = [];
            
            % fit cs as the IMU position
            trans_imu = trans + quatrotate(quatinv(quat), r');
            quat_interp = quatnormalize(interp1(kk, quat, kk2, 'pchip', 'extrap'));
            as = quatrotate(quatinv(quat_interp), (aimu - abias') * inv(Ta)');

            % position ==============================================
            A1 = obj.generate_predictor_matrix(bord, kk, nknot, 1);
            b1 = trans_imu;

            % acceleration ================================================
            A2 = obj.generate_predictor_matrix(bord, kk2, nknot, 3) / dt^2;
            b2 = as;

            % combine ==================================================
            A = [zeros(size(A1,1), 1), A1 
                 ones(size(A2,1), 1), A2];
            b = [b1; b2];
            W = blkdiag(speye(length(kk)) / trans_std, ...
                        speye(length(kk2)) / acc_std);  % weight
            
            rows = any(A ~= 0, 1);
            if all(rows)
                cs_grav = (A' * W * A) \ (A' * W * b);
            else
                warning("Not all B coefficients can be estimated...");
                cs_grav = (A(:,rows)' * W * A(:,rows)) \ (A(:,rows)' * W * b);
                cs_grav = [cs_grav(1,:)
                           interp1(find(rows(2:end)), cs_grav(2:end,:), 1:length(rows)-1, 'linear', 'extrap')];
            end
            
            grav = cs_grav(1,:);
            cs = cs_grav(2:end,:);

            if nargout == 0
                bh = A * cs_grav;
                b1h = bh(1:size(b1,1),:);
                b2h = bh(size(b1,1)+1:end,:);

                fit = 1 - [goodnessOfFit(b1h, b1, 'NMSE'), goodnessOfFit(b2h, b2, 'NMSE')];
                fprintf('Translation fitness (NMSE): [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f]\n', ...
                    fit(1), fit(2), fit(3), fit(4), fit(5), fit(6))

                figure
                h1 = subplot(221); plot(mtime, b1, '--'), title("trans, OMC")
                hold on, set(gca, 'ColorOrderIndex', 1), plot(mtime, b1h)
                h3 = subplot(223); plot(mtime, b1 - b1h), title('error')
                h2 = subplot(222); plot(itime, b2, '--'), title("acc, IMU")
                hold on, set(gca, 'ColorOrderIndex', 1), plot(itime, b2h)
                h4 = subplot(224); plot(itime, b2 - b2h), title('error'), sgtitle('Fit B-spline translation')
                linkaxes([h1, h2, h3, h4], 'x')
                
                
            end
        end

        function [Tw, Ta, wbias, abias, r, g] = calibrate_lsq(obj)
            fprintf("LSQ Calibration\n")
            
            % tuning params ============================
            sgolay_params = obj.sgolay_params_;
            
            % =========================================
            
            [time, quat, trans, ~, time_imu, w_imu, a_imu, gaps] = obj.get_data();
            Fs = 1 / median(diff(time));
            gaps = imdilate(gaps, strel('rectangle', [sgolay_params(2)*2, 1]));

            [wd, w] = angular_rates(quat, Fs, sgolay_params);
            [~, ~, a] = deriv_sgolay(trans, Fs, sgolay_params);
            
            % time delay ================================
            toffset = 0;
            
            time_imu_shift = time_imu - toffset;

            % Gyro =============================================
            w_imu_interp = interp1(time_imu_shift, w_imu, time);

            rows = time > time(1) + 0.1 & time < time(end) - 0.1 & ~gaps;  % remove transient errors from derivative
            lm1 = fitlm(w(rows,:), w_imu_interp(rows,1), 'RobustOpts', 'on');
            lm2 = fitlm(w(rows,:), w_imu_interp(rows,2), 'RobustOpts', 'on');
            lm3 = fitlm(w(rows,:), w_imu_interp(rows,3), 'RobustOpts', 'on');

            Tw = [lm1.Coefficients.Estimate(2:end)'
                   lm2.Coefficients.Estimate(2:end)'
                   lm3.Coefficients.Estimate(2:end)'];
            wbias = [lm1.Coefficients.Estimate(1), lm2.Coefficients.Estimate(1), lm3.Coefficients.Estimate(1)]';

            % Acc ================================================
            skew_func = @(p) [0, -p(3), p(2)
                              p(3), 0, -p(1)
                              -p(2), p(1), 0];

            a_imu_interp = interp1(time_imu_shift, a_imu, time);
            R = quat2rotm(quat);

            dlen = size(quat, 1);
            A = zeros(dlen*3, 18);
            b = zeros(dlen*3, 1);
            valid = true(size(b));
            for i = 10 : dlen-10  % avoid transients
                rows = (i-1)*3 + (1:3);
                A(rows,:) = [blkdiag(a_imu_interp(i,:), a_imu_interp(i,:), a_imu_interp(i,:)), ...
                             skew_func(wd(i,:)) + skew_func(w(i,:))*skew_func(w(i,:)), ...
                             R(:,:,i)', eye(3)];
                b(rows) = a(i,:) * R(:,:,i);
                valid(rows) = ~gaps(i);
            end

            lm = fitlm(A(valid,:), b(valid), 'RobustOpts', 'on', 'intercept', false);

            Ta = inv(reshape(lm.Coefficients.Estimate(1:9), [3, 3])');
            r = -lm.Coefficients.Estimate(10:12);
            g = lm.Coefficients.Estimate(13:15);
            abias = -Ta * lm.Coefficients.Estimate(16:18);
            
            % Plot results ========================================
            % reconstruct IMU measurements using OMC + calib
            wh = w * Tw' + wbias';
            ah = quatrotate(quat, (a - g')) * Ta' + abias';
            wh(gaps,:) = NaN;
            ah(gaps,:) = NaN;
            w_error = wh - interp1(time_imu_shift, w_imu, time);
            a_error = ah - interp1(time_imu_shift, a_imu, time);
            
            fprintf("w error: MAE (%.2f, %.2f, %.2f) rad/s\n", nanmean(abs(w_error)))
            fprintf("a error: MAE (%.2f, %.2f, %.2f) m/s^2\n", nanmean(abs(a_error)))
            
            figure
            h1 = subplot(221); plot(time_imu_shift, w_imu, 'linewidth', 2), grid on, title('\omega [rad/s]')
            hold on, set(gca, 'ColorOrderIndex', 1), plot(time, wh, '-')
            h2 = subplot(222); plot(time_imu_shift, a_imu, 'linewidth', 2), grid on, title('a [m/s^2]')
            hold on, set(gca, 'ColorOrderIndex', 1), plot(time, ah, '-')
            h3 = subplot(223); plot(time, w_error)
            grid on, title('error [rad/s]')
            h4 = subplot(224); plot(time, a_error)
            grid on, title('error [m/s^2]'), sgtitle('LSQ OMC-IMU Calibration')
            linkaxes([h1, h2, h3, h4], 'x')
        end
        
        function [mtime, quat, trans, mtrans, itime, w, a] = get_raw_data(obj)
            mtime = obj.mtime_;
            quat = obj.quat_;
            trans = obj.trans_;
            mtrans = obj.mtrans_;
            itime = obj.itime_;
            w = obj.w_;
            a = obj.a_;
        end
        
        function [mtime, quat, trans, mtrans, itime, w, a, gaps] = get_data(obj)
            [mtime, quat, trans, mtrans, itime, w, a] = obj.get_raw_data();
            trange = obj.trange_;
            
            rows = mtime >= trange(1) & mtime <= trange(end);
            gaps = any(isnan(trans(rows,:)), 2);
            mtime = mtime(rows);
            R = fillmissing(reshape(quat2rotm(quat(rows,:)), 9, [])', 'linear', 'EndValues', 'nearest');
            quat = rotm2quat(reshape(R', 3, 3, []));
            for i = 2 : size(quat, 1)
                if sum(abs(quat(i,:) - quat(i-1,:))) > sum(abs(quat(i,:) + quat(i-1,:)))
                    quat(i,:) = -quat(i,:);
                end
            end
%             quat = quatnormalize(fillmissing(quat(rows,:), 'linear', 'EndValues', 'nearest'));
            trans = fillmissing(trans(rows,:), 'linear', 'EndValues', 'nearest');
            mtrans = mtrans(rows,:,:);
            
            rows = itime >= mtime(1) & itime <= mtime(end);
            itime = itime(rows);
            w = w(rows,:);
            a = a(rows,:);
        end
        
        function set_data(obj, mtime, quat, trans, mtrans, itime, w, a)
           obj.mtime_ = mtime;
           obj.quat_ = quat;
           obj.trans_ = trans;
           obj.mtrans_ = mtrans;
           obj.itime_ = itime;
           obj.w_ = w;
           obj.a_ = a;
        end
        

    end
    
    methods (Static)
        function generate_functions()
            % Assume diagonal cov noise
            poly_order = 4;
            poly_order_bias = 2;

            % orientation representation ============================
            [rotm_fcn, w_fcn] = Calibration_OMC_IMU.generate_EulerYZX_transformation();

            % B-spline basis ============================================
            syms t dt
            % body pose, translated to the IMU center
            cq = sym('cq', [poly_order+1, 3]);
            cs = sym('cs', [poly_order+1, 3]);
            cwbias = sym('cwbias', [poly_order_bias+1, 3]);
            cabias = sym('cabias', [poly_order_bias+1, 3]);

            M = sym(Calibration_OMC_IMU.build_bspline_matrix(poly_order + 1));
            M_bias = sym(Calibration_OMC_IMU.build_bspline_matrix(poly_order_bias + 1));

            % states ===========================================
            % sensor pose
            qp = ((t.^(0:poly_order)) * M * cq).';
            qdp = diff(qp, t) / dt;
            wp = w_fcn(qp, qdp);

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
            qq = sym('qq', [3, 1]);
            T = Ry(qq(2)) * Rz(qq(3)) * Rx(qq(1));
%             quat = [sqrt(1 - qq.'*qq), qq.'];
%             T = quat2rotm_sym(quat);
%             Ta = sym('Ta', [3, 3]);
%             Tw = sym('Tw', [3, 3]);
            g13 = sym('g', [2, 1]);
            g = [g13(1); -sqrt(9.80136^2 + g13.'*g13); g13(2)];
            ri = sym('ri', [3, 1]);  % plate to marker, plate frame
            syms tshift

            % Model ==============================================
            % IMU
            R = rotm_fcn(qp);  % plate orientation

            w_imu = T * wp + wbias;
            a_imu = T * R.' * (as - g) + abias;

            sbody = ss - R * r;   % plate
            rj = sbody + R * ri;  % marker
            rj_shift = subs(rj, t, t+tshift);

            % % calibration
            vars_imu = [cq(:); cs(:); cwbias(:); cabias(:); qq; r; g13];
            vars_bias = [cwbias(:); cabias(:)];
            vars_omc = [cq(:); cs(:); r; tshift; ri];

            w_imu_jac = jacobian(w_imu, vars_imu);
            a_imu_jac = jacobian(a_imu, vars_imu);
            rj_jac = jacobian(rj_shift, vars_omc);
            wbiasd_jac = jacobian(wbiasd, vars_bias);
            abiasd_jac = jacobian(abiasd, vars_bias);

            matlabFunction(w_imu, a_imu, w_imu_jac, a_imu_jac, 'Vars', ...
                {vars_imu.', t, dt, t_bias, dt_bias}, ...
                'File', 'Calibration_OMC_IMU_cost_imu', 'Optimize', true)
            matlabFunction(wbiasd, abiasd, wbiasd_jac, abiasd_jac, 'Vars', ...
                {vars_bias.', t_bias, dt_bias}, ...
                'File', 'Calibration_OMC_IMU_cost_imubias', 'Optimize', true)
            matlabFunction(rj_shift, rj_jac, 'Vars', ...
                {vars_omc.', t, dt}, ...
                'File', 'Calibration_OMC_IMU_cost_omc', 'Optimize', true)
        end
        
        function [rotm_fcn, w_fcn] = generate_EulerYZX_transformation()
            % rotm_fcn: R <- (p)
            % w_fcn: w <- (p, pd)
            % pd_fcn: pd <- (quat, w)
            
            p = sym('p', [3, 1]);
            R = Ry(p(2)) * Rz(p(3)) * Rx(p(1));
            rotm_fcn = matlabFunction(R, 'Vars', {p});
            
            % =============================
            w = sym('w', [3, 1]);
            quat = sym('quat', [1, 4]);
            syms t
            
            assume(t, 'real'); assume(w, 'real'); assume(quat, 'real'); 
            quatt = str2sym('[qt0(t), qt1(t), qt2(t), qt3(t)]');
            assume(quatt, 'real'); 
            
            quatd = 0.5 * quatmultiply_sym(quat, [0, w.']);
            Rquat = quat2rotm_sym(quatt);
            p = [atan2(-Rquat(2,3), Rquat(2,2))
                 atan2(-Rquat(3,1), Rquat(1,1))
                 asin(Rquat(2,1))];
%             pd = diff([-atan(Rquat(2,3) / Rquat(2,2))
%                         -atan(Rquat(3,1) / Rquat(1,1))
%                         asin(Rquat(2,1))], t);  % YZX
            pd = diff(p, t);        
            p = subs(p, quatt, quat);
            pd = subs(subs(pd, diff(quatt), quatd), quatt, quat);
            
            matlabFunction(pd.', 'Vars', {quat, w.'}, 'File', 'Calibration_OMC_IMU_pd');
            matlabFunction(p.', 'Vars', {quat}, 'File', 'Calibration_OMC_IMU_p');
            
            % ==============================
            p = sym('p', [3, 1]);
            pd = sym('pd', [3, 1]);

%             R1 = eye(3);
%             R2 = Rx(p(1)).' * Rz(p(3)).';
%             R3 = Ry(p(2)).';
%             S = [R1(1,:); R2(2,:); R3(3,:)];
%             w = S * pd;
            
            wG = [0; pd(2); 0] + Ry(p(2)) * [0; 0; pd(3)] + Ry(p(2)) * Rz(p(3)) * [pd(1); 0; 0];
            R = Ry(p(2)) * Rz(p(3)) * Rx(p(1));
            w = R.' * wG;
            
            w_fcn = matlabFunction(w, 'Vars', {p, pd});
        end
        
        function M = build_bspline_matrix(bspline_order)
            M = 1;
            for k = 1 : bspline_order - 1
                M = ([M; zeros(1,k)] * ([diag(1:(k)), zeros(k,1)] + [zeros(k,1), diag((k-1):-1:0)]) + ...
                     [zeros(1,k); M] * ([-eye(k), zeros(k,1)] + [zeros(k,1), eye(k)]) ) / k;
            end
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
        
        function time2_fix = sync_quat_w(time, quat, time2, w)

            Fs_omc = 1 / median(diff(time));
            Fs_sd = 1 / median(diff(time2));

            % First pass, use LP signals ===========================================
            [bft, aft] = butter(5, 1/(Fs_omc/2));  % LP

            [~, wa] = angular_rates(quat, Fs_omc, [5, 11]);
            wa_abs = fillmissing(sqrt(sum(wa.^2, 2)), 'constant', 0);
            wa_filt = filtfilt(bft, aft, wa_abs);

            % downsample IMUs
            [bft, aft] = butter(5, 1/(Fs_sd/2));
            wb_abs = sqrt(sum(w.^2, 2));
            wb_filt_down = interp1(time2, filtfilt(bft, aft, wb_abs), time);

            d = finddelay(wa_filt, wb_filt_down, round(60*Fs_omc));
            time2_fix = time2 - d/Fs_omc;

            % Second pass use unfiltered signal =======================================

            dT = 20.0;  % comparison window

            if nargout == 0, figure, end
            for reps = 1 : 6
                wa_abs_up = interp1(time, wa_abs, time2_fix);

                tcomp = ((time(1)+dT) : dT : (time(end)-dT))';
                dd = zeros(size(tcomp));
                for i = 1 : length(tcomp)
                    rows = time2_fix > tcomp(i) - dT/2 & time2_fix < tcomp(i) + dT/2;

                    dd(i) = finddelay(wa_abs_up(rows), wb_abs(rows), round(1.0*Fs_sd));
                end

                lm = fitlm(tcomp, dd, 'RobustOpts', 'on');

                time2_fix = (time2_fix - lm.Coefficients.Estimate(1)/Fs_sd) * ...
                    (-lm.Coefficients.Estimate(2)/Fs_sd + 1);

                if nargout == 0
                    subplot(2,3,reps), plot(lm)
                end
            end

            assert(abs(lm.Coefficients.Estimate(1)/Fs_sd) < 1e-3 && ...
                   lm.Coefficients.Estimate(2)/Fs_sd < 1e-4 && ...
                   prctile(abs(wa_abs_up - wb_abs), 90) < 2, "Time sync failed")

            % nanstd(wa_abs_up - wb_abs)
            % prctile(abs(wa_abs_up - wb_abs), 90)

            if nargout == 0
                figure, 
                subplot(121),
                plot(time2_fix, wb_abs, '.-')
                hold on, plot(time, wa_abs, '.-')
                subplot(122), histogram(wa_abs_up - wb_abs, 100)
            end
        end
    end
        
end

