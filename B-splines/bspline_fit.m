function c = bspline_fit(bord, dt, t, y, t2, yd, mode)
% LSQ fit

    if nargin == 4
        t0 = t(1);
        nknot = ceil((t(end) - t0) / dt + 1e-10) + bord - 1;
        Fs = 1/(t(2)-t(1));
        
        c = bspline_fit_zero(bord, (t - t0) / dt, y, nknot);
        
    elseif nargin == 6 || (nargin == 7 && strcmp(mode, 'first'))
        t0 = min(t(1), t2(1));
        nknot = ceil((max(t(end), t2(end)) - t0) / dt + 1e-10) + bord - 1;
        Fs = max(1/(t(2)-t(1)), 1/(t2(2)-t2(1)));
        
        c = bspline_fit_first_doffset(bord, (t - t0) / dt, y, (t2 - t0) / dt, yd * dt, nknot);
    elseif nargin == 7 && strcmp(mode, 'second')
        t0 = min(t(1), t2(1));
        nknot = ceil((max(t(end), t2(end)) - t0) / dt + 1e-10) + bord - 1;
        Fs = max(1/(t(2)-t(1)), 1/(t2(2)-t2(1)));
        
        c = bspline_fit_second(bord, (t - t0) / dt, y, (t2 - t0) / dt, yd * dt^2, nknot);
    else
        error('Wrong outputs')
    end
    
    tt = (t(1) : 0.1/Fs : t(end))';
    yh = bspline_eval(c, (tt - t0) / dt, bord);
    
    figure, plot(tt, yh, t, y, 'o')
    legend([compose("fit %d", 1:size(y,2)), compose("data %d", 1:size(y,2))])
    grid on
end

function [c, A, y] = bspline_fit_zero(bord, kk, y, nknot)
    % minimum number of knots needed
    nsamples = length(kk);
    
    N = ceil(bord / mean(diff(kk)));
    aa = zeros(N, nknot);
    jj = zeros(N, nknot);
    ii = zeros(N, nknot);
    
    ci = zeros(nknot, 1);
    ci(1) = 1;
    
    for j = 1 : nknot

        k0 = find(kk >= j, 1, 'first');
        if isempty(k0)  % border cases
            k0 = nsamples + 1;
        elseif k0 - N < 1
            k0 = N + 1;
        end
        rows = k0 + (-N : -1);
        ii(:,j) = rows;
        jj(:,j) = j;
        aa(:,j) = bspline_eval(ci, kk(rows), bord);
        
        ci = circshift(ci, 1);
    end
    
    valid = ii(:) > 0 & jj(:) > 0;
    A = sparse(ii(valid), jj(valid), aa(valid), nsamples, nknot);
    c = A \ y;
end

% function [c, A, y] = bspline_fit_zero(bord, kk, y, nknot)
%     % minimum number of knots needed
%     nsamples = length(kk);
% 
%     % create predictor matrix
%     A = zeros(nsamples, nknot);
%       
%     ci = zeros(nknot, 1);
%     ci(1) = 1;
%     
%     for i = 1 : nknot
% 
%         rows = kk > (i - bord) & kk < i;
%         A(rows,i) = bspline_eval(ci, kk(rows), bord);
%         
%         ci = circshift(ci, 1);
%     end
%     c = A \ y;
% end

function [c, A2, yd] = bspline_fit_first(bord, kk, y, kk2, yd, nknot)

    % position ==============================================
    [~, A, b] = bspline_fit_zero(bord, kk, y, nknot);

    % velocity ================================================
    nsamples = length(kk2);
    N = ceil(bord / mean(diff(kk2)));
    aa = zeros(N, nknot);
    jj = zeros(N, nknot);
    ii = zeros(N, nknot);
    
    ci = zeros(nknot, 1);
    ci(1) = 1;

    for j = 1 : nknot
        k0 = find(kk2 >= j, 1, 'first');
        if isempty(k0)  % border cases
            k0 = nsamples + 1;
        elseif k0 - N < 1
            k0 = N + 1;
        end
        rows = k0 + (-N : -1);
        ii(:,j) = rows;
        jj(:,j) = j;
        [~, aa(:,j)] = bspline_eval(ci, kk2(rows), bord);
        
        ci = circshift(ci, 1);
    end

    valid = ii(:) > 0 & jj(:) > 0;
    A2 = sparse(ii(valid), jj(valid), aa(valid), nsamples, nknot);
    % combine ==================================================

    c = [A; A2] \ [b; yd];
end

function [c, A2, yd] = bspline_fit_first_doffset(bord, kk, y, kk2, yd, nknot)

    % position ==============================================
    [~, A, b] = bspline_fit_zero(bord, kk, y, nknot);

    % velocity ================================================
    nsamples = length(kk2);
    N = ceil(bord / mean(diff(kk2)));
    aa = zeros(N, nknot);
    jj = zeros(N, nknot);
    ii = zeros(N, nknot);
    
    ci = zeros(nknot, 1);
    ci(1) = 1;

    for j = 1 : nknot
        k0 = find(kk2 >= j, 1, 'first');
        if isempty(k0)  % border cases
            k0 = nsamples + 1;
        elseif k0 - N < 1
            k0 = N + 1;
        end
        rows = k0 + (-N : -1);
        ii(:,j) = rows;
        jj(:,j) = j;
        [~, aa(:,j)] = bspline_eval(ci, kk2(rows), bord);
        
        ci = circshift(ci, 1);
    end

    
    valid = ii(:) > 0 & jj(:) > 0;
    A2 = sparse(ii(valid), jj(valid), aa(valid), nsamples, nknot);
    % combine ==================================================

%     c = [A; A2] \ [b; yd];
    cb = [A, zeros(size(A,1), 1)
         A2, ones(size(A2,1), 1)] \ [b; yd];
     
    c = cb(1:end-1,:);
    intercept = cb(end,:) / 5e-3
end

function [c, A2, ydd] = bspline_fit_second(bord, kk, y, kk2, ydd, nknot)

    % position ==============================================
    [~, A, b] = bspline_fit_zero(bord, kk, y, nknot);

    % acc ================================================
    nsamples = length(kk2);
    N = ceil(bord / mean(diff(kk2)));
    aa = zeros(N, nknot);
    jj = zeros(N, nknot);
    ii = zeros(N, nknot);
    
    ci = zeros(nknot, 1);
    ci(1) = 1;

    for j = 1 : nknot
        k0 = find(kk2 >= j, 1, 'first');
        if isempty(k0)  % border cases
            k0 = nsamples + 1;
        elseif k0 - N < 1
            k0 = N + 1;
        end
        rows = k0 + (-N : -1);
        ii(:,j) = rows;
        jj(:,j) = j;
        [~, ~, aa(:,j)] = bspline_eval(ci, kk2(rows), bord);
        
        ci = circshift(ci, 1);
    end

    valid = ii(:) > 0 & jj(:) > 0;
    A2 = sparse(ii(valid), jj(valid), aa(valid), nsamples, nknot);
    % combine ==================================================
    c = [A; A2] \ [b; ydd];
end



