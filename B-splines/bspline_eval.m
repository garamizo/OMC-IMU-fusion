function [y, yd, ydd] = bspline_eval(c, kk, bord)
% kk = (t - t0) / dt

ndim = size(c, 2);  % signal dimensions

M = create_bspline_coefficients(bord);
idx = floor(kk);
tt = kk - idx;

ct = reshape(c((1 : bord) + idx,:), [], bord, ndim);
y = squeeze(sum((tt.^(0 : bord-1)) * M .* ct, 2));

if nargout > 1
    tmp = ((1 : bord-1) .* tt.^(0 : bord-2)) * M(2:end,:) .* ct;
    yd = squeeze(sum(tmp, 2));
end

if nargout > 2
    tmp = ((1 : bord-2) .* (2 : bord-1) .* tt.^(0 : bord-3)) * M(3:end,:) .* ct;
    ydd = squeeze(sum(tmp, 2));
end



