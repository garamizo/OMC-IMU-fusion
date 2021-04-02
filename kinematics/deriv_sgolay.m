function [yd, y, ydd] = deriv_sgolay(x, Fs, order)
% DERIV_SGOLAY  Derivative using Sarvitz Golay filter
%
%   [yd, yf, ydd] = DERIV_SGOLAY(y, Fs) reads derivative, filtered input,
%   and 2nd derivative

% if nargin == 1
%     Fs = 350;
% end

dt = 1/Fs;
xmean = nanmean(x);
% x = x - xmean;
x = bsxfun(@minus, x, xmean);

if nargin == 3
    [~,g] = sgolay(order(1), order(2));
else

% fprintf('hey!\n')
[~,g] = sgolay(11,15);  % default
% [~,g] = sgolay(15, 31);  % lower variance
% [~,g] = sgolay(5, 21);

% [~,g] = sgolay(3, 15);  % ankle angle
% [~,g] = sgolay(3, 7);  % ptorque, pforce

% [~,g] = sgolay(3, 5);  % COP
end

dx = zeros(size(x,1),size(x,2),3);
for c = 1 : size(x, 2)
    for p = 0:2
      dx(:,c,p+1) = conv(x(:,c), factorial(p)/(-dt)^p * g(:,p+1), 'same');
    end
end

% y = dx(:,:,1) + xmean;
y = bsxfun(@plus, dx(:,:,1), xmean);
yd = dx(:,:,2);
ydd = dx(:,:,3);

ord = order(1);
y(1:ord,:) = repmat(y(ord,:), [ord, 1]);
yd(1:ord,:) = repmat(yd(ord,:), [ord, 1]);
ydd(1:ord,:) = repmat(ydd(ord,:), [ord, 1]);
y(end-ord+1:end,:) = repmat(y(end-ord+1,:), [ord, 1]);
yd(end-ord+1:end,:) = repmat(yd(end-ord+1,:), [ord, 1]);
ydd(end-ord+1:end,:) = repmat(ydd(end-ord+1,:), [ord, 1]);



