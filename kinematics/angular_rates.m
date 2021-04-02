function [wd, w] = angular_rates(quat, varargin)
    % ANGULAR_RATES  read 2nd and 1st derivatives
    %
    %   [wd, w] = Mocap.ANGULAR_RATES(quat, Fs)
    % angular velocity on body frame
    quat = quatnormalize(quat);
    [quatd, ~, quatdd] = deriv_sgolay(quat, varargin{:});
    
    w = 2 * quatmultiply(quatinv(quat), quatd);  % body frame
    w(:,1) = 0;
    wd = quatmultiply(quatinv(quat), 2*quatdd - quatmultiply(quatd, w));
    w = w(:,2:4);
    wd = wd(:,2:4);

end