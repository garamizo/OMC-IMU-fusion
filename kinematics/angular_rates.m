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
    
    
    R = quat2rotm(quat);
    [rd, ~, rdd] = deriv_sgolay(reshape(R, 9, [])', varargin{:});
    Rd = reshape(rd', 3, 3, []);
    Rdd = reshape(rdd', 3, 3, []);
    for i = 1 : size(quat, 1)
        wskew = R(:,:,i)' * Rd(:,:,i);
        w(i,:) = [wskew(3,2) - wskew(2,3), wskew(1,3) - wskew(3,1), wskew(2,1) - wskew(1,2)] / 2;
        
        % this wd is incorrect. Use previous method
%         wskew = Rd(:,:,i)' * Rd(:,:,i) + R(:,:,i)' * Rdd(:,:,i);
%         wd(i,:) = [wskew(3,2) - wskew(2,3), wskew(1,3) - wskew(3,1), wskew(2,1) - wskew(1,2)] / 2;
    end
end