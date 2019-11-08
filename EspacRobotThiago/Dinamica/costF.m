function [Qe] = costF(err_x,err_y,err_theta)

%   Cost Function

Qe = sum(err_x.^2) + sum(err_y.^2) + sum(err_theta.^2);

end

