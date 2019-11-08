function [ubf] = ubfilter(ub)

%   This function filters ubase with the robot model
v0p = 0;
w0p = 0;
vp = 0;
wp = 0;

[m,n] = size(ub);

for l = 1:n
    
    v(l) = 0.2881*v0p + 0.7105*vp;
    vp = v(l);
    v0p = ub(1,n);

    w(l) = 0.2232*w0p + 0.755*wp;
    wp = w(l);
    w0p = ub(2,n);
    
end

ubf = [v;w];

end
