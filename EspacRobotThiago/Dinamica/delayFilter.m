function [ubf] = delayFilter(ub,dv,dw)

%   This function filters ubase with the robot model

% dv and dw are the delay for v and w, respectively

wp = zeros(1,dw);
vp = zeros(1,dv);

[m,n] = size(ub);

for l = 1:n
    
    % for v
    
    if dv > 0
        
        v(l) = vp(dv);
        
        if dv > 1
            
            for j = dv:2
                
                vp(j) = vp(j-1);
                
            end
            
        end
        
        vp(1) = ub(1,l);
        
    else
        
        v(l) = ub(1,l);
        
    end
    
    % for w
    
    if dw > 0
        
        w(l) = wp(dw);
        
        if dw > 1
            
            for j = dw:2
                
                wp(j) = wp(j-1);
                
            end
            
        end
        
        wp(1) = ub(2,l);
        
    else
        
        w(l) = ub(2,l);
        
    end
    
end

ubf = [v;w];

end
