function [G] = bdyChk(n,bdyFn)
    tolerance = 10^-14;
    G = zeros(length(n(:,1)),1);
    
    for i = 1:length(n(:,1))
        if ((n(i,2) - bdyFn(n(i,1))<= tolerance) && (abs(n(i,1)) < 1/2))
            G(i) = 1;
           
        end
        
    end
    
end