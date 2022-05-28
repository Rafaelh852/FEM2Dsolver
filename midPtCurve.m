function [xm, ym] = midPtCurve(f,fder,x1,x2)

    %arclength integrand
    h = @(x)(sqrt(1+fder(x)*fder(x)));
    
    %arclength function
    L =@(x)(gaussQuad1d(h,x1,x,5));
    
    % midpoint function to find zero of
    m = @(x)(L(x) - L(x2)/2);
    
    %finding the zero
    xm = fzero(m,x1);
    
    %returning the midpoint
    ym = f(xm);

end