function [val] = quad2dTri(f, x1, y1, x2, y2, x3, y3, noOfIntegPt)

    x = @(xi,eta)(x1 + (x2-x1)*xi + (x3-x1)*eta);
    y = @(xi,eta)(y1 + (y2-y1)*xi + (y3-y1)*eta);
    
    g = @(xi,eta)(f(x(xi,eta),y(xi,eta)));
    
    detA = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
    
    val = abs(detA)*quad2dTs(g, noOfIntegPt);
    
end