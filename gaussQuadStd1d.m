function [y] = gaussQuadStd1d(g,noOfIntegPt)
    
    if noOfIntegPt ==2 
        
        y = g(-sqrt(3)/3) + g(sqrt(3)/3);
        
       return 
    end

    if noOfIntegPt == 3
        y = (5/9)*g(-sqrt(3/5)) + (8/9)*g(0) + (5/9)*g(sqrt(3/5));
        
        return
    end
    
    if noOfIntegPt == 5
       y = (322-13*sqrt(70))*(1/900)*(g((1/3)*sqrt(5+2*sqrt(10/7)))+g((-1/3)*sqrt(5+2*sqrt(10/7))))...
           + (128/225)*g(0)...
           + (322+13*sqrt(70))*(1/900)*(g((1/3)*sqrt(5-2*sqrt(10/7)))+ g((-1/3)*sqrt(5-2*sqrt(10/7)))) ;
        
        return
    end
    
    "noOfIntegPt is not 2 or 3"
    y = 0;
    return
end
