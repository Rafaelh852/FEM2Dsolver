function [val] = quad2dTs(g, noOfIntegPt)
    w = 0;
    
    if noOfIntegPt == 4
        w = [-9/32; 25/96; 25/96; 25/96];
        xi = [1/3;3/5;1/5;1/5];
        eta = [1/3;1/5;1/5;3/5];
    end
    
    if noOfIntegPt == 6
       w(1:3) = 137/2492;
       w(4:6) = 1049/9392;
       xi(1) = 1280/1567; 
       eta(2) = xi(1);
       xi(2:3) = 287/3134;
       eta([1;3]) = xi(2);
       xi(4) = 575/5319;
       eta(5) = xi(4);
       xi(5:6) = 2372/5319; 
       eta([4;6]) = xi(5);
       
    end 
   
    if noOfIntegPt == 7
       w(1) = 9/80;
       w(2:4) = 352/5590;
       w(5:7) = 1748/26406;
       xi(1) = 1/3;
       eta(1) = 1/3;
       xi(2) = 248/311;
       eta(3) = xi(2);
       xi(3:4) =496/4897;
       eta([2;4]) = xi(3);
       xi(5) = 248/4153;
       eta(6) = xi(5);
       xi(6:7) = 496/1055;
       eta([5;7]) = xi(6);
       
    end
    
    val = 0;
    for i = 1:length(w)
        val = val + w(i)*g(xi(i),eta(i));  
    end
   
end