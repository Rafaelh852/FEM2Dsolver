function [L2norm,H1norm]= H1norm2d(f,gradf,g,nodes,triangles)
noOfIntegPt = 7;

func = @(x,y)( f(x,y)^2 + gradf(x,y)'*gradf(x,y));

ff = @(x,y)(g(x,y)^2);

L2norm = 0;
z = 0;

for e = 1:length(triangles)
    
    x1 = nodes(triangles(e,1),1);
    x2 = nodes(triangles(e,2),1);
    x3 = nodes(triangles(e,3),1);

    y1 = nodes(triangles(e,1),2);
    y2 = nodes(triangles(e,2),2);
    y3 = nodes(triangles(e,3),2);

    
   L2norm = L2norm + quad2dTri(ff,x1,y1,x2,y2,x3,y3,noOfIntegPt);
    
    z = z + quad2dTri(func,x1,y1,x2,y2,x3,y3,noOfIntegPt);
end

  L2norm = sqrt(L2norm);
  H1norm = sqrt(sum(z));
  
  
  
end
