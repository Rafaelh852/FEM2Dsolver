function [z] = keij2d(k,e,i,j,triangles,nodes,shapeFn,noOfIntegPt)

x1 = nodes(triangles(e,1),1);
x2 = nodes(triangles(e,2),1);
x3 = nodes(triangles(e,3),1);

y1 = nodes(triangles(e,1),2);
y2 = nodes(triangles(e,2),2);
y3 = nodes(triangles(e,3),2);



f = @(x,y)(k(x,y)*shapeFnGrad2d(i, x, y, x1, y1, x2, y2, x3, y3, shapeFn)'*shapeFnGrad2d(j, x, y, x1, y1, x2, y2, x3, y3, shapeFn));

z = quad2dTri(f, x1, y1, x2, y2, x3, y3, noOfIntegPt);

end