function [z] = fei2d(k,f,g,h,e,i,triangles,nodes,edges,triangleMidPts,midNodes,Gamma1Nodes,Gamma1MidNodes,Gamma2Edges,shapeFn,noOfIntegPt,noOfIntegPt1d)
    
    n = length(nodes);
    nt = length(triangles);
   
    
    x1 = nodes(triangles(e,1),1);
    x2 = nodes(triangles(e,2),1);
    x3 = nodes(triangles(e,3),1);

    y1 = nodes(triangles(e,1),2);
    y2 = nodes(triangles(e,2),2);
    y3 = nodes(triangles(e,3),2);

    %first term code 
    f1 = @(x,y)(f(x,y)*shapeFn2d(i, x, y, x1, y1, x2, y2, x3, y3, shapeFn));
    z1 = quad2dTri(f1, x1, y1, x2, y2, x3, y3, noOfIntegPt);


    %second term code    
    if shapeFn == 1

        j = find(Gamma1Nodes>0);
        z2 = 0;
        p1 = (triangles(e,1));
        p2 = (triangles(e,2));
        p3 = (triangles(e,3));
        s = 1;
        for p = [p1, p2, p3]
            if any(p == j)
               z2 = z2 + g(nodes(p,1),nodes(p,2))*keij2d(k,e,i,s,triangles,nodes,shapeFn,noOfIntegPt);
            end
            s = s+1;
        end 


    end
    if shapeFn == 2
        triangles = [triangles;triangleMidPts];
        gnodes = [Gamma1Nodes;Gamma1MidNodes];
        nodes = [nodes;midNodes];

        j = find(gnodes>0);

        z2 = 0;

        p1 = (triangles(e,1));
        p2 = (triangles(e,2));
        p3 = (triangles(e,3));

        mp1 = (triangles(e+nt,1));
        mp2 = (triangles(e+nt,2));
        mp3 = (triangles(e+nt,3));

        s = 1;
        for p = [p1, p2, p3, mp1, mp2, mp3]
            if any(p == j)  
               z2 = z2 + g(nodes(p,1),nodes(p,2))*keij2d(k,e,s,i,triangles,nodes,shapeFn,noOfIntegPt);
            end
            s = s+1;
        end

    end

 
    %third term code
    f3  = @(x,y)(h(x,y)*shapeFn2d(i, x, y, x1, y1, x2, y2, x3, y3, shapeFn));
    z3 = 0;
    
    edg1 = [p1 p2]' ;
    edg2 = [p2 p3]' ;    
    edg3 = [p1 p3]' ;
 
     edg11 = [p2 p1]' ;
     edg22 = [p3 p2]' ;
     edg33 = [p3 p1]' ;
    
    for edg = [edg1 edg2 edg3 edg11 edg22 edg33]
        
         ed = edg' == edges(Gamma2Edges>0,1:2);
      if norm(ed(:,1) + ed(:,2),inf) == 2
         x1 = nodes(edg(1),1);
         y1 = nodes(edg(1),2);

         x2 = nodes(edg(2),1);
         y2 = nodes(edg(2),2);

         fn = @(s)(f3((1-s)*x1 + s*x2,(1-s)*y1+s*y2));
         z3 = z3 + sqrt((x2-x1)^2 +(y2-y1)^2)*gaussQuad1d(fn,0,1, noOfIntegPt1d);
      end
    end

        
    z = z1 - z2 + z3;
end
