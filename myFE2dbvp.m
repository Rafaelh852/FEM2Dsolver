function [uh, graduh, hMax, nodes, triangles] = myFE2dbvp(k, f, g, h, maxh, nodes, triangles, edges, bdyNode, bdyEdge, curveEdge, bdyFn, bdyFnder, shapeFn)
    
    noOfIntegPt = 7;
    noOfIntegPt1d = 5;
    
     while(  findhmax(edges,nodes) > maxh)
        [nodes,triangles,edges,bdyNode,bdyEdge,curveEdge]=refine(nodes,triangles,edges,bdyNode,bdyEdge,curveEdge,bdyFn,bdyFnder);
     end
     
    hMax = findhmax(edges,nodes);

   [midNodes, triangleMidPts] = midNT(nodes,triangles);
   
    if shapeFn==2 
        
        Gamma1Nodes =bdyChk(nodes,bdyFn);
        
        Gamma1MidNodes = zeros(length(midNodes),1);
        
        for e = edges(curveEdge>0,1:2)'
            e1 = e(1)';
            e2 = e(2)';
            em = 0.5*(nodes(e1,1:2) + nodes(e2,1:2));
            Gamma1MidNodes = Gamma1MidNodes +(sum(em == midNodes,2)==2);
        end

    end
    if shapeFn == 1
        Gamma1Nodes =bdyChk(nodes,bdyFn);
         Gamma1MidNodes = bdyChk(midNodes,bdyFn);
    end 
    

    Gamma2Edges =  bdyEdge - curveEdge;
     
    
    indVec = indVect(Gamma1Nodes,Gamma1MidNodes, shapeFn);
    
    K = stiffK2d(k, triangles, nodes, triangleMidPts, indVec, shapeFn, noOfIntegPt);
    
    F = loadF2d(k,f,g,h,triangles,nodes,edges,triangleMidPts,midNodes,indVec,Gamma1Nodes,Gamma1MidNodes,Gamma2Edges,shapeFn,noOfIntegPt,noOfIntegPt1d);
    
    
    [Kbar,Q] = permRCM(K);
    
    Fbar = Q*F;

    [p, q] = bandwidth(Kbar);
    LU = bandLU(Kbar,p,q);
    L = tril(LU,-1) + diag(ones(1,length(LU)));
    U = triu(LU);
    
    y = bandforward(L,Fbar,p);
    Wbar = bandbackward(U,y,p);
    
    W = Q.'*Wbar;
    
    uh = approxSol2d(W,g,indVec,nodes,triangles,midNodes,triangleMidPts,shapeFn);
    graduh = approxSolGrad2d(W,g,indVec,nodes,triangles,midNodes,triangleMidPts,shapeFn);
   
end