function [F] = loadF2d(k,f,g,h,triangles,nodes,edges,triangleMidPts,midNodes,indVec,Gamma1Nodes,Gamma1MidNodes,Gamma2Edges,shapeFn,noOfIntegPt,noOfIntegPt1d)

freenodes = find(indVec);

F = zeros(length(freenodes), 1);

totaltriangles = [triangles , triangleMidPts];

for e = 1:length(triangles)
    for i = 1:3*shapeFn
         if ( ismember(totaltriangles(e,i),freenodes))
            F(indVec(totaltriangles(e,i))) =  F(indVec(totaltriangles(e,i))) + fei2d(k,f,g,h,e,i,triangles,nodes,edges,triangleMidPts,midNodes,Gamma1Nodes,Gamma1MidNodes,Gamma2Edges,shapeFn,noOfIntegPt,noOfIntegPt1d);
         end
    end
end

end