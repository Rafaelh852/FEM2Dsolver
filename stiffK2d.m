function [K] = stiffK2d(k, triangles, nodes, triangleMidPts, indVec, shapeFn, noOfIntegPt)

%indvec gives the free node # for node j, returns 0 if its a free node. .
freenodes = find(indVec); % find() removes the zero elements of a vector

% Initialize K 
K = zeros(length(freenodes), length(freenodes));
% Append the midpoint triangles. 
totaltriangles = [triangles , triangleMidPts];

% Calculate Keij for each triangle, then put it in global K
for e = 1 : size(triangles,1) %for each triangle
    for i = 1 : 3*shapeFn  %for each row 
        for j = 1 : 3*shapeFn %for each column
            %Check to see if local node i and j are free nodes. 
            if ( ismember(totaltriangles(e,i),freenodes) && ismember(totaltriangles(e,j), freenodes) )                   
                K(indVec(totaltriangles(e,i)), indVec(totaltriangles(e,j)))  = K(indVec(totaltriangles(e,i)), indVec(totaltriangles(e,j))) + keij2d(k,e,i,j,triangles,nodes,shapeFn,noOfIntegPt);                    
            end
        end
    end
end

K = sparse(K);
end