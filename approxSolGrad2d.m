% this function takes the values of the nodes and gives the
% approximated solution on the domain
% 
% [graduh] = gradApproxSol2d(w,g,indVec,nodes,triangles,midNodes,triangleMidPts,shapeFn)
% 
% where w is the values of the nodes, g is the Dirichlet BC on Gamma_1, 
% nodes are the corrdinates of nodes, triangles records the vertices of the
% triangles, indVec gives the index between the order of all nodes and the 
% order of free nodes, triangleMidPts records the node number of midedge 
% points, midNodes has the coordinates of the midedge points, shapeFn 
% indicates the choice of linear or quadratic shape functions 
% Last update: Chung-min Lee May 11, 2021

function [graduh] = approxSolGrad2d(w,g,indVec,nodes,triangles,midNodes,triangleMidPts,shapeFn)

if (shapeFn == 1) % if linear shape functions are used
    graduh = @(x,y) gradLinearApprox2d(x,y,w,g,indVec,nodes,triangles);
else    % if quadratic shape functions are used
    graduh = @(x,y) gradQuadraticApprox2d(x,y,w,g,indVec,nodes,triangles,midNodes,triangleMidPts);
end


%%---------
% assemble the gradient of the linear shape functions on the whole domain

function gradz = gradLinearApprox2d(x,y,w,g,indVec,nodes,triangles)

noNodes = length(nodes(:,1));

% put all node values together
uu = zeros(noNodes,1);
uu(logical(indVec)) = w;
gamma1s = (indVec==0);
uu(gamma1s) = g(nodes(gamma1s,1),nodes(gamma1s,2));

% vertices of triangles
N1 = triangles(:,1);
N2 = triangles(:,2);
N3 = triangles(:,3);

% no of triangles
notri = length(N1);

[lxr,lxc] = size(x);
lx = lxr*lxc;
gradz = zeros(lx,2);

if (lx > notri)
    ckNan = zeros(lx,1);
    for i = 1:notri
        gradpsi1 = reshape(shapeFnGrad2d(1,x(:),y(:),nodes(N1(i),1),nodes(N1(i),2),nodes(N2(i),1),nodes(N2(i),2),nodes(N3(i),1),nodes(N3(i),2),1),lx,2);
        gradpsi2 = reshape(shapeFnGrad2d(2,x(:),y(:),nodes(N1(i),1),nodes(N1(i),2),nodes(N2(i),1),nodes(N2(i),2),nodes(N3(i),1),nodes(N3(i),2),1),lx,2);
        gradpsi3 = reshape(shapeFnGrad2d(3,x(:),y(:),nodes(N1(i),1),nodes(N1(i),2),nodes(N2(i),1),nodes(N2(i),2),nodes(N3(i),1),nodes(N3(i),2),1),lx,2);
        ckNan = ckNan+logical(sum((gradpsi1 | gradpsi2 | gradpsi3),2));
        gradz = gradz + (uu(N1(i))*ones(lx,2)).*gradpsi1 + (uu(N2(i))*ones(lx,2)).*gradpsi2 + (uu(N3(i))*ones(lx,2)).*gradpsi3;
    end
    gradz(logical(ckNan),:) = gradz(logical(ckNan),:)./(ckNan(logical(ckNan))*ones(1,2));
    gradz(~ckNan,:) = nan(sum(~ckNan),2);
else
     for i = 1:lx
         gradpsi1 = reshape(shapeFnGrad2d(1,x(i),y(i),nodes(N1,1),nodes(N1,2),nodes(N2,1),nodes(N2,2),nodes(N3,1),nodes(N3,2),1),notri,2);
         gradpsi2 = reshape(shapeFnGrad2d(2,x(i),y(i),nodes(N1,1),nodes(N1,2),nodes(N2,1),nodes(N2,2),nodes(N3,1),nodes(N3,2),1),notri,2);
         gradpsi3 = reshape(shapeFnGrad2d(3,x(i),y(i),nodes(N1,1),nodes(N1,2),nodes(N2,1),nodes(N2,2),nodes(N3,1),nodes(N3,2),1),notri,2);
         ckNan = max(sum((gradpsi1 | gradpsi2 | gradpsi3),1));
         if (ckNan)
            gradz(i,:) = sum((uu(N1)*ones(1,2)).*gradpsi1 + (uu(N2)*ones(1,2)).*gradpsi2 + (uu(N3)*ones(1,2)).*gradpsi3,1)./ckNan;
         else
             gradz(i,:) = [nan nan];
         end
     end
end


%%------------------
% assemble gradient of the quadratic shape functions on the whole domain

function gradz = gradQuadraticApprox2d(x,y,w,g,indVec,nodes,triangles,midNodes,triangleMidPts)

noVer = length(nodes(:,1));
noNodes = noVer +length(midNodes(:,1));

% put all node values together
uu = zeros(noNodes,1);
uu(logical(indVec)) = w;
gamma1s = find(indVec==0);
gamma1nodes = find(gamma1s <= noVer);
gamma1midNodes = setdiff([1:length(gamma1s)],gamma1nodes);
uu(gamma1s(gamma1nodes)) = g(nodes(gamma1s(gamma1nodes),1),nodes(gamma1s(gamma1nodes),2));
uu(gamma1s(gamma1midNodes)) = g(midNodes(gamma1s(gamma1midNodes)-noVer,1),midNodes(gamma1s(gamma1midNodes)-noVer,2));

% vertices and midpoints of triangles
N1 = triangles(:,1);
N2 = triangles(:,2);
N3 = triangles(:,3);
N4 = triangleMidPts(:,1);
N5 = triangleMidPts(:,2);
N6 = triangleMidPts(:,3);


% no of triangles
notri = length(N1);

[lxr,lxc] = size(x);
lx = lxr*lxc;
gradz = zeros(lx,2);

if (lx > notri)
    ckNan = zeros(lx,1);
    for i = 1:notri
        gradpsi1 = reshape(shapeFnGrad2d(1,x(:),y(:),nodes(N1(i),1),nodes(N1(i),2),nodes(N2(i),1),nodes(N2(i),2),nodes(N3(i),1),nodes(N3(i),2),2),lx,2);
        gradpsi2 = reshape(shapeFnGrad2d(2,x(:),y(:),nodes(N1(i),1),nodes(N1(i),2),nodes(N2(i),1),nodes(N2(i),2),nodes(N3(i),1),nodes(N3(i),2),2),lx,2);
        gradpsi3 = reshape(shapeFnGrad2d(3,x(:),y(:),nodes(N1(i),1),nodes(N1(i),2),nodes(N2(i),1),nodes(N2(i),2),nodes(N3(i),1),nodes(N3(i),2),2),lx,2);
        gradpsi4 = reshape(shapeFnGrad2d(4,x(:),y(:),nodes(N1(i),1),nodes(N1(i),2),nodes(N2(i),1),nodes(N2(i),2),nodes(N3(i),1),nodes(N3(i),2),2),lx,2);
        gradpsi5 = reshape(shapeFnGrad2d(5,x(:),y(:),nodes(N1(i),1),nodes(N1(i),2),nodes(N2(i),1),nodes(N2(i),2),nodes(N3(i),1),nodes(N3(i),2),2),lx,2);
        gradpsi6 = reshape(shapeFnGrad2d(6,x(:),y(:),nodes(N1(i),1),nodes(N1(i),2),nodes(N2(i),1),nodes(N2(i),2),nodes(N3(i),1),nodes(N3(i),2),2),lx,2);
        ckNan = ckNan+logical(sum((gradpsi1 | gradpsi2 | gradpsi3 | gradpsi4 | gradpsi5 | gradpsi6),2));
        gradz = gradz + (uu(N1(i))*ones(lx,2)).*gradpsi1 + (uu(N2(i))*ones(lx,2)).*gradpsi2 + (uu(N3(i))*ones(lx,2)).*gradpsi3 +...
            (uu(N4(i))*ones(lx,2)).*gradpsi4 + (uu(N5(i))*ones(lx,2)).*gradpsi5 + (uu(N6(i))*ones(lx,2)).*gradpsi6;
    end
    gradz(logical(ckNan),:) = gradz(logical(ckNan),:)./(ckNan(logical(ckNan))*ones(1,2));
    gradz(~ckNan,:) = nan(sum(~ckNan),2);
else
     for i = 1:lx
         gradpsi1 = reshape(shapeFnGrad2d(1,x(i),y(i),nodes(N1,1),nodes(N1,2),nodes(N2,1),nodes(N2,2),nodes(N3,1),nodes(N3,2),2),notri,2);
         gradpsi2 = reshape(shapeFnGrad2d(2,x(i),y(i),nodes(N1,1),nodes(N1,2),nodes(N2,1),nodes(N2,2),nodes(N3,1),nodes(N3,2),2),notri,2);
         gradpsi3 = reshape(shapeFnGrad2d(3,x(i),y(i),nodes(N1,1),nodes(N1,2),nodes(N2,1),nodes(N2,2),nodes(N3,1),nodes(N3,2),2),notri,2);
         gradpsi4 = reshape(shapeFnGrad2d(4,x(i),y(i),nodes(N1,1),nodes(N1,2),nodes(N2,1),nodes(N2,2),nodes(N3,1),nodes(N3,2),2),notri,2);
         gradpsi5 = reshape(shapeFnGrad2d(5,x(i),y(i),nodes(N1,1),nodes(N1,2),nodes(N2,1),nodes(N2,2),nodes(N3,1),nodes(N3,2),2),notri,2);
         gradpsi6 = reshape(shapeFnGrad2d(6,x(i),y(i),nodes(N1,1),nodes(N1,2),nodes(N2,1),nodes(N2,2),nodes(N3,1),nodes(N3,2),2),notri,2);
         ckNan = max(sum((gradpsi1 | gradpsi2 | gradpsi3 | gradpsi4 | gradpsi5 | gradpsi6 ),1));
         if (ckNan)
             gradz(i,:) = sum((uu(N1)*ones(1,2)).*gradpsi1 + (uu(N2)*ones(1,2)).*gradpsi2 + (uu(N3)*ones(1,2)).*gradpsi3 + ...
                (uu(N4)*ones(1,2)).*gradpsi4 + (uu(N5)*ones(1,2)).*gradpsi5 + (uu(N6)*ones(1,2)).*gradpsi6,1)./ckNan;
         else
             gradz(i,:) = [nan nan];
         end
     end
end
 