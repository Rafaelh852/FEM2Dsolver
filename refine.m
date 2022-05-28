% this program refines a triangulation using standard refinement
% [nodes,triangles,edges,bdyNode,bdyEdge,curveEdge]=refine(nodes,triangles,
% edges,bdyNode,bdyEdge,curveEdge,bdyFn,bdyFnder)
% 
% input and output:
% nodes - 2-column array that stores the (x,y)-coordinates of nodes
% trangles - 3-column array that stores the nodes in each triangle
% edges - 2-column array that stores the end nodes of edges
% bdyNode - 1-column vector that indicates if a node is a boundary node (1
% = true, 0 = false)
% bdyEdge - 1-column array that indicates if an edge is a boundary edge (
% 1= true, 0 = false)
% curveEdge - 1-column array that indicates if an edge is an approximation
% of a curve (1 = true, 0 = false)
% bdyFn - a function handle represents the function that describ the
% boundary
% bdyFnder - the derivative of the boundary function
%
% Chung-min Lee Mar 31, 2011

function [nodes,triangles,edges,bdyNode,bdyEdge,curveEdge]=refine(nodes,triangles,edges,bdyNode,bdyEdge,curveEdge,bdyFn,bdyFnder)

[noTri] = length(triangles(:,1)); % no of triangles
[noNode] = length(nodes(:,1)); % no of nodes
parentNodes = [];   % track where the midpoint is from


for k = 1:noTri
    node1 = triangles(k,1);
    node2 = triangles(k,2);
    node3 = triangles(k,3);
    
    % check the first edge
    inEdges = (sum((edges==node1)+(edges==node2),2)==2);
    if (sum(inEdges))
        edgeNo = find(inEdges); % which edge it is
        
        % find midpoint
        if (curveEdge(edgeNo) == 1) % if it is a curve edge
            [x1m,y1m] = midPtCurve(bdyFn,bdyFnder,nodes(node1,1),nodes(node2,1));
        else
            x1m = 0.5*(nodes(node1,1)+nodes(node2,1));
            y1m = 0.5*(nodes(node1,2)+nodes(node2,2));
        end
        [nodes,edges,bdyNode,bdyEdge,curveEdge] = updateTri(x1m,y1m,node1,node2,edgeNo,nodes,edges,bdyNode,bdyEdge,curveEdge);
        [numNode] = length(nodes(:,1));
        nodeNo4MidPt1 = numNode; % record the node number of the 1st midpoint
        parentNodes = [parentNodes;node1 node2]; % record the parents that produce this midpoints
    else % not in existing edges
        % find the node number for the 1st midpoint
        nodeNo4MidPt1 = noNode+find(sum((parentNodes==node1)+(parentNodes==node2),2)==2);
    end
    clear x1m y1m inEdges edgeNo
    
    % check the second edge
    inEdges = (sum((edges==node2)+(edges==node3),2)==2);
    if (sum(inEdges))
        edgeNo = find(inEdges); % which edge it is
        
        % find midpoint
        if (curveEdge(edgeNo) == 1) % if it is a curve edge
            [x2m,y2m] = midPtCurve(bdyFn,bdyFnder,nodes(node2,1),nodes(node3,1));
        else
            x2m = 0.5*(nodes(node2,1)+nodes(node3,1));
            y2m = 0.5*(nodes(node2,2)+nodes(node3,2));
        end
        [nodes,edges,bdyNode,bdyEdge,curveEdge] = updateTri(x2m,y2m,node2,node3,edgeNo,nodes,edges,bdyNode,bdyEdge,curveEdge);
        [numNode] = length(nodes(:,1));
        nodeNo4MidPt2 = numNode; % record the node number of the 2nd midpoint
        parentNodes = [parentNodes;node2 node3]; % record the parents that produce this midpoints
    else % not in existing edges
        % find the node number for the 2nd midpoint
        nodeNo4MidPt2 = noNode+find(sum((parentNodes==node2)+(parentNodes==node3),2)==2);
    end
    clear x2m y2m inEdges edgeNo
    

    % check the third edge
    inEdges = (sum((edges==node1)+(edges==node3),2)==2);
    if (sum(inEdges))
        edgeNo = find(inEdges); % which edge it is
        
        % find midpoint
        if (curveEdge(edgeNo) == 1) % if it is a curve edge
            [x3m,y3m] = midPtCurve(bdyFn,bdyFnder,nodes(node1,1),nodes(node3,1));
        else
            x3m = 0.5*(nodes(node1,1)+nodes(node3,1));
            y3m = 0.5*(nodes(node1,2)+nodes(node3,2));
        end
        [nodes,edges,bdyNode,bdyEdge,curveEdge] = updateTri(x3m,y3m,node1,node3,edgeNo,nodes,edges,bdyNode,bdyEdge,curveEdge);
        [numNode] = length(nodes(:,1));
        nodeNo4MidPt3 = numNode; % record the node number of the 3rd midpoint
        parentNodes = [parentNodes;node1 node3]; % record the parents that produce this midpoints
    else % not in existing edges
        % find the node number for the 3rd midpoint
        nodeNo4MidPt3 = noNode+find(sum((parentNodes==node1)+(parentNodes==node3),2)==2);
    end
    clear x3m y3m inEdges edgeNo
    
    % put the 3 new edges in (cnnecting 3 midpoints)
    edges = [edges;nodeNo4MidPt1 nodeNo4MidPt2;nodeNo4MidPt2 nodeNo4MidPt3;nodeNo4MidPt1 nodeNo4MidPt3];
    bdyEdge = [bdyEdge;0;0;0];
    curveEdge = [curveEdge;0;0;0];
    
    % replace this triangle in the triangle list by 4 new triangles
    triangles(k,:) = [node1 nodeNo4MidPt1 nodeNo4MidPt3];
    triangles = [triangles; node2 nodeNo4MidPt2 nodeNo4MidPt1];
    triangles = [triangles; node3 nodeNo4MidPt3 nodeNo4MidPt2];
    triangles = [triangles; nodeNo4MidPt1 nodeNo4MidPt2 nodeNo4MidPt3];
end


%%------------
% the function that updates nodes, edges, bdyNode, bdyEdge, curveEdge
%
% Chung-min Lee Oct 28, 2008

function [nodes,edges,bdyNode,bdyEdge,curveEdge] = updateTri(xm,ym,parentNode1,parentNode2,edgeNo,nodes,edges,bdyNode,bdyEdge,curveEdge)

% add the midpoint to nodes
nodes = [nodes; xm ym];
[row] = length(nodes(:,1));

% replace the original edge by 2 new edges
edges(edgeNo,:) = [parentNode1 row];
edges = [edges;parentNode2 row];

% if the parent nodes are on a bdyEdge, then this midpoint is a boundary
% node and both new edges are boundary edges
if (bdyEdge(edgeNo) == 1)
    bdyNode = [bdyNode;1];
    bdyEdge = [bdyEdge;1];
else
    bdyNode = [bdyNode;0];
    bdyEdge = [bdyEdge;0];
end

% update curveEdge for two new edges
if (curveEdge(edgeNo) == 1)
    curveEdge = [curveEdge;1];
else
    curveEdge = [curveEdge;0];
end


