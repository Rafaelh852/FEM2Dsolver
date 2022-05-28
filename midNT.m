function [midNodes, triangleMidPts] = midNT(nodes,triangles)

    p1 = nodes(triangles(:,1),1:2);
    p2 = nodes(triangles(:,2),1:2);
    p3 = nodes(triangles(:,3),1:2);
    
    mp1 = .5*(p1 + p2);
    mp2 = .5*(p2 + p3);
    mp3 = .5*(p1 + p3);
    
    %all mid nodes, repeats included
    amn = [mp1; mp2; mp3];
    %unique midNodes i.e midNodes
    umn = unique(amn,"rows");
    
    
    %find(sum(mp3(1,1:2)==umn,2)==2) finds the position in the midnodes
    %[mp1 mp2 mp3] contains the nodes of the midpoints of each tri
    
    m = 1;
    for t = [mp1 mp2 mp3]'
        T(m,:) = [find(sum(t(1:2)'==umn,2)==2), find(sum(t(3:4)'==umn,2)==2),find(sum(t(5:6)'==umn,2)==2)];
        m = m+1;
    end
    
    midNodes = umn;
    triangleMidPts = T + length(nodes(:,1));  


end
