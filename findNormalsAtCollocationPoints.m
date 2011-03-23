function [ collocNormals ] = findNormalsAtCollocationPoints(nPts, ne, elRange, bsFnConn, collocPts, dispConn, controlPts)

global knotVec

collocNormals=zeros(nPts,2,2);    % normals of collocation points

tol =eps;

for c=1:nPts
    srcXi_param=collocPts(c);
    for element=1:ne
        range=elRange(element,:);
        glbBsFnConn=bsFnConn(element,:);
        dElConn=dispConn(element,:);        % element connectivity for displacement nodes
        elcoords=controlPts(dElConn,1:2);        % coordinates of element nodes
        
        if (srcXi_param <= range(2)) && (srcXi_param >= range(1) || (element==ne && srcXi_param==0))
            
            % find the normals 'just' before and after the collocation point
            if abs(srcXi_param-range(2)) < tol
                collocNormals(c,:,1)=getNormal(elcoords, srcXi_param-tol, glbBsFnConn);
            elseif abs(srcXi_param-range(1)) < tol
                collocNormals(c,:,2)=getNormal(elcoords, srcXi_param, glbBsFnConn);
            elseif element==ne && srcXi_param==0
                srcXi_param=1-tol;
                collocNormals(c,:,1)=getNormal(elcoords, max(knotVec)-tol, glbBsFnConn);
            else
                collocNormals(c,:,1)=getNormal(elcoords, srcXi_param, glbBsFnConn);
                collocNormals(c,:,2)=collocNormals(c,:,1);
            end
        end
    end
    
end

