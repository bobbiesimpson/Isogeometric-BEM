function [ collocNormals ] = findNormalsAtCollocationPoints(nPts, collocPts, NURBScurve)

global knotVec
global controlPts

collocNormals=zeros(nPts,2,2);    % normals of collocation points
numCurves = size(NURBScurve,2);

tol =eps;

for c=1:nPts                        % loop over all collocation poitns
    srcXi_param=collocPts(c);       % get the collocation point in parameter space
    
    for curve=1:numCurves  % loop over all curves
        
        ne = NURBScurve(curve).ne;
        knotVec = NURBScurve(curve).knotVec;
        elRange = NURBScurve(curve).elRange;
        bsFnConn = NURBScurve(curve).bsFnConn;
        controlPts = NURBScurve(curve).controlPts;
        for element=1:ne            % loop over all elements in the curve
            
            range=elRange(element,:);
            glbBsFnConn=bsFnConn(element,:);
            dElConn=bsFnConn(element,:); % element connectivity for displacement nodes
            elcoords=controlPts(dElConn,1:2);        % coordinates of element nodes
            
            if (srcXi_param <= range(2)) && (srcXi_param >= range(1) || (element==ne && srcXi_param==0 && curve==numCurves))
                
                % find the normals 'just' before and after the collocation
                % point
                if abs(srcXi_param-range(2)) < tol
                    collocNormals(c,:,1)=getNormal(elcoords, srcXi_param-tol, glbBsFnConn);
                elseif abs(srcXi_param-range(1)) < tol
                    collocNormals(c,:,2)=getNormal(elcoords, srcXi_param, glbBsFnConn);
                elseif element==ne && srcXi_param==0 && curve == numCurves
                    srcXi_param=1-tol;
                    collocNormals(c,:,1)=getNormal(elcoords, max(knotVec)-tol, glbBsFnConn);
                else
                    collocNormals(c,:,1)=getNormal(elcoords, srcXi_param, glbBsFnConn);
                    collocNormals(c,:,2)=collocNormals(c,:,1);
                end
            end
        end
    end
    
end

end

