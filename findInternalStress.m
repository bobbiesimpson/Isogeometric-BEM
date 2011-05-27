function [ stress ] = findInternalStress( displacements, tractions, collocCoords )

% find the stress at an internal point

global p knotVec controlPts elRange bsFnConn tracConn dispConn

stress = zeros(2,2);
ngp = 50;
[gpt gwt]=lgwt(ngp,-1,1);

for element=1:ne
    range=elRange(element,:);
    glbBsFnConn=bsFnConn(element,:);    % connectiviy of NURBS basis fns
    dElConn=dispConn(element,:);        % element connectivity for displacement nodes
    tElConn=tracConn(element,:);        % element connectivity for traction nodes
    
    elcoords=controlPts(dElConn,1:2);        % coordinates of element nodes
    dispSctrX = 2 * dElConn - 1;    dispSctrY = 2 * dElConn;
    tracSctrX = 2 * tElConn - 1;    tracSctrY = 2 * tElConn;
    
    jacob_param=(range(2)-range(1)) / 2;    % jacobian from parent to parameter space
    
    for pt=1:ngp    % integrate using Gaussian quadrature
        xi_param=convertToParamSpace( gpt(pt), range);    % get the gauss point in parameter space
        
        N=zeros(1,length(glbBsFnConn)); dN = zeros(1,length(glbBsFnConn);
        for lclBasis=1:length(glbBsFnConn)
            i=glbBsFnConn(lclBasis);
            [N(lclBasis) dN(lclBasis)] = NURBSbasis(i, p, xi_param, knotVec, controlPts(:,3)' );
        end
        keyboard
        
        trac = N * [tractions(tracSctrX) tractions(tracSctrY)];
        disp = N * [displacements(dispSctrX) displacements(dispSctrY)];
        
        [jacob_xi, normal, r, dr, drdn] = getKernelParameters( elcoords, collocCoords, N, dN);
        jacob=jacob_xi*jacob_param;     % the final jacobian we use
        
        tempS = zeros(2,2,2); tempD = zeros(2,2,2);
        
        %             for i=1:2
        %                 for j=1:2
        %                     for k=1:2
        %                         [tempS(i,j,k) tempD(i,j,k)] =  TBIEkernels(i,j, k, r, dr, drdn, normal );
        %                     end
        %                     stress(i,j) = stress(i,j) + (tempD(i,j,1) * trac(1) +  tempD(i,j,2) * trac(2)) ...
        %                                               - (tempS(i,j,1) * disp(1) + tempS(i,j,2) * disp(2)) * gwt(pt) * jacob;
        %
        %                 end
        %             end
        
        
    end

end

end

