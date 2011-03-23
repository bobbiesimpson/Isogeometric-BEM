function [ Hsubmatrix ] = integrateHsubmatrixSST( ngp, elcoords, bsFnConn, collocCoords, srcXi_param, range, jumpTerm)
% A function which calculate the H matrix for the singular case using the
% subtraction of singularity technique by Guiggiani and Casalini

% We pass in ....

% ngp:              # of gauss points
% elcoords:         the coordinates of the element we are on
% bsFnConn          the global basis Fns which are non zero in this element
% collocCoords:     the coordinates of the collocation points (do we need this?)
% srcXi_param:      the coordinate of our source point in parameter space
% range:            the range of our 'element' [xi_i xi_i+1]

global const3 const4 p knotVec controlPts

[gpt gwt]=lgwt(ngp,-1,1);               % get the gauss points

Hsubmatrix=zeros(2,6);

numBasisFns=length(bsFnConn);
srcN=zeros(1,numBasisFns); N=zeros(1,numBasisFns);
srcdN=zeros(1,numBasisFns); dN=zeros(1,numBasisFns);

if srcXi_param==range(1)            % Annoying, but I need to evaluate the parameters 
    nudgedXi=srcXi_param+eps;    % at the source point a small distance away from the actual point
elseif srcXi_param==range(2)
    nudgedXi=srcXi_param-eps;
else 
    nudgedXi=srcXi_param;
end
srcXi=convertToParentCoordSpace(srcXi_param, range);

for c=1:length(bsFnConn)
    [srcN(c) srcdN(c)]=NURBSbasis( bsFnConn(c), p, nudgedXi, knotVec, controlPts(:,3)' );
end

hterm=[0 -const4*const3; const4*const3 0];
htermMatrix=[hterm.*srcN(1) hterm.*srcN(2) hterm.*srcN(3)];

jacob_param=(range(2)-range(1)) / 2;    % jacobian from parent to parameter space

for pt=1:ngp    % integrate using Gaussian quadrature
    xi=gpt(pt);
    xi_param=convertToParamSpace( xi, range);    % get the gauss point in parameter space
            
    for lclBasis=1:numBasisFns
        i=bsFnConn(lclBasis);
        [N(lclBasis) dN(lclBasis)]=NURBSbasis(i, p, xi_param, knotVec, controlPts(:,3)' );
    end
    
    [jacob_xi,normals, r, dr, drdn] = getKernelParameters( elcoords, collocCoords, N, dN );
    jacob=jacob_xi*jacob_param;     % the final jacobian we use
    
    Ttemp=zeros(2,2);
    
    for i=1:2
        for j=1:2
            Ttemp(i,j)=DBIEkernels(i,j,r,dr,drdn,normals);
        end
    end
    
%     figure(3)
%     hold on
%     plot(xi, N(1)*Ttemp(1,2)*jacob, 'ko-')
%     hold off
    
    tempMatrix=[N(1)*Ttemp N(2)*Ttemp N(3)*Ttemp]...
        *(xi-srcXi)*jacob*gwt(pt);
    tempMatrix=tempMatrix-htermMatrix;
    tempMatrix=tempMatrix./(xi-srcXi);
    Hsubmatrix=Hsubmatrix + tempMatrix;
        
end

% and add the analytical integral on at the end
if abs(srcXi)<(1-100*eps)
    Hsubmatrix=Hsubmatrix+htermMatrix*log(abs((1-srcXi)/(1+srcXi)));
else
    jacob_xi=getKernelParameters(elcoords, collocCoords, srcN, srcdN);
    jacob_s=jacob_xi*jacob_param;
    beta_m=1/jacob_s;
    Hsubmatrix=Hsubmatrix+htermMatrix*log(abs(2/beta_m))*sign(xi-srcXi);
end

% and add the jump terms at the end
if abs(srcXi_param-range(2)) > 100*eps
    jumpMatrix=[jumpTerm.*srcN(1) jumpTerm.*srcN(2) jumpTerm.*srcN(3)];
    Hsubmatrix=Hsubmatrix+jumpMatrix;
end

end

