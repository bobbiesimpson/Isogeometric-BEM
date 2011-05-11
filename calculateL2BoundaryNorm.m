function [ L2relNorm ] = calculateL2BoundaryNorm( displacement, dispConn, bsFnConn, traction )

% calculate the L2 norm for disps around the boundary

global controlPts knotVec p elRange

ne=size(dispConn,1);

ngp=6;
[gpt gwt]=lgwt(ngp,-1,1);

L2Norm=0;
L2exactNorm=0;

numBasisFns=size(bsFnConn,2);
N=zeros(1,numBasisFns); dN=zeros(1,numBasisFns);

for element=1:ne
    elConn=dispConn(element,:);
    elCoords=controlPts(elConn,1:2);
    uxDof=elConn*2-1;
    uyDof=elConn*2;
    
    range=elRange(element,:);
    jacob_param=(range(2)-range(1)) / 2;    % jacobian from parent to parameter space
    
    for pt=1:ngp
        xi_param=convertToParamSpace( gpt(pt), range);
        
        for lclBasis=1:numBasisFns
            i=bsFnConn(element,lclBasis);
            [N(lclBasis) dN(lclBasis)]=NURBSbasis(i, p, xi_param, knotVec, controlPts(:,3)' );
        end
        
        dxydxi=dN*elCoords;             % the geometry derivatives
        jacobXi=norm(dxydxi);             % jacobian
        jacob=jacobXi*jacob_param;

        pointCoords=N*elCoords;
        r=norm(pointCoords);
        theta=atan2(pointCoords(2),pointCoords(1));
        exactDisp=getExactDisplacements(r,theta,traction);
        
        approxDisp=N*[displacement(uxDof) displacement(uyDof)];
        
        L2Norm=L2Norm + norm(approxDisp-exactDisp)*gwt(pt)*jacob;
        L2exactNorm=L2exactNorm + norm(exactDisp) * gwt(pt)*jacob;
    end
end

L2relNorm=L2Norm/L2exactNorm;

end

