function [ n jacob ] = getNormal( elcoords, xi_param, bsFnConn)
% simply calculate the normal at a local coordinate on an element

global knotVec controlPts p

numBasisFns=length(bsFnConn);
N=zeros(1,numBasisFns); dN=zeros(1,numBasisFns);

for lclBasis=1:numBasisFns
    i=bsFnConn(lclBasis);
    [N(lclBasis) dN(lclBasis)]=NURBSbasis(i, p, xi_param, knotVec, controlPts(:,3)' );
end

dxydxi=dN*elcoords;             % the geometry derivatives
jacob=norm(dxydxi);             % jacobian
n=1/jacob * [ dxydxi(2) -dxydxi(1) ];


end
