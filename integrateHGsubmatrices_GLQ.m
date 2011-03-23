function [ Hsubmatrix, Gsubmatrix ] = integrateHGsubmatrices_GLQ( ngp, elcoords, bsFnConn, collocCoords, range)
    % integrate the H and G submatrices using Gauss Legendre quadrature
    
    global controlPts p knotVec
    
    [gpt gwt]=lgwt(ngp,-1,1);
    
    Hsubmatrix=zeros(2,6);
    Gsubmatrix=zeros(2,6);
    
    numBasisFns=length(bsFnConn);
    N=zeros(1,numBasisFns); dN=zeros(1,numBasisFns);
    
    jacob_param=(range(2)-range(1)) / 2;    % jacobian from parent to parameter space
    
    for pt=1:ngp    % integrate using Gaussian quadrature
        
        xi_param=convertToParamSpace( gpt(pt), range);    % get the gauss point in parameter space

        for lclBasis=1:numBasisFns
            i=bsFnConn(lclBasis);
            [N(lclBasis) dN(lclBasis)]=NURBSbasis(i, p, xi_param, knotVec, controlPts(:,3)' );
        end
        
        [jacob_xi, normals, r, dr, drdn] = getKernelParameters( elcoords, collocCoords, N, dN );
        jacob=jacob_xi*jacob_param;     % the final jacobian we use
        
        Ttemp=zeros(2,2); Utemp=zeros(2,2);
        for i=1:2
            for j=1:2
                [Ttemp(i,j) Utemp(i,j)]=DBIEkernels(i,j,r,dr,drdn,normals);
            end
        end
        
        Hsubmatrix=Hsubmatrix + [N(1)*Ttemp N(2)*Ttemp N(3)*Ttemp]*jacob*gwt(pt);
        Gsubmatrix=Gsubmatrix + [N(1)*Utemp N(2)*Utemp N(3)*Utemp]*jacob*gwt(pt);
        
    end
end


