function [ zterm ] = integrateExactTractionTerm(ngp, elcoords, bsFnConn, collocCoords, srcXi_param, range, traction)

global controlPts knotVec p

% This is used to calculate the exact traction terms in the infinite plate
% model. We simply return a 2x1 vector that is put into the zvector
zterm=zeros(2,1);
[gpt,gwt]=lgwt(ngp,-1,1);
jacob_param=(range(2)-range(1)) / 2;    % jacobian from parent to parameter space

numBasisFns=length(bsFnConn);
N=zeros(1,numBasisFns); dN=zeros(1,numBasisFns);
    
% do we have a singular integral?

if (srcXi_param <= range(2)) && (srcXi_param >= range(1))
    
    % --------------------------------------
    % -------- Telles' transformation ------
    % --------------------------------------
    
    srcXi=convertToParentCoordSpace(srcXi_param, range);
    xiStar=srcXi^2 - 1;
    gamBar=nthroot( srcXi*xiStar + abs(xiStar) ,3) + nthroot(srcXi*xiStar - abs(xiStar),3) + srcXi;
    
    for pt=1:ngp    % integrate using Gaussian quadrature
        gXi=gpt(pt);
        xi=((gXi-gamBar).^3 + gamBar*(gamBar.^2+3)) / (1+3*gamBar.^2);
        xi_param=convertToParamSpace( xi, range);    % get the gauss point in parameter space
        
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
        jacobTelles=3* ((gXi-gamBar).^2) / (1+3*gamBar.^2);
        coords=N*controlPts(bsFnConn,1:2);
        exactTraction=getExactStressTensor( traction, coords )*normals';
        zterm=zterm + Utemp*exactTraction*jacob*jacobTelles*gwt(pt);
        
    end
    
else
    
    % --------------------------------------
    % -------- Gaussian Quadrature ---------
    % --------------------------------------
    
    for pt=1:ngp                            % loop over gauss points
        xi=gpt(pt); wt=gwt(pt);
        xi_param=convertToParamSpace( xi, range);    % get the gauss point in parameter space
        
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
        coords=N*controlPts(bsFnConn,1:2);
        exactTraction=getExactStressTensor( traction, coords )*normals';
        zterm=zterm + Utemp*exactTraction*jacob*wt;
    end
end


end

