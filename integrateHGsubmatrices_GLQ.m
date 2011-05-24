function [ Hsubmatrix, Gsubmatrix ] = integrateHGsubmatrices_GLQ( ngp, NURBScurve, element, collocCoords)
    % integrate the H and G submatrices using Gauss Legendre quadrature
    
    [gpt gwt]=lgwt(ngp,-1,1);
    
    Hsubmatrix=zeros(2,6);
    Gsubmatrix=zeros(2,6);
    
    range = NURBScurve.elRange(element,:);
    bsFnConn = NURBScurve.bsFnConn(element,:);
    elcoords = NURBScurve.controlPts(bsFnConn,1:2);
    
    jacob_param=(range(2)-range(1)) / 2;    % jacobian from parent to parameter space
    
    for pt=1:ngp    % integrate using Gaussian quadrature
        
        xi_param=convertToParamSpace( gpt(pt), range);    % get the gauss point in parameter space
        
        [dispShape geomShape geomShapeDeriv] = getDispAndGeomBasis(xi_param, NURBScurve, element);
        
        [jacob_xi, normals, r, dr, drdn] = getKernelParameters( elcoords, collocCoords, geomShape, geomShapeDeriv );
        jacob=jacob_xi*jacob_param;     % the final jacobian we use
        
        Ttemp=zeros(2,2); Utemp=zeros(2,2);
        for i=1:2
            for j=1:2
                [Ttemp(i,j) Utemp(i,j)]=DBIEkernels(i,j,r,dr,drdn,normals);
            end
        end
        
        Hsubmatrix=Hsubmatrix + [dispShape(1)*Ttemp dispShape(2)*Ttemp dispShape(3)*Ttemp]*jacob*gwt(pt);
        Gsubmatrix=Gsubmatrix + [dispShape(1)*Utemp dispShape(2)*Utemp dispShape(3)*Utemp]*jacob*gwt(pt);
        
    end
end


