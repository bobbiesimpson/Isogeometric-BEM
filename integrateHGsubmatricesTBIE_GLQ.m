function [ Hsubmatrix, Gsubmatrix ] = integrateHGsubmatricesTBIE_GLQ( ngp, NURBScurve, element, collocCoords, collocNormal)
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
        
        [jacob_xi, normal, r, dr, drdn] = getKernelParameters( elcoords, collocCoords, geomShape, geomShapeDeriv );
        jacob=jacob_xi*jacob_param;     % the final jacobian we use
        
        Sjk=zeros(2,2); Djk=zeros(2,2);
        tempS = zeros(2,2,2); tempD = zeros(2,2,2);
        
        for k=1:2
            for j=1:2
                for i=1:2
                    [tempS(j,k,i) tempD(j,k,i)] =  TBIEkernels(i,j, k, r, dr, drdn, normal );
                end
                Sjk(j,k) = tempS(j,k,1) * collocNormal(1) + tempS(j,k,2) * collocNormal(2);
                Djk(j,k) = tempD(j,k,1) * collocNormal(1) + tempD(j,k,2) * collocNormal(2);
                
            end
        end
         
        Hsubmatrix=Hsubmatrix + [dispShape(1)*Sjk dispShape(2)*Sjk dispShape(3)*Sjk]*jacob*gwt(pt);
        Gsubmatrix=Gsubmatrix + [dispShape(1)*Djk dispShape(2)*Djk dispShape(3)*Djk]*jacob*gwt(pt);
        
    end
end


