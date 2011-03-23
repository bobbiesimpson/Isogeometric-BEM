function [ Tij Uij ] = DBIEkernels(i,j,r, dr, drdn, normal )
% calculate the traction kernel 

global const1 const2 const3 const4

Tij= 1/r * const3 * ( -drdn *  (const4*(i==j) + 2*( dr(i)*dr(j) )) ...
            + const4*( dr(i)*normal(j) - dr(j)*normal(i) ));
        
Uij=const1*( dr(i)*dr(j) + (i==j) * ( const2*log(1/r) ) );        
    
end

