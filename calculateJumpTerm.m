function [ Cij ] = calculateJumpTerm( n1, n2)
% Caclulate the jump term expliciting at a collocation point.
% It depends entirely on the geometry and the following code is based on
% the formula given in the paper by Guiggiani and Casalini

% n1: the normal of the 'first' element
% n2: the normals of the 'second' element

global mu

Cij=zeros(2,2);
t=zeros(2,2);

tempTang=cross([n1 0],[0 0 1]);
t(1,:)=tempTang(1:2);
tempTang=cross([n2 0],[0 0 1]);
t(2,:)=tempTang(1:2);

dotProd=dot(n1,n2);
if (dotProd - 1) >= eps
    theta = 0;
else
    theta=acos(dot(n1,n2));
end

if n1(1) * n2(2) > n1(2) * n2(1)
    thetaBar=pi-theta;
else
    thetaBar=pi + theta;
end

xVector=[1 0];
elmTheta=zeros(1,2);

for i=1:2
    theta=acos(dot(xVector,t(i,:)));
    if xVector(1)*t(i,2) > xVector(2)*t(i,1)
        elmTheta(i)=theta;
    else
        elmTheta(i)=2*pi-theta;
    end
end

term1 = 1 / ( 8 * pi * ( 1 - mu ) );
term2 = 4 * ( 1 - mu ) * thetaBar;
term3 = sin( 2 * elmTheta(1) ) - sin( 2 * elmTheta(2) );
term4 = cos( 2 * elmTheta(2) ) - cos( 2 * elmTheta(1) );

Cij(1,1) = ( term2 + term3 ) * term1;		
Cij(1,2) =  term1 * term4;
Cij(2,1)=Cij(1,2);
Cij(2,2) = ( term2 - term3 ) * term1;
    

end

