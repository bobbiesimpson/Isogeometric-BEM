function plotDeformedProfile( displacement, nPts, controlPts, infTraction )

% pass in a vector of control point displacements, and we can plot the
% deformed profile using NURBS 

global p knotVec 

factor=10; % for exaggerating displacements

newControlPts=controlPts(1:nPts,1:2)+factor*[displacement(1:2:end) displacement(2:2:end)];
newControlCoords=[newControlPts; newControlPts(1,:)];   % add first entry again (for plotting)

numPts=400;
xi=linspace(0,max(knotVec),numPts);

extrapPoints = zeros(numPts,2); x = zeros(numPts,1); y=zeros(numPts,1);
for point=1:numPts
    [extrapPoints(point,1)] = NURBSinterpolation(xi(point), p, knotVec, newControlPts(:,1)', controlPts(:,3)');
    [extrapPoints(point,2)] = NURBSinterpolation(xi(point), p, knotVec, newControlPts(:,2)', controlPts(:,3)');  
    x(point) = NURBSinterpolation(xi(point), p, knotVec, controlPts(:,1)', controlPts(:,3)');
    y(point) = NURBSinterpolation(xi(point), p, knotVec, controlPts(:,2)', controlPts(:,3)');
end

% N=zeros(numPts,numBasisFns);
% for i=1:numBasisFns
% 
%     for point=1:numPts
%         N(point,i)=NURBSbasis(i,p,xi(point),knotVec,controlPts(:,3));
%     end
% 
% end
% 
% x=N*controlPts(:,1);
% y=N*controlPts(:,2);
r=sqrt(x.^2+y.^2)';
theta=atan2(y,x)';

% Bit of hardcoding here - dimensions of plate
% a=1;

%disp = getExactDisplacements( r, theta, infTraction );

% points=zeros(numPts,2);
% points(:,1)=x+disp(:,1).*factor;
% points(:,2)=y+disp(:,2).*factor;

figure(2); hold on
plot(extrapPoints(:,1), extrapPoints(:,2), 'k--')
hold off

%% Plate with hole plotting

% figure(4); hold on
% axis([-4.5 0.5 -0.5 4.5])
% axis equal
% plot(points(:,1), points(:,2), 'k:', extrapPoints(:,1), extrapPoints(:,2),'k-', newControlCoords(:,1), newControlCoords(:,2), 'ko')
% legend('exact','isoBEM', 'isoBEM control points')

save 'dat_files/deformedPoints.dat' newControlCoords -ASCII
save 'dat_files/deformedProfile.dat' extrapPoints -ASCII
%save 'dat_files/exactDeformedProfile.dat' points -ASCII

end

