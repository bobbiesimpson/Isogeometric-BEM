% This routine uses the C-function to construct the B-spline functions

close all
clear all
clc


knot=[0 0 0 1 2 2 2 ];
xi=0:0.01:max(knot);
numPts=numel(xi);
m = numel(knot)-1;
p = 2;
n = m - p;
weights = ones(1,n);
weights(2) = 1;

points = [0 10; 2 7; 7 7; 10 8];
interpolatedPoints = zeros(numPts,2);
% weights(4) =1;
% weights(6) =1;

BsplineVals = zeros(numPts,n);
NURBSderivs = zeros(numPts,n);

tic
for i=1:n
    for c=1:numPts
        [BsplineVals(c,i+1) NURBSderivs(c,i+1)] = NURBSbasis(i, p, xi(c), knot, weights);
    end
end
toc

figure(1)
plot(xi, BsplineVals)

figure(2)
plot(xi, NURBSderivs)

for c=1:numPts
    [interpolatedPoints(c,1)] = NURBSinterpolation(xi(c), p, knot, points(:,1)', weights);
    [interpolatedPoints(c,2)] = NURBSinterpolation(xi(c), p, knot, points(:,2)', weights);      
end

figure(3)
plot(interpolatedPoints(:,1), interpolatedPoints(:,2), 'k-', points(:,1), points(:,2), 'ko')


%% Now try some knot insertion to create two separate curves

% controlPts = [points weights'];
% knotInsert = [1 1];
% [ newKnot newControlPts ] = knotInsertion( knot, controlPts, knotInsert, p );
% 
% newWeights = newControlPts(:,3);
% 
% % create a crack
% 
% crackEntrance = [NURBSinterpolation(1, p, newKnot, newControlPts(:,1)', newWeights)...
%                  NURBSinterpolation(1, p, newKnot, newControlPts(:,2)', newWeights)];
% %newControlPts(4,1:2) = 6;
% firstKnotEnd = find(newKnot==1, 1, 'last' );
% secondKnotStart = find(newKnot==1, 1 );
% 
% firstKnot = newKnot(1:firstKnotEnd); 
% index  = newKnot(firstKnotEnd);
% 
% crackPoints = [[newControlPts(3,1:2) 1];
%                 [newControlPts(3,1:2) + 2 1];
%                 [newControlPts(3,1:2) + 4 1] ];
% crackPoints = [crackPoints;
%                 crackPoints(end-1:-1:1,:)];
%             
% newControlPts = [newControlPts(1:3,:);
%                  crackPoints;
%                  newControlPts(4:end,:)];
%              
% newKnot = [newKnot(1:firstKnotEnd) (index + 1)*ones(1,2) newKnot(secondKnotStart:end) + 2];
% 
% xi=0:0.01:max(newKnot);
% numPts=numel(xi);
% interpolatedPoints = zeros(numPts,2);
% for c=1:numPts
%     [interpolatedPoints(c,1)] = NURBSinterpolation(xi(c), p, newKnot, newControlPts(:,1)', newWeights);
%     [interpolatedPoints(c,2)] = NURBSinterpolation(xi(c), p, newKnot, newControlPts(:,2)', newWeights);      
% end
% 
% figure(3)
% hold on
% plot(interpolatedPoints(:,1), interpolatedPoints(:,2), 'r--', newControlPts(:,1), newControlPts(:,2), 'ro')
% hold off
% 
% n = size(newControlPts,1) + p + 1;
% tic
% for i=1:n
%     for c=1:numPts
%         [BsplineVals(c,i+1) NURBSderivs(c,i+1)] = NURBSbasis(i, p, xi(c), newKnot, newWeights);
%     end
% end
% toc
% 
% figure(4)
% plot(xi, BsplineVals)

