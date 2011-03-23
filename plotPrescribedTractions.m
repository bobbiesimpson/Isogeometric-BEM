function  plotPrescribedTractions( nonZeroXTracDOFs, nonZeroYTracDOFs, tracDispConnDOF,presXTracs, presYTracs)

% plot the tractions that are prescribed
global controlPts

xPresTracDof = tracDispConnDOF(nonZeroXTracDOFs);
yPresTracDof = tracDispConnDOF(nonZeroYTracDOFs);

xPresPts = (xPresTracDof + 1) ./2;
yPresPts = (yPresTracDof) ./2;

%xAndYPresPts = intersect(xPresPts, yPresPts);   % the points which have both x and y prescribed tractions

figure(1); hold on
quiver(controlPts(xPresPts,1), controlPts(xPresPts,2), presXTracs, presYTracs);
hold off

end

