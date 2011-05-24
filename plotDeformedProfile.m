function plotDeformedProfile( displacement, NURBScurve)

% ------------------------------------------------------------
% ---------------- Plot deformed shape -----------------------
% ------------------------------------------------------------

% ----------------------------
% --------- inputs -----------
% ----------------------------

% displacement      :        a vector of all displacements at DISPLACEMENT NODES
% NURBScurve        :        a struct containing all the NURBS curve information

factor=10; % for exaggerating displacements
numPtsPerEl = 50;

numCurves = size(NURBScurve,2);

figure(2); hold on
for curve=1:numCurves
    ne = NURBScurve(curve).ne;
    for element=1:ne
        range = NURBScurve(curve).elRange(element,:);
        dispConn = NURBScurve(curve).dispConn(element,:);
        bsFnConn = NURBScurve(curve).bsFnConn(element,:);
        elCoords = NURBScurve(curve).controlPts(bsFnConn,1:2);
        
        xi_param = linspace(range(1), range(2),numPtsPerEl);
        
        sctrX = dispConn * 2 - 1;
        sctrY = dispConn * 2;
        
        deformedPts = zeros(numPtsPerEl,2);
        for pt=1:numPtsPerEl
            % -- interpolate the displacements -------
            [ shape geomShape ] = getDispAndGeomBasis( xi_param(pt), NURBScurve(curve), element );
            deformedDisp = shape * [displacement(sctrX) displacement(sctrY)];
            geometryCoord = geomShape * elCoords;

            deformedPts(pt,:) = geometryCoord + factor .* deformedDisp;
        end
        
        plot(deformedPts(:,1), deformedPts(:,2), 'k--')
    end
end
hold off

%save 'dat_files/exactDeformedProfile.dat' points -ASCII

end

