function [ newKnotVec newControlPts ] = knotInsertion( knotVec, controlPts, newKnots, p )
% Insert a knot into a knot vector and modify the control Points
% accordigly

    weightedPts=[controlPts(:,1).*controlPts(:,3) controlPts(:,2).*controlPts(:,3) controlPts(:,3)];
    oldKnotVec=knotVec;
    
    % and now redefine the control points
    for knot=1:length(newKnots)                                % loop over all the new knots
        newControlPts=zeros(length(weightedPts)+1,3);        % increase the array of control points
        knot_bar=newKnots(knot);
        kval=find(knotVec>knot_bar, 1 )-1;
        for i=1:length(newControlPts)
            if i<= (kval - p)
                alpha=1;
            elseif i>=(kval-p+1) && i<= kval
                alpha= ( newKnots(knot) - oldKnotVec(i) ) / ( oldKnotVec(i+p) - oldKnotVec(i) );
            else
                alpha=0;
            end
            newPoints=zeros(1,3);
            if i~=1
                newPoints=(1-alpha).*weightedPts(i-1,:);
            end
            if i~= length(newControlPts)
                newPoints=newPoints + alpha * weightedPts(i,:);
            end
            newControlPts(i,:)=newPoints;
        end
        weightedPts=newControlPts;
        knotVec=sort([knotVec knot_bar]);
        oldKnotVec=knotVec;
    end
    
    newControlPts=[weightedPts(:,1)./weightedPts(:,3) weightedPts(:,2)./weightedPts(:,3) weightedPts(:,3)];
    newKnotVec = knotVec;

end

