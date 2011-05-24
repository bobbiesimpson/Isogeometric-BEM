function [ mirrorRange] = getMirrorElement( elRange, crackRange )

% -------- Get the range of the mirror element if we are on a crack ----

% ----- The purpose of this function is to pass in the crack boundary ----
% ----- range and the current element range. If we find the current ------
% ----- element lies on the crack, we return a range corresponding -------
% ----------------- to the element on the opposite face ------------------

% ---- First check if element lies outside crack --------
if ~(elRange(1) >= crackRange(1,1) && elRange(2) <= crackRange(2,2))
    mirrorRange=[0 0];

% ---- Otherwise we are on the crack. We can find a valid mirror element
else
    crackTip = crackRange(1,2);
    % --- we are on the upper crack surface -------
    if elRange(1) >= crackRange(1,1) && elRange(2) <= crackRange(1,2)
        crackTip = crackRange(1,2);
        lowerLimit = elRange(2) + 2 * (crackTip - elRange(2));
        upperLimit = lowerLimit + elRange(2) - elRange(1);
        
    % --- we are on the lower crack surface ------
    else
        upperLimit = elRange(1) - 2*(elRange(1) - crackTip);
        lowerLimit = upperLimit - (elRange(2) - elRange(1));
    end
    mirrorRange = [lowerLimit upperLimit];
end


