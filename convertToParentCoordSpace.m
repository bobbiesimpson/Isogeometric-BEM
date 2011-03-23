function [ xi ] = convertToParentCoordSpace( xi_param, range )
% Take in a coordinate in the paramater space and convert it to the
% parental coordinate space

xi= ( 2 * xi_param - (range(2) + range(1) ) ) / ( range(2) - range(1) );

end

