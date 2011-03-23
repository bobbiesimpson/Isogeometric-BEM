function [ xi_param ] = convertToParamSpace( xi, range)
% Take in a coordinate in the parent coordinate system xi in [-1,1]
% and convert it to the parameter space

xi_param=( (range(2) - range(1)) * xi  + (range(2) + range(1) ) ) / 2;


end

