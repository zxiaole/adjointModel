function destinationData = fillInData(destinationData,originData)
% destinationData = fillInData(destinationData,originData)
%   Xiaole Zhang, Nov. 15, 2021
% fill in the valid data into the final data
% destinationData and originData shoud have the same shape
id = originData > 0;
destinationData(id) = originData(id);
end

