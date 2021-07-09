%% DistanceMatrix
% Author: Jan Glaubitz 
% Date: June 22, 2021 
%
% Compute matrix containing the distances between points 
%
%  INPUT:
%  X :          data points 
%  evalpoints : evaluation points 
%
%  OUTPUT:
%  DM : distance matrix

function DM = DistanceMatrix( X, evalpoints )

    [N, dim] = size(X);
    [K, dim] = size(evalpoints);
    DM = zeros(N,K);

    for n = 1:N 
        for k = 1:K
            DM(n,k) = norm( X(n,:) - evalpoints(k,:) ); 
        end
    end

end