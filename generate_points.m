%% generate_points
% Author: Jan Glaubitz 
% Date: June 22, 2021 
%
% Generates data points in [a,b]^dim
%
%  INPUT:
%  dim :    dimension of the domain
%  a, b :   left and right boundary of domain \Omega = [a,b]^dim
%  N :      number of data points  
%  points : type of data points 
%
%  OUTPUT:
%  X : data points 
%%
function X = generate_points( dim, a, b, N, points )
   
    X = zeros(N,dim); % array of data points 

    %% equidistant data points 
    if strcmp( points, 'equid')  
    
        n = ceil( N^(1/dim) ); % number of points in every direction 
        coord_aux = linspace(0, 1, n); 

        % different dimensions 
        if dim == 1 
            X = coord_aux';
        elseif dim == 2 
            [PointsX, PointsY] = meshgrid(coord_aux, coord_aux); 
            X = [ reshape(PointsX, numel(PointsX),1)'; 
                  reshape(PointsY, numel(PointsY),1)']';
        else 
            error('Desired dimension not yet implemented!') 
        end    
        
    %% Halton points 
    elseif strcmp(points,'Halton')
        p = haltonset(dim); % generate Halton point set in [0,1]^dim
        p = scramble(p,'RR2'); % scramble point set 
        X = net(p,N); % generate the first N points 
    
    %% random points
    elseif strcmp(points,'random')
        rng('default'); % for reproducibility
        rng(1,'twister'); % for reproducibility
        X = rand(N,dim); % generate random points
    
    %% else
    else
        error('Desired points not yet implemented!')    
    end
    
    % if in one dimension, order points
    if dim==1 
        X = sort(X); % order points
    end
    % scale points to [a,b]^dim
    X = (b-a)*X + a; % scale 
    
end

