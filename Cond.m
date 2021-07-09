%% Cond
% Author: Jan Glaubitz 
% Date: June 22, 2021 
%
% Compute the condition number of the original Vandermonde matrix 
% That is, w.r.t the 'kernel basis' of the space S_{N,d}
%
%  INPUT: 
%  a, b :   left and right boundary of the domain [a,b]^dim
%  rbf :    radial basis function 
%  ep :     shape parameter 
%  X :      data points 
%  d : polynomial degree 
%
%  OUTPUT:
%  cond_number : condition number 

%%
function cond_number = Cond( a, b, rbf, ep, X, d )
    
    [N,dim] = size(X); % number of data points and dimension
    DM = DistanceMatrix(X,X); % matrix with distances between points 
    V_rbf = rbf(ep',DM); % Vandermonde matrix
    
    %% no polynomials included
    if d < 0 
        
        % Vandermonde matrix 
        A = V_rbf;
        
    %% polynomials included     
    else
        
        if dim==1 
        
            % define the polynomial basis p_0,...,p_d 
            if d == 0
                p = @(x) x.^0; 
            elseif d == 1
                p = @(x) [x.^0 , x]; 
            elseif d == 2
                p = @(x) [x.^0 , x , x.^2]; 
            elseif d == 3
                p = @(x) [x.^0 , x , x.^2 , x.^3]; 
            else 
                'd to large';
            end
            % Polynomial matrix P 
            P = p(X); % polynomial matrix
            % Vandermonde Matrix A = ( V P ; P^T 0 )
            A = [ V_rbf P; P' zeros(d+1,d+1) ];
        
        elseif dim==2 
            
            K = nchoosek(dim+d,dim); % number of basis functions 
            % define the polynomial basis p_1,...,p_K 
            if d == 0
                p = @(x,y) x.^0 + y.^0; 
            elseif d == 1
                p = @(x,y) [x.^0 , x, y]; 
            elseif d == 2
                p = @(x,y) [x.^0 , x , y, x.^2, x.*y, y.^2]; 
            elseif d == 3
                p = @(x,y) [x.^0 , x , y, x.^2, x.*y, y.^2, ... 
                    x.^3, (x.^2).*(y.^1), (x.^1).*(y.^2), y.^3]; 
            elseif d == 4
                p = @(x,y) [x.^0 , x , y, x.^2, x.*y, y.^2, ... 
                    x.^3, (x.^2).*(y.^1), (x.^1).*(y.^2), y.^3, ... 
                    x.^4, (x.^3).*(y.^1), (x.^2).*(y.^2), (x.^1).*(y.^3), y.^4]; 
            else 
                'd to large';
            end
            % Polynomial matrix P 
            P = p( X(:,1), X(:,2) ); % polynomial matrix
            % Vandermonde Matrix A = ( V P ; P^T 0 )
            A = [ V_rbf P; P' zeros(K,K) ];
          
        end
            
    end

    %% Condition number 
    cond_number = cond(A);
    
end