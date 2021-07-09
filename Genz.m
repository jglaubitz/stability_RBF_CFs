%% Genz
% Author: Jan Glaubitz 
% Date: June 22, 2021 
%
% Compute the the RBF-CF weights 
%
%  INPUT: 
%  a, b :   left and right boundary of the domain 
%  X :      data points 
%  w :      cubature weights 
%  CC :     number of experiments 
%  noise_level : amount of white uniform noise
%
%  OUTPUT:
%  error :  average errors for the Genz test functions 

%%
function error = Genz( a, b, X, w, CC, noise_level )
    
    [N,dim] = size(X); % number of data points 
    error = [0, 0, 0, 0]; % list for errors 
    
    for c=1:CC 

        %% set up Genz test functions 
        a = rand(1,dim); b = rand(1,dim); % random parameters 
        if dim == 1  
            genz1 = @(x) cos( 2*pi*b(1) + a(1)*x ); % osciallatory 
            genz2 = @(x) ( a(1)^(-2) + (x-b(1)).^2 ).^(-1); % product peak 
            genz3 = @(x) ( 1 + a(1)*x ).^(-2); % corner peak 
            genz4 = @(x) exp( -a(1)^2*(x-b(1)).^2 ); % Gaussian
            I1 = integral(genz1,0,1,'AbsTol',1e-14); 
            I2 = integral(genz2,0,1,'AbsTol',1e-14); 
            I3 = integral(genz3,0,1,'AbsTol',1e-14); 
            I4 = integral(genz4,0,1,'AbsTol',1e-14); 
        elseif dim == 2  
            genz1 = @(x,y) cos( 2*pi*b(1) + a(1)*x + a(2)*y ); % osciallatory 
            genz2 = @(x,y) ( a(1)^(-2) + (x-b(1)).^2 ).^(-1).*( a(2)^(-2) + (y-b(2)).^2 ).^(-1); % product peak 
            genz3 = @(x,y) ( 1 + a(1)*x + a(2)*y ).^(-3); % corner peak 
            genz4 = @(x,y) exp( -a(1)^2*(x-b(1)).^2 - a(2)^2*(y-b(2)).^2 ); % Gaussian 
            I1 = integral2(genz1,0,1,0,1,'AbsTol',1e-14); 
            I2 = integral2(genz2,0,1,0,1,'AbsTol',1e-14); 
            I3 = integral2(genz3,0,1,0,1,'AbsTol',1e-14); 
            I4 = integral2(genz4,0,1,0,1,'AbsTol',1e-14); 
        elseif dim == 3  
            genz1 = @(x,y,z) cos( 2*pi*b(1) + a(1)*x + a(2)*y + a(3)*z ); % osciallatory 
            genz2 = @(x,y,z) ( a(1)^(-2) + (x-b(1)).^2 ).^(-1).*( a(2)^(-2) + (y-b(2)).^2 ).^(-1).*( a(3)^(-2) + (z-b(3)).^2 ).^(-1); % product peak 
            genz3 = @(x,y,z) ( 1 + a(1)*x + a(2)*y + a(3)*z).^(-4); % corner peak 
            genz4 = @(x,y,z) exp( -a(1)^2*(x-b(1)).^2 - a(2)^2*(y-b(2)).^2 - a(3)^2*(z-b(3)).^2); % Gaussian 
            I1 = integral3(genz1,0,1,0,1,0,1,'AbsTol',1e-14); 
            I2 = integral3(genz2,0,1,0,1,0,1,'AbsTol',1e-14); 
            I3 = integral3(genz3,0,1,0,1,0,1,'AbsTol',1e-14); 
            I4 = integral3(genz4,0,1,0,1,0,1,'AbsTol',1e-14); 
        else 
            error('Desired dimension not yet implemented!') 
        end

        %% function values 
        if dim == 1  
            genz1_values = genz1(X);
            genz2_values = genz2(X);
            genz3_values = genz3(X);
            genz4_values = genz4(X);
        elseif dim == 2  
            genz1_values = genz1(X(:,1),X(:,2)); 
            genz2_values = genz2(X(:,1),X(:,2));
            genz3_values = genz3(X(:,1),X(:,2));
            genz4_values = genz4(X(:,1),X(:,2));
        elseif dim == 3  
            genz1_values = genz1(X(:,1),X(:,2),X(:,3));
            genz2_values = genz2(X(:,1),X(:,2),X(:,3));
            genz3_values = genz3(X(:,1),X(:,2),X(:,3));
            genz4_values = genz4(X(:,1),X(:,2),X(:,3));
        else 
            error('Desired dimension not yet implemented!') 
        end
        
        %% noise 
        if noise_level ~= 0 
            noise = 10^(-noise_level)*(2*rand(N,1)-1); % generate uniform noise 
        else 
            noise = 0; % no noise 
        end
            
        %% sum errors       
        error(1) = error(1) + abs(I1-dot(w,genz1_values + noise)); % absolute error
        error(2) = error(2) + abs(I2-dot(w,genz2_values + noise)); % absolute error
        error(3) = error(3) + abs(I3-dot(w,genz3_values + noise)); % absolute error
        error(4) = error(4) + abs(I4-dot(w,genz4_values + noise)); % absolute error
        
        
    end
    
    % compute average errors 
    error = error/CC;
    
end