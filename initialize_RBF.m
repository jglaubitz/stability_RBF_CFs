%% initialize_RBF
% Author: Jan Glaubitz 
% Date: June 22, 2021 
%
% Define the radial basis function 
%
%  INPUT:
%  kernel : kernel function 
%  dim : dimension of the domain 
%  order : smoothness paramter in case Wendland's compactly supported RBFs
%  are used 
%
%  OUTPUT:
%  rbf : radial basis function 

%%
function rbf = initialize_RBF( kernel, dim, order )
    
    %% Gauss, MQ, and IQ
    rbf_G = @(ep,r) exp(-(ep*r).^2); % Gaussians
    rbf_MQ = @(ep,r) sqrt(1 + (ep*r).^2); % multiquadrics
    rbf_IQ = @(ep,r) 1./(1 + (ep*r).^2); % inverse quadrics 
    
    %% Wendland 
    Wendland_fun = Wendland(dim,order); % Wendland function in dimension "dim and" smoothness order "order"
    rbf_W = @(ep,r) Wendland_fun(ep.*r); % RBF 
    
    %% PHS 
    rbf_TPS = @(ep,r) r.*log( r.^r ); % thin plate spline (TPS)
    rbf_cubic = @(ep,r) (ep*r).^3; % cubic 
    rbf_quintic = @(ep,r) (ep*r).^5; % quintic 
    
    % Choose the RBF
    if strcmp(kernel,'G')
        rbf = rbf_G; 
    elseif strcmp(kernel,'MQ')
        rbf = rbf_MQ; 
    elseif strcmp(kernel,'IQ')
        rbf = rbf_IQ; 
    elseif strcmp(kernel,'Wendland')
        rbf = rbf_W; 
    elseif strcmp(kernel,'TPS') 
        rbf = rbf_TPS;
    elseif strcmp(kernel,'cubic')
        rbf = rbf_cubic;
    elseif strcmp(kernel,'quintic')
        rbf = rbf_quintic; 
    else 
        error('Desired kernel not yet implemented!') 
    end 

end