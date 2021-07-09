%% RBF_moments
% Author: Jan Glaubitz 
% Date: June 22, 2021 
%
% Compute the RBF's moments
%
%  INPUT: 
%  a, b :   left and right boundary of the domain 
%  kernel : kernel 
%  rbf :    radial basis function 
%  ep :     shape parameter 
%  X :      data points 
%
%  OUTPUT:
%  m_RBF :  moments of the translated RBFs

%%
function m_RBF = RBF_moments( a, b, kernel, rbf, ep, X )
    
    [N,dim] = size(X); % number of data points 
    m_RBF = zeros(N,1);
    
    %% One dimensional
    if dim==1 
        
        %% Gaussian 
        if strcmp(kernel,'G') 
            m_RBF = (0.5/ep)*sqrt(pi)*( erf( ep*(b-X) ) - erf( ep*(a-X) ) ); % moments 
        elseif strcmp(kernel,'TPS') 
            k = 2; 
            m_RBF = (X-a).^(k+1).*( log(X-a+10^(-14))/(k+1) - 1/(k+1)^2 ) + ... 
                (b-X).^(k+1).*( log(b-X+10^(-14))/(k+1) - 1/(k+1)^2 ); % moments
        elseif strcmp(kernel,'cubic')
            k = 3; 
            m_RBF = ( (a-X).^(k+1) + (b-X).^(k+1) )/(k+1); % moments
        elseif strcmp(kernel,'quintic')
            k = 3; 
            m_RBF = ( (a-X).^(k+1) + (b-X).^(k+1) )/(k+1); % moments 
        elseif strcmp(kernel,'Wendland')  
            for n=1:N 
               % Support of kernel with center x_n is [c,d] with 
               c = max(a,X(n)-1/ep(n)); d = min(X(n)+1/ep(n),b); 
               int = @(x) rbf(ep(n),abs( x-X(n) )); % integrand 
               m_RBF(n) = integral( @(x) int(x), c, d ); % nth moment
            end
        else    
            rbf_basis = @(x) rbf(ep,abs(x-X')); % RBF basis
            syms x 
            m_RBF = integral( @(x) rbf_basis(x), a, b, 'ArrayValued', true )'; % moments 
        end
    
    %% Two dimensional
    elseif dim==2
        if strcmp(kernel,'G') 
            for n=1:N 
                mx = (0.5/ep)*sqrt(pi)*( erf( ep*(b-X(n,1)) ) - erf( ep*(a-X(n,1)) ) ); % component in x direction
                my = (0.5/ep)*sqrt(pi)*( erf( ep*(b-X(n,2)) ) - erf( ep*(a-X(n,2)) ) ); % component in x direction
                m_RBF(n) = mx*my; % moments 
            end
        elseif strcmp(kernel,'TPS') 
            I_tr = @(u,v) (u/144)*( ... 
                    24*u^3*atan(v/u) + 6*v*(3*u^2+v^2)*log(u^2+v^2) - 33*u^2*v - 7*v^3 ...
                ); % reference integral 
            % compute moments 
            for n=1:N 
                % shifted edges of the rectangle 
                c = a; d = b; % we assume the domain [a,b]^2
                a_tilde = abs( a - X(n,1) ); 
                b_tilde = abs( b - X(n,1) ); 
                c_tilde = abs( c - X(n,2) ); 
                d_tilde = abs( d - X(n,2) ); 
                
                % partition rectangle in 8 right triangles and compute the
                % corresponding integrals 
                I(1) = I_tr(b_tilde,d_tilde); 
                I(2) = I_tr(d_tilde,b_tilde); 
                I(3) = I_tr(d_tilde,a_tilde); 
                I(4) = I_tr(a_tilde,d_tilde); 
                I(5) = I_tr(a_tilde,c_tilde); 
                I(6) = I_tr(c_tilde,a_tilde); 
                I(7) = I_tr(c_tilde,b_tilde); 
                I(8) = I_tr(b_tilde,c_tilde); 
                I(isnan(I))=0; % set all NaN values to zero;
                % sum these up to get the moment 
                m_RBF(n) = ( b_tilde*d_tilde ~= 0 )*(I(1)+I(2)) + ... 
                    ( a_tilde*d_tilde ~= 0 )*(I(3)+I(4)) + ...
                    ( a_tilde*c_tilde ~= 0 )*(I(5)+I(6)) + ...
                    ( b_tilde*c_tilde ~= 0 )*(I(7)+I(8));
            end
        elseif strcmp(kernel,'cubic')
            I_tr = @(u,v) (u/40)*( ... 
                    3*u^4*asinh(v/u) + ... 
                    v*( 5*u^2 + 2*v^2 )*sqrt( u^2 + v^2 ) ...
                ); % reference integral
            % compute moments 
            for n=1:N 
                % shifted edges of the rectangle 
                c = a; d = b; % we assume the domain [a,b]^2
                a_tilde = abs( a - X(n,1) ); 
                b_tilde = abs( b - X(n,1) ); 
                c_tilde = abs( c - X(n,2) ); 
                d_tilde = abs( d - X(n,2) );
                % partition rectangle in 8 right triangles and compute the
                % corresponding integrals 
                I(1) = I_tr(b_tilde,d_tilde); 
                I(2) = I_tr(d_tilde,b_tilde); 
                I(3) = I_tr(d_tilde,a_tilde); 
                I(4) = I_tr(a_tilde,d_tilde); 
                I(5) = I_tr(a_tilde,c_tilde); 
                I(6) = I_tr(c_tilde,a_tilde); 
                I(7) = I_tr(c_tilde,b_tilde); 
                I(8) = I_tr(b_tilde,c_tilde); 
                I(isnan(I))=0; % set all NaN values to zero;
                % sum these up to get the moment 
                m_RBF(n) = ( b_tilde*d_tilde ~= 0 )*(I(1)+I(2)) + ... 
                    ( a_tilde*d_tilde ~= 0 )*(I(3)+I(4)) + ...
                    ( a_tilde*c_tilde ~= 0 )*(I(5)+I(6)) + ...
                    ( b_tilde*c_tilde ~= 0 )*(I(7)+I(8));
            end
        elseif strcmp(kernel,'quintic')
            I_tr = @(u,v) (u/336)*( ... 
                    15*u^6*asinh(v/u) + ... 
                    v*( 33*u^4 + 26*u^2*v^2 + 8*v^4 )*sqrt( u^2 + v^2 ) ...
                ); % reference integral
            % compute moments 
            for n=1:N 
                % shifted edges of the rectangle 
                c = a; d = b; % we assume the domain [a,b]^2
                a_tilde = abs( a - X(n,1) ); 
                b_tilde = abs( b - X(n,1) ); 
                c_tilde = abs( c - X(n,2) ); 
                d_tilde = abs( d - X(n,2) );
                % partition rectangle in 8 right triangles and compute the
                % corresponding integrals 
                I(1) = I_tr(b_tilde,d_tilde); 
                I(2) = I_tr(d_tilde,b_tilde); 
                I(3) = I_tr(d_tilde,a_tilde); 
                I(4) = I_tr(a_tilde,d_tilde); 
                I(5) = I_tr(a_tilde,c_tilde); 
                I(6) = I_tr(c_tilde,a_tilde); 
                I(7) = I_tr(c_tilde,b_tilde); 
                I(8) = I_tr(b_tilde,c_tilde); 
                I(isnan(I))=0; % set all NaN values to zero;
                % sum these up to get the moment 
                m_RBF(n) = ( b_tilde*d_tilde ~= 0 )*(I(1)+I(2)) + ... 
                    ( a_tilde*d_tilde ~= 0 )*(I(3)+I(4)) + ...
                    ( a_tilde*c_tilde ~= 0 )*(I(5)+I(6)) + ...
                    ( b_tilde*c_tilde ~= 0 )*(I(7)+I(8));
            end
        elseif strcmp(kernel,'septic')
            I_tr = @(u,v) (u/3346)*( ... 
                    105*u^8*asinh(v/u) + ... 
                    v*( 279*u^6 + 326*u^4*v^2 + 200*u^2*v^4 + 48*v^6 )*sqrt( u^2 + v^2 ) ...
                ); % reference integral
            % compute moments 
            for n=1:N 
                % shifted edges of the rectangle 
                c = a; d = b; % we assume the domain [a,b]^2
                a_tilde = abs( a - X(n,1) ); 
                b_tilde = abs( b - X(n,1) ); 
                c_tilde = abs( c - X(n,2) ); 
                d_tilde = abs( d - X(n,2) );
                % partition rectangle in 8 right triangles and compute the
                % corresponding integrals 
                I(1) = I_tr(b_tilde,d_tilde); 
                I(2) = I_tr(d_tilde,b_tilde); 
                I(3) = I_tr(d_tilde,a_tilde); 
                I(4) = I_tr(a_tilde,d_tilde); 
                I(5) = I_tr(a_tilde,c_tilde); 
                I(6) = I_tr(c_tilde,a_tilde); 
                I(7) = I_tr(c_tilde,b_tilde); 
                I(8) = I_tr(b_tilde,c_tilde); 
                I(isnan(I))=0; % set all NaN values to zero;
                % sum these up to get the moment 
                m_RBF(n) = ( b_tilde*d_tilde ~= 0 )*(I(1)+I(2)) + ... 
                    ( a_tilde*d_tilde ~= 0 )*(I(3)+I(4)) + ...
                    ( a_tilde*c_tilde ~= 0 )*(I(5)+I(6)) + ...
                    ( b_tilde*c_tilde ~= 0 )*(I(7)+I(8));
            end
        else   
            for n=1:N 
               int = @(x,y) rbf( ep, sqrt( (x-X(n,1)).^2 + (y-X(n,2)).^2 ) ); % integrand  
               m_RBF(n) = integral2( int, a,b, a,b ); % moment
            end
        end
    
    %% Higher dimensional 
    else
        error('Desried dimension not yet implemented')
        
 end