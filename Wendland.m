% Wendland - a small script to iteratively compute Wendland's compact
% support radial basis function (RBF).
%
% References:
% Wendland, H.  (1995)  Piecewise polynomial, positive definite and compactly
%                       supported radial functions of minimal degree. 
%                       Adv. Comput. Math., 4: 389–396.
% Wendland, H.  (2004)  Scattered Data Approximation.
%                       Chapter 9.
% Fasshauer, G. (2007)  Meshfree Approximation Methods With MATLAB.
%                       Chapter 11. 
%
% Inputs:
% d = dimension of state space (positive integer)
% k = smoothness parameter (non-negative integer)
%
% Outputs:
% out = function handle to the Wendland RBF
%       accepts values in [0,inf)
%       normalised such that out(0) = 1
%
% Background:
% Wendland(d,k) is a polynomial on [0,1] and vanishes on (1,inf)
% Wendland(d,k) is C^2k(R); i.e. it has 2k continuous derivatives
% As a RBF, Wendland(d,k) is positive definite on R^d
%
% Example:
% Wendland(d,0) = (1-r)_+^l
% Wendland(d,1) = (1-r)_+^(l+1) * ((l+1)r + 1)
% Wendland(d,2) = (1-r)_+^(l+2) * ((l^2+4l+3)r^2 + (3l+6)r + 3)
%                 in each case up to a multiplicative constant
%                 z_+^l := max(0,z)^l
%                 l = floor(d/2) + k + 1
%
% © Chris Oates 2017.

function out = Wendland(d,k)
l = floor(d/2) + k + 1;
phi = @(r) (1-r)^l; % power function (Sec 11.2 in Fasshauer)
for i = 1:k
    phi = I(phi); % iterative integration (Defn 11.2 in Fasshauer)
end
nor = phi(0);
syms r
out = phi(r) / nor; % normalise at r = 0
out = heaviside(r) * heaviside(1-r) * out; % truncation to compact support
out = matlabFunction(out);
end

% Integration operator
% f(r) -> \int_r^\infty t f(t) dt
function out = I(in)
syms r
integrand = @(t) t * in(t);
out = int(integrand,r,1);
out = matlabFunction(out);
end





