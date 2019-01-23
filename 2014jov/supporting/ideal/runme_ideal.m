% runme_ideal.m  Simulate ideal observer in lighting direction discrimination
%                task from Morgenstern, Geisler, and Murray (2014)

clear; clc;

% set mean lighting direction in polar coordinates
mtheta = 60;   % slant
mphi = 0;      % tilt

% reflectance distribution
rmu = 0.6;     % reflectance mean
rsigma = 0.2;  % reflectance nominal standard deviation
rclip = 2.0;   % number of nominal standard deviations from mean where reflectance distribution is clipped

% viewing direction
vdir = [ 0 0 1 ];

% ambient-to-total illuminance ratio, i.e., Ea/(Ep+Ea)
atratio = 0.20;

% load normal vectors for faces of sphere-like polyhedron
nvec = load('spherenormals.txt','-ascii');

% use normals only for visible faces
% (this gives the same faces used in the experiment with human observers)
nvec = nvec((nvec*vdir')>0,:);

% set lighting direction perturbation
% i.e., observer discriminates between lighting directions
%       ( mtheta+dtheta, mphi+dphi ) and ( mtheta-dtheta, mphi-dphi )
dtheta = 0;  % slant perturbation
dphi = 1;    % tilt perturbation

% calculate lighting direction vectors
lvec1 = [ sind(mtheta+dtheta)*sind(mphi+dphi) sind(mtheta+dtheta)*cosd(mphi+dphi) cosd(mtheta+dtheta) ];
lvec2 = [ sind(mtheta-dtheta)*sind(mphi-dphi) sind(mtheta-dtheta)*cosd(mphi-dphi) cosd(mtheta-dtheta) ];

% find ideal observer's performance
[ dprime, hit, fa ] = idealdprime( rmu, rsigma, rclip, lvec1, lvec2, atratio, nvec )
