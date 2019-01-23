%
% main1.m  Generate a diffusely lit stimulus for the glow experiment.
% 
%  Part of the Supplementary Information for
%
%            PERCEIVED 3D SHAPE TOGGLES PERCEIVED GLOW
%       Minjung Kim, Laurie M. Wilcox, and Richard F. Murray
%

% *** Be sure to set the path to your installation of RADIANCE
% *** on line 25 (variable rayPath), below.

clear; clc; close all;

%% 1.  Setup

% Rendering can take more than an hour.  Set fast = 1 for a faster,
% low-resolution, low-quality rendering, if you just want to see whether
% the code works on your computer.  With fast = 1, the time estimates
% printed to the command window are invalid, and the whole script should
% finish in a few minutes.
fast = 1;

% add RADIANCE folders to system paths
rayPath = '/usr/local/ray';  % typical Radiance path
if isempty( strfind( getenv('PATH'), rayPath ) )
    setenv( 'PATH',    [ fullfile(rayPath,'bin') ':' getenv('PATH') ] );
    setenv( 'RAYPATH', [ '.:' fullfile(rayPath,'lib') ':' fullfile(pwd,'tools') ] );
end

% add tools folder to MATLAB path
addpath( fullfile( pwd, 'tools' ) );

%% 2.  Create the scene description for RADIANCE

% Set some parameters.
lengthM = 0.10;      % Size of square (meters)
sigmaM = 0.0068;     % Standard deviation of shape noise amplitude (meters), i.e., waviness of the surface
cutoff = 9;          % Upper cutoff frequency for shape noise (cycles/shape)
reflectance = .3;    % Surface reflectance
if fast
    nGrid = 64;      % Size of coordinate matrices (nGrid x nGrid)
else
    nGrid = 512;
end

% Make the coordinate matrices for the mesh.
xy = linspace(-lengthM/2,lengthM/2,nGrid);  % Coordinate vector
[x,y] = meshgrid(xy,-xy);                   % x and y coordinate matrices
z = sigmaM*bpnoise2d(nGrid,0,cutoff);       % z coordinate matrix (depth map), a sample of bandpass noise

% Write the mesh as a .rad file (RADIANCE text file).
fprintf('Writing the square mesh as a .rad file (takes about a minute) ...');
rad_surf('tmp/square',x,y,z,reflectance);

% Combine .rad files into an .oct file (RADIANCE binary).
%     square.rad       -- the blobby square that we just generated
%     light.rad        -- completely diffuse lighting
fprintf('\nConverting .rad files to an .oct file (takes a few seconds) ... ');
unix('oconv tmp/square.rad tools/light.rad > tmp/scene.oct');

%% 3. Render the scene in RADIANCE

% Set some parameters.
if fast
    renderRes = 128;                % Size of rendered image (renderRes x renderRes)
    ambAcc = .15;                   % RADIANCE quality setting; lower is better, 0 = best
else
    renderRes = 1024;
    ambAcc = 0;
end
finalRes  = floor(renderRes/2);     % Size of anti-aliased image (finalRes x finalRes)

% Render the scene, to produce scene_highres.hdr (an image file)
fprintf('\nRendering (May take more than an hour!  Check scene.log for progress.) ... ');
unix(sprintf('rpict -vtl -x %d -y %d -vv %.4f -vh %.4f -vu 0 1 0 -vp 0 0 200 -vd 0 0 -1 -ab 1 -aa %.4f -t 1 -e scene.log tmp/scene.oct > tmp/scene_highres.hdr', ...
    renderRes, renderRes, lengthM, lengthM, ambAcc));

% Anti-alias and downsample scene_highres.hdr, to produce scene.hdr (another image file)
unix(sprintf('pfilt -1 -x %d -y %d -e 1 -r .6 tmp/scene_highres.hdr > tmp/scene.hdr',...
    finalRes, finalRes));

fprintf('\nRendering done!\n');

%% 4.  Save coordinate matrices and rendered image in a .mat file

% save as .mat file
lum = hdrread('tmp/scene.hdr');
lum = double(lum(:,:,1));
save('scene.mat','x','y','z','reflectance','lum');

% show result
imshow( lum, [ 0 max(lum(:)) ] ); colormap gray; title 'rendered image'
