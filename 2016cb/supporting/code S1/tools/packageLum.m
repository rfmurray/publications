function tex = packageLum(win, lum)

% packageLum Package the pre-rendered luminance map to be used for OpenGL.
% 
%     usage:  packageLum( lum )
% 
%     input arguments
%         
%         lum -- luminance map. Must be m x n, where m and n are powers
%                of two.
%
%     output variables
%         win -- window to use OpenGL in
%         tex -- a struct with all texture-related information, 
%                including the luminance map and pointer to the luminance
%                map

% Correct for monitor gamma.
gamma = [60 10 2 0]; % Parameters of the gamma function
lum2rgb = @(lum, gamma) round((255-gamma(2))*power((lum-gamma(4))/gamma(1),1/gamma(3))+gamma(2)); % inverse gamma function

maxLumDisp = gamma(1)+gamma(4);
minLumDisp = gamma(4);

lum = max(maxLumDisp * lum./max(lum(:)), minLumDisp); % Scale to maximum value of display
lum = lum2rgb(lum, gamma);                            % Gamma correction

[m,n] = size(lum);

% Reorient luminance matrix to OpenGL's expectations.
lum = repmat(lum, [1 1 3]);
lum = permute(lum, [3 2 1]);
tex.lum.data = uint8(lum);

% Make luminance look-up table. This will be used to texture-map the
% luminance map onto the mesh.

tx = repmat( linspace(0,1,n),[ m 1 ]);
ty = repmat( linspace(0,1,m)',[ 1 n ]);
txy = cat(3,tx,ty);
tex.ind.data = double(permute( txy, [ 3 1 2 ] ));

% Bind luminance map
global GL;
Screen('BeginOpenGL', win); 

tex.lum.id = glGenTextures(1);
glBindTexture( GL.TEXTURE_2D, tex.lum.id );
glTexImage2D(GL.TEXTURE_2D,0,GL.RGB,size(tex.lum.data,2),size(tex.lum.data,3),0,GL.RGB,GL.UNSIGNED_BYTE,tex.lum.data);

glTexParameterfv(GL.TEXTURE_2D,GL.TEXTURE_WRAP_S,GL.REPEAT);
glTexParameterfv(GL.TEXTURE_2D,GL.TEXTURE_WRAP_T,GL.REPEAT);
glTexParameterfv(GL.TEXTURE_2D,GL.TEXTURE_MAG_FILTER,GL.NEAREST);
glTexParameterfv(GL.TEXTURE_2D,GL.TEXTURE_MIN_FILTER,GL.NEAREST);
glTexEnvfv(GL.TEXTURE_ENV,GL.TEXTURE_ENV_MODE,GL.REPLACE);

Screen('EndOpenGL', win); 
end
