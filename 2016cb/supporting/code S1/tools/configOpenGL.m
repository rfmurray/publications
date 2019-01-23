function configOpenGL(win)

% configOpenGL  Configure OpenGL
% 
%     usage:  configOpenGL( win, tex )
% 
%     input arguments
%         win -- window to configure

global GL;

% Configuration for depth, texture.
Screen('BeginOpenGL',win);

% configure depth, etc.
glEnable(GL.DEPTH_TEST);
glShadeModel(GL.SMOOTH);
glEnableClientState(GL.VERTEX_ARRAY);

% configure for texture
glEnable(GL.TEXTURE_2D);
glEnableClientState(GL.TEXTURE_COORD_ARRAY);

Screen('EndOpenGL',win);

end