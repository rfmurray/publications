function main2motion()

%
% MAIN2MOTION.M  Run a demonstration for the perception of glow,
%                using motion as the 3D depth cue
%
%  Part of the Supplementary Information for
%
%            PERCEIVED 3D SHAPE TOGGLES PERCEIVED GLOW
%       Minjung Kim, Laurie M. Wilcox, and Richard F. Murray
%
%

addpath( fullfile( pwd, 'tools' ) );

% Open Screen
InitializeMatlabOpenGL();
AssertOpenGL();
HideCursor();
ListenChar(2);
Screen('Preference','SyncTestSettings', .03);
win = Screen('OpenWindow', max(Screen('Screens')), [0 0 0], [], [], [], 0);

% get access to OpenGL constants
global GL

% null the colour lookup table
BackupCluts;
clut = repmat(linspace(0,1,256), [3 1])';
Screen('LoadNormalizedGammaTable',win,clut);

% Configure OpenGL
configOpenGL(win);

% Set some parameters
lengthM   = 0.10;         % length of side of square in metres
viewDistM = .74;          % viewing distance
thetaMaxDeg = 45;              % rotation amplitude
thetaSpeedDegPerS = 900/pi;    % rotation speed
[ screenWidthMM, screenHeightMM ] = Screen('DisplaySize',win);  % screen width and height (mm)

% Load the scene data.
tmp = load('scene.mat');

x =  tmp.x;
y =  tmp.y;
z =  tmp.z;
lum = tmp.lum;

[m,n] = size(x);

% Get the vertices for the dark-means-deep and bright-means-deep stimuli.
% For the bright-means-deep stimuli, set all z to -z.
v.left.data  = double(permute( cat( 3, x(:), y(:),  z(:) ), [ 3 1 2 ] ));
v.right.data = double(permute( cat( 3, x(:), y(:), -z(:) ), [ 3 1 2 ] ));
v.ind.data = mesh2ind(m,n); 

% Get the luminance map and prepare it for use in OpenGL.
tex = packageLum(win, lum);

% Display until key press.
timeOn = GetSecs;
while 1

    % calculate the current rotation angle
    eTime = GetSecs - timeOn;
    theta = thetaMaxDeg*cosd(thetaSpeedDegPerS*eTime);
    
    % get ready to do some OpenGL operations
    Screen('BeginOpenGL', win);
    glClear();
    
    % Set projection matrix
    glMatrixMode(GL.PROJECTION);
    glLoadIdentity;
%     gluPerspective(2*ATand2(.5*screenHeightM, viewDistM), screenWidthM/screenHeightM, .5*viewDistM, 1.5*viewDistM);
    gluPerspective(2*ATand2(.5*screenHeightMM/1000, viewDistM), screenWidthMM/screenHeightMM, .5*viewDistM, 1.5*viewDistM);
% *** changed a few variables
    
    % Set modelview matrix
    glMatrixMode(GL.MODELVIEW);
    glLoadIdentity;
    glTranslated(0,0,-viewDistM); % needed for glvStereoView
    
    % Put the luminance map onto the mesh.
    % The same luminance map is used for both stimuli.
    glBindTexture(GL.TEXTURE_2D, tex.lum.id);
    glTexCoordPointer(2,GL.DOUBLE,0,tex.ind.data);
    
    % Draw the dark-means-deep stimulus on the left
    glVertexPointer(size(v.left.data,1), GL.DOUBLE, 0, v.left.data);
    glPushMatrix();
    glTranslated(-.75*lengthM,0,0);
    glRotated(theta,0,1,0);
    glDrawElements(GL.TRIANGLES, numel(v.ind.data), GL.UNSIGNED_INT, v.ind.data);
    glPopMatrix();
    
    % Draw the bright-means-deep stimulus on the right
    glVertexPointer(size(v.right.data,1), GL.DOUBLE, 0, v.right.data);
    glPushMatrix();
    glTranslated( .75*lengthM,0,0);
    glRotated(theta,0,1,0);
    glDrawElements(GL.TRIANGLES, numel(v.ind.data), GL.UNSIGNED_INT, v.ind.data);
    glPopMatrix();
   
    % done with OpenGL operations
    Screen('EndOpenGL', win);
    
    % Show some helpful text
    if IsWindows
        Screen('TextFont',  win, 'Courier New');
    else
        Screen('TextFont',  win, 'Monaco');
    end
    Screen('TextColor', win, [255 255 255]);
    Screen('TextSize',  win, 16);
    Screen('DrawText',  win, 'Left:  dark-means-deep',   100, 100);
    Screen('DrawText',  win, 'Right: bright-means-deep', 100, 125);
    Screen('DrawText',  win, 'Click the mouse to exit.', 100, 175);
    Screen('Flip', win);
    
    % On mouse click, exit the loop.
    [~, ~, buttons] = GetMouse;
    if any(buttons)
        pause(.2);
        break
    end
    
end

% Close Screen
Screen('CloseAll');
RestoreCluts;
ListenChar(0);
ShowCursor();

end
