function glvStereoView( righteye, screendist, screensize, iodist, neardist, fardist )

% glvStereoView  Set a stereo view projection matrix
% 
%     usage:  glvStereoView( righteye, screendist, screensize, iodist, neardist, fardist )
% 
%     input arguments
%         righteye   -- flag (0 or 1) that indicates whether the projection matrix is for the left eye or right eye
%         screendist -- distance from the cyclopean eye to the screen
%         screensize -- 1 x 2 matrix [ h v ] that gives the horizontal and vertical size of the screen
%         iodist     -- interocular distance
%         neardist   -- distance to the near clipping plane
%         fardist    -- distance to the far clipping plane

% find eye displacement
if righteye
    deltae = iodist/2;
else
    deltae = -iodist/2;
end

% find x coordinates of near end of frustum, in a coordinate system
% where the eye is at (deltae,0,0)
xleft  = deltae + (neardist/screendist)*( (-screensize(1)/2)-deltae );
xright = deltae + (neardist/screendist)*( ( screensize(1)/2)-deltae );

% find y coordinates of near end of frustum, in a coordinate system
% where the eye is at (deltae,0,0)
ytop    = (neardist/screendist)*( screensize(2)/2);
ybottom = (neardist/screendist)*(-screensize(2)/2);

% set projection matrix, in a coordinate system where the eye is at (0,0,0)
glFrustum( xleft-deltae, xright-deltae, ybottom, ytop, neardist, fardist );

% shift coordinate system so that the eye is effectively at (deltae,0,0)
glTranslated( -deltae, 0, 0 );

end
