-------------------------------------------------------------------

             PERCEIVED 3D SHAPE TOGGLES PERCEIVED GLOW
        Minjung Kim, Laurie M. Wilcox, and Richard F. Murray

                         SUPPLEMENTAL CODE

       For any questions, please contact minjung.kim@nyu.edu.

-------------------------------------------------------------------

Our code comes in two parts: (1) rendering the luminance map L[i,j], and (2) putting L[i,j] onto shape meshes, Z[i,j] and -Z[i,j].


(1). Rendering the luminance map.
We used Radiance, an open source 3D graphics package that closely mimics the physics of light transport, to render the luminance images. Please install Radiance from 
     http://www.radiance-online.org/ 

To run our code, please call     >> main1For details of the rendering parameters, please see:
     http://radsite.lbl.gov/radiance/man_html/rpict.1.html
     http://radsite.lbl.gov/radiance/man_html/pfilt.1.html

This code creates a scene.mat, a Matlab file containing the mesh structure (variables x, y, z), and the luminance map (variable lum). 


(2). Putting the luminance map onto the meshes.
We used the Psychophysics Toolbox’s interface to OpenGL to put the luminance maps onto the depth maps, then to display them either in stereo or in motion. Please install the Psychophysics Toolbox from 
     http://psychtoolbox.org/ 

To run our code, please call     >> main2stereo % stereo demo     >> main2motion % motion demoThough the actual experiment was run on a mirror stereoscope, main2stereo uses temporally interleaved stereo, to be viewed with LCD shutter glasses, for demonstration purposes. 
