

This package contains supporting data and MATLAB code for:

     Morgenstern, Y., Geisler, W. S., & Murray, R. F. (2014).  Human vision is attuned to the diffuseness of natural light.  Journal of Vision, 14(9):15, 1-18.

This package can be downloaded from http://purl.org/NET/rfm/jov2014 .



The directory 'lightprobes' contains data and code relating to the light probes that Morgenstern et al. measured using a multidirectional photometer.

- light_probes.txt contains the complex-valued coefficients of a spherical harmonic representation of the light probes.  Comments in the header of this file explain how the light probes are represented.

- runme_lightprobe.m illustrates how to load and use the data in light_probes.txt.  It loads the spherical harmonic coefficients of a randomly chosen light probe, reconstructs the light probe, shows an image of it, and calculates the light probe's illuminance contrast energy (ICE).

- sphharm.m generates spherical harmonics.  It is used by runme_lightprobe.m.

- ice.m calculates illuminance contrast energy (ICE).  It is used by runme_lightprobe.m.



The directory 'ideal' contains code that implements the ideal observer used in Morgenstern et al.'s experiment 2.

- runme_ideal.m calculates ideal performance at discriminating between two images of a polyhedral approximation to a sphere, illuminated from two different lighting directions.

- idealdprime.m does the actual calculation of ideal performance.  It is used by runme_ideal.m.

- spherenormals.txt contains the surface normals of a polyhedral approximation to a sphere.  It is used by runme_ideal.m.



This code was developed in MATLAB 8.3 (R2014a).  It does not require add-on toolboxes.
