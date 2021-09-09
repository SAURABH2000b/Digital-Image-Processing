the image_filtering.m file applies different 2D convolution filters on the input image. the default filter width is set to 5. 
each of these filters are derived from six different filtering methods. and the sampling rate (for certain methods only) is chosen to be 1/3. 
the filtering methods implemented are:
box filter (natural sampling rate)
tent filter (sampling rate embedded in equation to make area under tent equal to 1)
gaussian filter
B-Spline cubic filter
Catmull-Rom cubic filter
Mitchell-Netravali cubic filter

application of filtering: to reduce aliasing or jagginess in the input image. one drawback of antialiasing filtering is increased blurriness in the image.
the level of blurriness depends on the antialiasing filter used. worst filter is box filter, best is Mitchell Netravali filter.
