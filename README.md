# MPLocalizer

This repository contains stimulus code that was used to map the M and P subdivisions of human LGN using fMRI in:

Denison, R. N., Vu, A. T., Yacoub, E., Feinberg, D. A., & Silver, M. A. (2014). 
Functional mapping of the magnocellular and parvocellular subdivisions of human LGN. 
NeuroImage, 102 Pt 2, 358–369. http://doi.org/10.1016/j.neuroimage.2014.07.019

Please cite our paper if you use this code.
***
The main function is **runMPLocalizer.m**

The isoluminance point is set in **mpLocalizerColorParamsStim.m**

To display alternating left-right flickering checkerboards to localize the LGN, **runHemifieldMapping.m**

The code runs on Matlab with Psychtoolbox. Your graphics card must be compatible with the following Psychtoolbox command:
Screen('BlendFunction',win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
***
For more information: www.racheldenison.com
