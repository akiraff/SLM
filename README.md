## General Info
This software is used to generate and display phase pattern on a spatial light modulator (SLM). The algorithm used here is phase fixed weighted Gerchberg-Saxton (WGS).
The detailed information is included in the reference: https://www.osapublishing.org/ol/abstract.cfm?uri=ol-44-12-3178. 
## UI Interface
    This software has a UI interface, run by UI_SLM.py. This interface has two mode:
        1 Calculate mode: SLM disconnected, calculate phase pattern with phase fixed WGS.
        2 Display mode: SLM connected, choose caculated phase pattern and send it to SLM screen.
    To toggle between the two mode, one needs to change display mode in the UI_SLM.py script. 
#### UI input
    pixel pitch: this number is the pixel pitch of the SLM.
    Array spacing: Specify the array spacing in the focal plane along x and y direction, in unit of micron.
    Distance from origine: The closest foci point on the focal plane to the zero order spot, in unit of micron.
    Geometry: Foci array geometry. Currently you can choose from square, triangle and kagome. You can add more
              by some straight forward modification in tweezer class in IMG.py.
    rot angle: Rotate the foci array around the zero order spot. You need to check the checkbox Rotate? if you
               want this rotation.
    Array size: Array size along x and along y
    FFT grid (bit): This number is the size of FFT grid to perform WGS, in unit of bit. Larger number will be
                    beneficial for more complicated foci intensity pattern, and also slower. For a regular
                    lattice pattern, 12 is good enough.
    Threshold: This number is for phase fixed method. When the uniformity of the foci array is below this value,
               WGS will stop evolving the phase. This method speeds up the convergence dramatically. Usually 
               0.01-0.05 will work.
    Loop: The number of loops used in WGS. 20 is a good input.
    
