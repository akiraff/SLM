## General Info
This software is used to generate and display phase pattern on a spatial light modulator (SLM). The algorithm used here is phase fixed weighted Gerchberg-Saxton (WGS).
The detailed information is included in the reference: https://www.osapublishing.org/ol/abstract.cfm?uri=ol-44-12-3178. 
## UI Interface
    This software has a UI interface, run by UI_SLM.py. This interface has two mode:
        1 Calculate mode: SLM disconnected, calculate phase pattern with phase fixed WGS.
        2 Display mode: SLM connected, choose caculated phase pattern and send it to SLM screen.
    To toggle between the two mode, one needs to change display mode in the UI_SLM.py script. 
