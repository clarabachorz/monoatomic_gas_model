README

Two modules should be present :
-Ball.py
-Gas.py

Both must be uploaded, these are the two classes necessary for the code to run. 
The user can then open a new worksheet. It is necessary to import Gas in the new worksheet for the code to work.
Then the user can have access to all methods defined in these classes. It is recommended to 
create a new Gas object with parameters wished, or the default parameters are called by using Gas().

Gas takes the following parameters :

N- number of molecules (default 20)
T- scaling temperature (default 293)
radius_N -radius of the balls (default 0.5)
container_radius - radius of the container (default 30)
mass_N -mass of the molecule

Use My_Gas.animation(Number_of_frames_wished) to set off an animation. Recommended frame number is 1000.
Use My_Gas.information(time_to_run) to access only the data. Time_to_run is the time in seconds for 
which the gas system will be simulated. 
N.b.This method can take a long time to run if the number of molecules in the gas is very big.

Known bug: Two balls may sometimes stick on collision. 
