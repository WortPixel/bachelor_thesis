## "PaDIF - A Particle Detector Interaction Framework for cosmic rays" - a bachelor thesis

Within the scope of this bachelor thesis a new framework for simulations of particle-detector-interactions and based upon that two programs emulating first measurement processes are created. While the framework enables users to build various kinds of interaction settings including custom detectors, the developed programs focus on measuring cosmic particles in ice.

## Required programs (and versions that work for me):  
* GNU Make 3.81  
* cmake 2.8.10.2  
* CINT/ROOT C/C++ Interpreter 5.18.00  

## Used libraries:  
* boost 1.53.0  
* ROOT 5.34/07  

## How to use:
If you want to test EMSA or iceplay make sure, that all the listed programs and libraries are working.
It might be possible, that you need to update the CMakeLists.txt with library locations.

It is recommended to create a new folder for each build. You might want to do something similar to this:  
`$ mkdir build`  
`$ cd build`  
`$ cmake ../.`  
`$ make`  
`$ ./main`  

The last command starts a measurement. For further information have a look in the (soon to be published) thesis or the documentation included in the source code. If that still leaves you with unanswered questions feel free to contact me. My mail address has the structure givenname@surname.de whereas my name is included in the source files.

## Licence:
The software published in the repository is available under the GPL v3 licence.
