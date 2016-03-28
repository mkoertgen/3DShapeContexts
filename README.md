# 3D Shape Matching with 3D Shape Contexts 

Accompanying source code for my 2003 paper [3D Shape Matching with 3D Shape Contexts](http://cg.cs.uni-bonn.de/en/publications/paper-details/koertgen-2003-3d/).
For historical reasons i made the code freely available here on GitHub. 
I did some tweaks here and there to ensure it builds & runs fine with current compilers (VS2015 C++).
Other than that the original source code from 2003 was not modified.

## Build & Run

Given you have installed Visual Studio 2015, open a command line and run the build script

    build.bat [/t:Clean;Build /v:m]

The project should build fine. Then change to the output directory

    cd 3DShapeContexts\Release

and run 

    3DShapeContexts.exe ANT.obj ANT.obj

This will match the 3D model `ANT.obj` with itself and display the results.
