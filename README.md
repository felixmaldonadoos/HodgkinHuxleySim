Author: by Felix A. Maldonado Osorio
Date: January 6, 2024

# Introduction

### To do: 
- check action potential model. AP is not reaching baseline. 
- voltage seems off. need to recheck equations and what the injection current looks like. maybe indexing is off.
- rendering issues. have seen some bugs online regarding zink and nvidia drivers. 
- uninstalling ```gnuplot-qt``` does not fix rendering issues.

### Notice:
- You may get a weird zink driver issue (I believe with QT dependency for matplotlib). Make sure you have this driver:
```
apt-file search zink_dri.so
libgl1-mesa-dri: /usr/lib/x86_64-linux-gnu/dri/zink_dri.so
```
Error: 
Solved by installing ```xorg-fonts-75dpi and xorg-fonts-100dpi```
|Error|My Solution|
|---|---|
|Error of failed request: BadName (named color or font does not exist) <br> Major opcode of failed request: 45 (X_OpenFont) Serial number of failed request: 21 <br> Current serial number in output stream: 25| Install ```xorg-fonts-75dpi``` and ```xorg-fonts-100dpi```
| err | expl |
# Setup
This is a work in progress. Currently focusing on Ubuntu 22 build in WSL2. 
### My system
```Linux username 5.15.133.1-microsoft-standard-WSL2 #1 SMP Thu Oct 5 21:02:42 UTC 2023 x86_64 x86_64 x86_64 GNU/Linux```

### Installation:
- Clone this repo recursively to include all submodules.

``` 
git clone --recurse-submodules https://github.com/felixmaldonadoos/HodgkinHuxleySim.git 
sudo apt-get install libboost-all-dev
apt install gnuplot-qt
apt install cmake
apt install ninja-build
```

### Dependencies
- CMake (minimum 3.0.0)
- Ninja
- Imgui 
- Implot
- Boost
- gnuplot-qt (required by sciplot for use in tests/HH.cpp)

Note: To run imgui's opengl and glfw implementation on WSL2 Ubuntu 22 but I still had to install ``` Windows SDK 8.1 ``` (found in ```resources/sdksetup.exe``` folder) and ```MSVC v140 - Vs 2015 C++ build tools (v14.00)```.
# Results
### Current version has: 
- HH.cpp creates plots showing gating values and action potential. Model is still off. 

### Demos (to do)

|Demo|Description|Image|
|---|---|---|
|`demo/main.cpp`|Main demo. Displays the ImPlot (and ImGui) library supplied demo windows.| |
|`tests/HH.cpp`|Simple filter toy for educational purposes. Displays time domain input/output signals, and the frequency domain transfer function, amplitude spectrum, etc.|![filter](https://raw.githubusercontent.com/epezent/implot_demos/master/screenshots/filter.png)|


### Resources/Acknowledgements
- https://github.com/Daniel-M/Hodgking-Huxley
- https://github.com/epezent/implot_demos

