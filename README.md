Author: by Felix A. Maldonado Osorio
Date: January 6, 2024

# Introduction

### To do: 
- [PRIORITY] Unusual and non-reproducible rendering issue when i run on my personal pc wsl2 ubuntu 22.04.
- check action potential model. AP is not reaching baseline. 
- voltage seems off. need to recheck equations and what the injection current looks like. maybe indexing is off.
- rendering issues. have seen some bugs online regarding zink and nvidia drivers. 
- uninstalling ```gnuplot-qt``` does not fix rendering issues.

# Setup
This is a work in progress. Currently focusing on Ubuntu 22 build in WSL2. 
### My system
```Linux username 5.15.133.1-microsoft-standard-WSL2 #1 SMP Thu Oct 5 21:02:42 UTC 2023 x86_64 x86_64 x86_64 GNU/Linux```

### Installation:
1. Clone this repo recursively to include all submodules.

``` 
git clone --recurse-submodules git@github.com:felixmaldonadoos/HodgkinHuxleySim.git
```
2. Install all packages with one of the following:
#### Method A: Installs all packages I have while developing (not the best but recommended for now)

This Installs all packages I have while developing (not the best but recommended for now). I installed a fresh WSL2 subsystem on my PC, ran the 
```sudo apt-get install etc etc``` line below and all the resulting packages are found in ```installed-packages-with-versions.txt```

```
xargs -a installed-packages-with-versions.txt sudo apt-get install -y
```

#### Method B: Install all packages I added directly one by one.

```
sudo apt-get install libboost-all-dev gnuplot-qt cmake ninja-build xorg-dev build-essential libgtk-3-dev libssl-dev
```


### Notice:
|Error|My Solution|
|---|---|
| CMake Error of failed request: BadName (named color or font does not exist) <br> Major opcode of failed request: 45 (X_OpenFont) Serial number of failed request: 21 <br> Current serial number in output stream: 25| Install ```xorg-fonts-75dpi``` and ```xorg-fonts-100dpi```
|CMake Error at /usr/share/cmake-3.22/Modules/FindPackageHandleStandardArgs.cmake:230 (message): Could NOT find X11 (missing: X11_X11_INCLUDE_PATH X11_X11_LIB)<br> /usr/share/cmake-3.22/Modules/FindPackageHandleStandardArgs.cmake| ```apt install xorg-dev``` |
|CMake Error at /usr/share/cmake-3.22/Modules/FindPackageHandleStandardArgs.cmake:230 (message): Could NOT find OpenSSL, try to set the path to OpenSSL root <br>folder in the system variable OPENSSL_ROOT_DIR (missing: OPENSSL_CRYPTO_LIBRARY OPENSSL_INCLUDE_DIR) /usr/share/cmake-3.22/Modules/FindPackageHandleStandardArgs.cmake |```apt-get install libssl-dev```|
|Checking for module 'gtk+-3.0' 1> [CMake] -- No package 'gtk+-3.0' found |```apt-get install build-essential libgtk-3-dev```|
|gdp not found|```sudo apt-get install build-essential gdb```|
|You may get a weird zink driver issue (I believe with QT dependency for matplotlib). | ```apt-file search zink_dri.so libgl1-mesa-dri: /usr/lib/x86_64-linux-gnu/dri/zink_dri.so```|



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

