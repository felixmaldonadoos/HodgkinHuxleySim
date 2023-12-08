#!/bin/sh

cd implot_demos

# Check if the build directory exists
if [ -d "build" ]; then
    # Remove it if it does
    echo "Removing existing build directory..."
    rm -rf build
fi

# Create a fresh build directory
echo "Creating new build directory..."
mkdir build
cd build

# Run CMake commands
cmake ..
cmake --build . --config Release
