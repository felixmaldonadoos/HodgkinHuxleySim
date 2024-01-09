cmake_minimum_required(VERSION 3.0.0)
project(implot_demos VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)
add_compile_options(-O3 -march=native)

###############################################################################
# 3RD PARTY LIBS
###############################################################################
# FETCH CONTENT LIBS
include(FetchContent) 
# set(FETCHCONTENT_FULLY_DISCONNECTED ON)

FetchContent_Declare(fmt GIT_REPOSITORY https://github.com/fmtlib/fmt)
FetchContent_MakeAvailable(fmt)

FetchContent_Declare(glfw GIT_REPOSITORY https://github.com/glfw/glfw) 
FetchContent_MakeAvailable(glfw)

FetchContent_Declare(iir GIT_REPOSITORY https://github.com/berndporr/iir1) 
FetchContent_MakeAvailable(iir)

if (WIN32)
  set(CMAKE_USE_SCHANNEL ON)
endif()
FetchContent_Declare(curl
URL                    https://github.com/curl/curl/releases/download/curl-7_75_0/curl-7.75.0.tar.xz
URL_HASH               SHA256=fe0c49d8468249000bda75bcfdf9e30ff7e9a86d35f1a21f428d79c389d55675 # the file hash for curl-7.75.0.tar.xz
USES_TERMINAL_DOWNLOAD TRUE) 
FetchContent_MakeAvailable(curl)
add_library(curl_int INTERFACE)
target_link_libraries(curl_int INTERFACE libcurl)
target_include_directories(curl_int INTERFACE ${curl_SOURCE_DIR}/include ${curl_BINARY_DIR}/include/curl)
add_library(CURL::libcurl ALIAS curl_int)
# LOCAL LIBS

# GLAD
add_subdirectory(3rdparty/glad)

# NFD
add_subdirectory(3rdparty/nativefiledialog-extended)

# KISS FFT
set(KISS_SRC "3rdparty/kissfft/kiss_fft.h" 
             "3rdparty/kissfft/kiss_fftr.h" 
             "3rdparty/kissfft/kiss_fft.c" 
             "3rdparty/kissfft/kiss_fftr.c")

# EVERYTHING ELSE
include_directories("3rdparty" "3rdparty/kissfft")

###############################################################################
# IMGUI
###############################################################################

set(LIBS_PATH "${CMAKE_CURRENT_SOURCE_DIR}/libs")
set(IMGUI_PATH ${LIBS_PATH}/imgui)
set(IMGUI_HEADERS
    ${IMGUI_PATH}/imconfig.h
    ${IMGUI_PATH}/imgui_internal.h
    ${IMGUI_PATH}/imgui.h
    ${IMGUI_PATH}/imstb_rectpack.h
    ${IMGUI_PATH}/imstb_textedit.h
    ${IMGUI_PATH}/imstb_truetype.h    
    ${IMGUI_PATH}/backends/imgui_impl_glfw.h
    ${IMGUI_PATH}/backends/imgui_impl_opengl3.h
    ${IMGUI_PATH}/backends/imgui_impl_opengl3_loader.h
    ${IMGUI_PATH}/misc/cpp/imgui_stdlib.h
)

set(IMGUI_SRC
    ${IMGUI_PATH}/imgui_demo.cpp
    ${IMGUI_PATH}/imgui_draw.cpp
    ${IMGUI_PATH}/imgui_tables.cpp
    ${IMGUI_PATH}/imgui_widgets.cpp
    ${IMGUI_PATH}/imgui.cpp
    ${IMGUI_PATH}/backends/imgui_impl_glfw.cpp
    ${IMGUI_PATH}/backends/imgui_impl_opengl3.cpp
    ${IMGUI_PATH}/misc/cpp/imgui_stdlib.cpp
)

find_package(OpenGL REQUIRED)

add_library(imgui ${IMGUI_HEADERS} ${IMGUI_SRC})
if(MSVC)
  target_compile_options(imgui PRIVATE /W4 /WX /arch:AVX2 /fp:fast)
endif()

if (WIN32)
  target_link_libraries(imgui PUBLIC glfw glad OpenGL::GL imm32)
else()
  target_link_libraries(imgui PUBLIC glfw glad OpenGL::GL)
endif()
target_compile_definitions(imgui PRIVATE IMGUI_DLL_EXPORT)

include_directories(${IMGUI_PATH} ${IMGUI_PATH}/examples ${IMGUI_PATH}/examples/libs/gl3w ${IMGUI_PATH}/backends ${IMGUI_PATH}/misc/cpp)

# imnodes
add_library(imnodes 3rdparty/imnodes/imnodes.h 3rdparty/imnodes/imnodes_internal.h 3rdparty/imnodes/imnodes.cpp)
target_link_libraries(imnodes imgui)
include_directories(3rdparty/imnodes)

###############################################################################
# IMPLOT
###############################################################################
set(IMPLOT_PATH ${LIBS_PATH}/implot)


set(IMPLOT_HEADERS ${IMPLOT_PATH}/implot.h ${IMPLOT_PATH}/implot_internal.h)
set(IMPLOT_SRC     ${IMPLOT_PATH}/implot.cpp ${IMPLOT_PATH}/implot_items.cpp ${IMPLOT_PATH}/implot_demo.cpp) 

if(EXISTS "${IMPLOT_PATH}/backends")
list(APPEND IMPLOT_HEADERS ${IMPLOT_PATH}/backends/implot_backend.h ${IMPLOT_PATH}/backends/implot_impl_opengl3.h)
list(APPEND IMPLOT_SRC ${IMPLOT_PATH}/backends/implot_backend.h ${IMPLOT_PATH}/backends/implot_impl_opengl3.cpp)
endif()

add_library(implot ${IMPLOT_HEADERS} ${IMPLOT_SRC})
target_link_libraries(implot PUBLIC imgui)
target_compile_definitions(implot PUBLIC IMPLOT_DEBUG IMPLOT_DLL_EXPORT IMPLOT_BACKEND_ENABLE_OPENGL3 IMGUI_IMPL_OPENGL_LOADER_GLAD)
set_property(TARGET implot PROPERTY CXX_STANDARD 11)
if(MSVC)
  target_compile_options(implot PRIVATE /W4 /WX /arch:AVX2 /fp:fast /permissive-) 
else()
  target_compile_options(implot PRIVATE -Wall -Wextra -pedantic -Werror -mavx2 -Ofast)
endif()

include_directories(${IMPLOT_PATH}/)

###############################################################################
# APP FRAMEWORK
###############################################################################

add_library(app 
  common/App.h
  common/App.cpp
  common/Shader.h
  common/Native.h
  common/Native.cpp
  common/Image.h
  common/Fonts/Fonts.h
  common/Fonts/FontAwesome5.cpp
  common/Fonts/RobotoRegular.cpp
  common/Fonts/RobotoBold.cpp
  common/Fonts/RobotoItalic.cpp
  common/Fonts/RobotoMonoRegular.cpp
  common/Fonts/RobotoMonoBold.cpp
  common/Fonts/RobotoMonoItalic.cpp
)
target_include_directories(app PUBLIC common)
target_link_libraries(app implot nfd)
target_compile_features(app PRIVATE cxx_std_17)

###############################################################################
# DEMO APPS
###############################################################################

# main ImPlot demo
add_executable(demo "demos/demo.cpp")
target_link_libraries(demo app)

# FIND BOOST 
find_package(Boost 1.50 REQUIRED)

set(TESTS_PATH "${CMAKE_CURRENT_SOURCE_DIR}/tests")
add_executable(lorenz "${TESTS_PATH}/lorenz.cpp")
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(lorenz ${Boost_LIBRARIES})
target_compile_features(lorenz PRIVATE cxx_std_17)

# LOOK FOR HELPERS
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/helpers)

add_executable(HH "${TESTS_PATH}/HH.cpp" 
    ${CMAKE_CURRENT_SOURCE_DIR}/helpers/gatingfunctions.h 
    ${CMAKE_CURRENT_SOURCE_DIR}/helpers/gatingfunctions.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/helpers/Models.h
    ${CMAKE_CURRENT_SOURCE_DIR}/helpers/Models.cpp
)
target_link_libraries(HH ${Boost_LIBRARIES} imgui)
target_compile_features(HH PRIVATE cxx_std_17)

# ===== main demo ====

add_executable(HodkingHuxleySim 
    "demos/main.cpp" 
    "common/Plot.cpp" 
    "common/Plot.h" 
    ${KISS_SRC} 
    ${Boost_LIBRARIES})              
target_link_libraries(HodkingHuxleySim app)        
target_compile_features(HodkingHuxleySim PRIVATE cxx_std_17)