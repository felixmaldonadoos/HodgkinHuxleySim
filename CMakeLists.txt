cmake_minimum_required(VERSION 3.0.0)
project(implot_demos VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)
add_compile_options(-O3 -march=native)

###############################################################################
# 3RD PARTY LIBS
###############################################################################

# FETCH CONTENT LIBS
include(FetchContent)

FetchContent_Declare(glfw GIT_REPOSITORY https://github.com/glfw/glfw)

set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(glfw)

# LOCAL LIBS

# GLAD
add_subdirectory(3rdparty/glad)

# NFD
add_subdirectory(3rdparty/nativefiledialog-extended)

# EVERYTHING ELSE
include_directories("3rdparty")

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
set(OpenGL_GL_PREFERENCE GLVND)

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
