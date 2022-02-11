This is an implementation of the Paper "Automatic Generation of Load-Balancing-Aware Block-Structured Grids for Complex Ocean Domains" presented at the International Meshing Roundtable 2022.

The provided example uses the ADCIRC mesh of the Bight of Abaco in the Bahamas region as an input. [[Grenier Jr. et al.]](https://doi.org/10.1029/95JC00841)

## Requirements Linux:
(install in order as listed)

### C++17

### CMAKE 
* cmake version ≥ 3.10
* follow installation guide: https://cmake.org/install/

### CUDA
[Requirements](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#system-requirements)
* CUDA-capable GPU
* A supported version of Linux with a gcc compiler and toolchain
* NVIDIA CUDA Toolkit

[Installation](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html)
* download [CUDA](https://developer.nvidia.com/cuda-downloads)
* To set path in .bashrc: get cuda installation path with 
    ```
    <cuda_path>=ldconfig -p | grep cuda, e.g. /usr/local/cuda-11.0/  
    export CUDA_ROOT=<cuda_path>/bin/
    export LD_LIBRARY_PATH=<cuda_path>/lib64/
    ```
* Select nvidia graphics card with 
    ```
    sudo prime-select nvidia
    ```
* reboot
* check if OpenGL uses nvidia GPU: 
    ```
    glxinfo | grep -i opengl
    ```


[Verify CUDA Installation](https://xcat-docs.readthedocs.io/en/stable/advanced/gpu/nvidia/verify_cuda_install.html)
* Verify driver version: 
    ```
    cat /proc/driver/nvidia/version
    ```
* Verify the CUDA Toolkit version: 
    ```
    nvcc -V
    ```
* Verify running CUDA GPU jobs (Running `deviceQuery` and `bandwidthTest`)
    ```
    cuda-install-samples-[11.0].sh .
    cd NVIDIA_CUDA-[11.0]_Samples/
    make
    ./bin/[x86_64]/[linux]/release/deviceQuery
    ./bin/[x86_64]/[linux]/release/bandwidthTest
    ```


### LIBS to install
1. [glog](https://github.com/google/glog)
2. [OpenMesh 8.1](https://www.graphics.rwth-aachen.de/software/openmesh/)

* Create a directory for building the libraries (non-root should note the dependency)
    ```
    mkdir 00_libs
    cd 00_libs
    ```
* Clone the repositories and/or download the archives
    ```
    git clone https://github.com/google/glog.git
    wget https://www.graphics.rwth-aachen.de/media/openmesh_static/Releases/8.1/OpenMesh-8.1.tar.gz
    ```

* [Build+Test glog](https://github.com/google/glog#building-from-source)
    (Root users may disregard the -DCMAKE_INSTALL_PREFIX flag)
    ```
    mkdir glog/build
    cd glog/build/
    cmake .. -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX="."
    cmake --build .
    cmake --build . --target test
    cmake --build . --target install
    cd ../..
    ```
* [Build OpenMesh 8.1](https://www.graphics.rwth-aachen.de/media/openmesh_static/Documentations/OpenMesh-Doc-Latest/a04315.html)
     (Root users may disregard the -DCMAKE_INSTALL_PREFIX flag)
    ```
    tar -xvf OpenMesh-8.1.tar.gz
    rm OpenMesh-8.1.tar.gz
    mkdir OpenMesh-8.1/build
    cd OpenMesh-8.1/build
    cmake .. -DCMAKE_INSTALL_PREFIX="."
    cmake --build .
    cmake --build . --target install
    cd ../..
    ```

* Note real path where LIBS have been installed (maybe add it to your bashrc extension if you plan on recompiling the project more often)
    ```
    realpath glog/build/ googletest/build/ OpenMesh-8.1/build/
    export TEMP_LIB_location=$(echo "$(realpath glog/build/);$(realpath OpenMesh-8.1/build)")
    ```
### Build HPMeshGen2
* Clone the repository
    ```
    git clone https://github.com/DanielZint/hpmeshgen.git
    ```
* Building
    ```
    mkdir hpmeshgen2/build
    cd hpmeshgen2/build
    cmake .. -DCMAKE_PREFIX_PATH="$TEMP_LIB_location"
    cmake --build .
    cd ../..
    ```
* Executing with parameters described in ./hpmeshgen2/HPMeshGenParameter.txt
     ```
    ./hpmeshgen2/build/bin/HPMeshGen2
     ```


## Requirements Windows:
(install in order as listed)


### CUDA
* [download](https://developer.nvidia.com/cuda-downloads)
* normal installation, restart device afterwards


### CMAKE
* [download](https://cmake.org/download/)
* select binary with installer


### VCPKG
* [download](https://github.com/Microsoft/vcpkg)
* follow along the [video](https://www.youtube.com/watch?v=b7SdgK7Y510&t=751s)
* libs to install: 
    1. [glog](https://github.com/google/glog)
    2. [OpenMesh 8.1](https://www.graphics.rwth-aachen.de/software/openmesh/)

### OPTIONAL

#### VISUAL STUDIO
* [download](https://visualstudio.microsoft.com/de/downloads/)
* select everything that has to do with C++ or MSVC
* select english language pack


## Usage:
Use file *HPMeshGenParameter.txt* to configure *HPMeshGen2*. The file marks comments with '!'.

```
workingDirectory        ! is the project root folder by default but can be changed here
cacheFolder             ! storage folder for cached files (cache reduces time for recomputing meshes significantly)
nBlocks                 ! number of blocks

meshFileName            ! input mesh, usually a .14-file
sizegridSizeX           ! size grid dimensions (integer)
sizegridSizeY           ! the size grid stores the depth and edge length information

nRefinementSteps        ! the number of nodes on a fragment edge
nPatches                ! number of fragments

outputFolder            ! all output files go in here
initialMeshOutput       ! 0/1 initial mesh as .off file
reductionOutput         ! 0/1 intermediate steps of reduction as .off file
blockMeshOutput         ! 0/1 fragment mesh as .off file
sizegridOutput          ! 0/1 print size grid in .vtk format

forceFragmentNumber     ! 0/1 force fragment number, 0: change fragment number for better quality
convexHullDecimation    ! 0/1 activate masks
```


