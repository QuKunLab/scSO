

- [MATLAB GUI](#matlab-gui)
    - [Installation](#installation)
    - [Tutorial](#tutorial)
- [Python versions](#python-versions)
    - [Installation in linux](#installation-in-linux)
    - [Installation in Windows](#installation-in-windows)
- [Help](#help)

scSO
===
![](./images/Fig.jpg)

scSO is a new algorithm for scRNA-seq data clustering based on Sparse Optimization and low-rank matrix factorization. 
## MATLAB GUI
### Installation

Because scSO relies on the [cvx](http://cvxr.com/cvx/download/) optimization toolbox, you need to download the [cvx](http://cvxr.com/cvx/download/) toolbox and then the extraction path is set to "C: \ Software \ cvx". After downloading the project, first set the "icon" folder to the MATLAB search path. Then just double-click "scSO_GUI.mlapp" to start running scSO

### Tutorial

Note: If it is the first time to cluster scRNA-seq data, it is recommended to select "Coarse". If you want to further subdivide the cell population, we recommend choosing "Fine". If using scSO to clustering 10x data, one should select three files matrix.mtx , barcodes.tsv and genes.tsv at the same time, and upload them to scSO.

<div align=center>
<img src= './images/scSO_app.png'  width="95%" height="50%"  />
</div>

## Python versions

### Installation in linux

Because scSO relies on the Gaussion mixture model in MATLAB, you need to download the [MATLAB Runtime](https://ssd.mathworks.com/supportfiles/downloads/R2019b/Release/7/deployment_files/installer/complete/glnxa64/MATLAB_Runtime_R2019b_Update_7_glnxa64.zip) (The MATLAB Runtime is a freely standalone set of shared libraries that enables the execution of compiled MATLAB applications or components). The detailed installation process is as follows：

    1. unzip MATLAB_Runtime_R2019b_Update_7_glnxa64.zip 
    2. cd MATLAB_Runtime_R2019b_Update_7_glnxa64
    3. sudo ./install
    4. setup MATLAB_Runtime_R2019b_Update_7_glnxa64 in "/usr/local/MATLAB/MATLAB_Runtimecd" 
    5. add "export LD_LIBRARY_PATH=/usr/local/MATLAB/MATLAB_Runtime/v97/runtime/glnxa64:/usr/local/MATLAB/MATLAB_Runtime/v97/bin/glnxa64:/usr/local/MATLAB/MATLAB_Runtime/v97/sys/os/glnxa64:/usr/local/MATLAB/MATLAB_Runtime/v97/sys/opengl/lib/glnxa64" to ~/.bashrc

Next, install the python package that scSO depends on

    pip3 install numpy==1.19.3
    pip3 install pandas
    pip3 install scipy
    pip3 install cvxpy
    pip3 install matplotlib
    pip3 install sklearn

If the above operations are successfully completed, you can execute the sample file of scSO: [test_scSO.ipynb](./scSO_py/Linux/test_scSO.ipynb)

    Attention: scSO was compiled under the environment of Ubuntu20.04 and python3.8.6. If you want to use scSO in other environments, please do the following before using scSO:
        1. pip3 install pybind11
        2. Modify the path in the file "CMakeLists.txt" based on the path of pybind11
        3. mkdir build
        4. cd build
        5. cmake ..
        6. make


### Installation in Windows 

Because scSO relies on the Gaussion mixture model in MATLAB, you need to download the [MATLAB Runtime](https://ssd.mathworks.com/supportfiles/downloads/R2019b/Release/7/deployment_files/installer/complete/win64/MATLAB_Runtime_R2019b_Update_7_win64.zip) (The MATLAB Runtime is a freely standalone set of shared libraries that enables the execution of compiled MATLAB applications or components.) and install with default settings. 

Next, install the python package that scSO depends on

    pip3 install numpy==1.19.3
    pip3 install pandas
    pip3 install scipy
    pip3 install cvxpy
    pip3 install matplotlib
    pip3 install sklearn

If the above operations are successfully completed, you can execute the sample file of scSO: [test_scSO.ipynb](./scSO_py/Windows/test_scSO.ipynb)

    Attention: scSO was compiled under the environment of Windows10, vs2017 and python 3.8.6 If you want to use scSO in other environments, please do the following before using scSO:
        1. pip3 install pybind11
        2. Install Visual Studio(https://visualstudio.microsoft.com/zh-hans/vs)
        3. Modify the path in the file "setup.py"
        4. Execute "python .\setup.py build_ext --inplace" to recompile scSO.

## Help
Click [here](https://git.ustc.edu.cn/hyl2016/scso_testdata/-/raw/master/test.zip?inline=false) to download the test data in the article. If you have any questions or require assistance using scSO, please contact us: hyl2016@mail.ustc.edu.cn .