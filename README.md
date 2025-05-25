# Introduction

This repository contains Python wrappers for running the  [Solstice](https://www.meso-star.com/projects/solstice/solstice.html) ray-tracing software, which can be used for the simulation of concentrating solar-thermal power (CSP) systems. The repository also contains copies of the [post-processing tools](https://www.meso-star.com/projects/solstice/solstice-resources.html) for Solstice, which are required by our Python wrappers. Finally, this repository contains code to generate an easy-to-use Windows installer for Solstice.

[![github](https://readthedocs.org/projects/solsticepy/badge/?version=latest)][1]

# 1. solsticepy

See the [documentation](https://solsticepy.readthedocs.io/en/latest/?badge=latest) for `solsticepy` at Read The Docs.

# 2. Install Solstice

## 2.1 Windows

* For Windows (Windows 10, 64-bit), use our 64-bit [binary installer](https://github.com/anustg/solstice-scripts/releases) to install Solstice. The resulting files will be located in `c:\Program Files\solstice-0.9.0` (with the appropriate version number). Note that the Windows installer **also includes** the post-processing tools.

## 2.2 Linux

* In Linux system (eg Ubuntu 18.04 64-bit), you need

(1) download the binary tarball from the [Solstice homepage](https://www.meso-star.com/projects/solstice/solstice.html). 

(2) set up the post-processing tools. 

### 2.2.1 Install solstice

```bash
### download and extract:
cd ~
wget "https://www.meso-star.com/projects/solstice/downloads/Solstice-0.9.0-GNU-Linux64.tar.gz"
tar zxvf Solstice-0.9.0-GNU-Linux64.tar.gz
### add solstice to your PATHs:
export PATH=$PATH:~/Solstice-0.9.0-GNU-Linux64/bin
export LD_LIBRARY_PATH=$LB_LIBRARY_PATH:~/Solstice-0.9.0-GNU-Linux64/lib
### check that it runs:
solstice --help
```

Note that you can add the two `export` commands to your `~/.profile` file in order that `solstice` remains available next time you log in.

### 2.2.2 Install the post-processing program

Perform the following comments to clone the solsticepy package, and compile the post-processing program by scons

```bash
cd ~
git clone https://github.com/anustg/solstice-scripts.git
cd solstice-scripts/postproc
scons
```

To install the post-processing programs in a user-defined path, replace the last step by the following comment, and then add the /user/defined/path to system path if needed.

```bash
scons INSTALL_PREFIX=/user/defined/path install
export PATH=$PATH:/user/defined/path

```

# 3. Install `solsticepy` 

The instructions below give easy instructions for end-users. If you are interested in hacking/developing/changing `solsticepy` code, see [HACKING](HACKING.md) instead.

## Windows

* First install Python if you have not done so already. Choose the latest "Windows x86-64 executable installer" of the latest Python 3 [release for Windows](https://www.python.org/downloads/windows/).

* If you added Python to your PATH during installation you can just open a command prompt and type `pip3 install solsticepy` and that should download and install everything you need.

* If you didn't add Python to your PATH, you can do the following:

  * Create a text file `install.py` in your home directory with the following content:

```python
# put this file in your home directory and name it install.py
import sys, subprocess
subprocess.check_call([sys.executable,"-m","pip","install","solsticepy"])
```
* Open a command prompt (type 'cmd' in the Start menu), then type

```
install.py
```
* You should see Python downloading and installing 'solsticepy' from the PyPI servers.

# 4. Install Paraview

Pavaview is a software to visualise Solstice simulation scene.

* **Be sure to install a version 4 release of Paraview**. Version 5 was not stable on Windows when we tested it (Apr 2020).
* On Windows:

  * Try downloading and installing [ParaView-4.4.0-Qt4-Windows-64bit.exe](https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v4.4&type=binary&os=Windows&downloadFile=ParaView-4.4.0-Qt4-Windows-64bit.exe).

* On Linux (Ubuntu 18.04):

  * `sudo apt install paraview` should be all you need.

# 5. Running an example:

* **This section yet to be completed**

* Download our Zip file from [this page](https://github.com/anustg/solstice-scripts/releases/)

* Extract the folder `example` to your home directory

* Open a command prompt, navigate to your example directory, then run the `run.py` script:

```
cd example
run.py
```
* You should see various output, followed by "Completed successfully". Also note the output what says "Case directory is...".

* In the case directory, you will file output files including `.csv` files that can be opened using Excel, and `.vtk` files that contain 3D graphics can be opened in Paraview.

# References

* **Solstice**: https://www.meso-star.com/projects/solstice/solstice.html

* Wang, Y., Potter, D., Asselineau, C.-A., Corsi, C., Wagner, M., Caliot, C., Piaud, B., Blanco, M., Kim, J.-S., Pye, J., 2019. [Verification of Optical Modeling on Sunshape and Surface Slope Error](https://www.researchgate.net/publication/337636543). Solar Energy 195. doi:10.1016/j.solener.2019.11.035

[1]: https://solsticepy.readthedocs.io/en/latest/?badge=latest




