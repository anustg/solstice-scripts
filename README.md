# Introduction
This repository contains the wrappers for running the Solstice ray-tracing software.

(The development is in progress ...)

# Solstice Installation 
* In Linux system (e.g. Ubuntu 16.04):

```bash
wget 'https://www.meso-star.com/projects/solstice/downloads/Solstice-0.8.2-GNU-Linux64.tar.gz'
tar xzf ~/Solstice-0.8.2-GNU-Linux64.tar.gz
rm ~/Solstice-0.8.2-GNU-Linux64.tar.gz
```

* In Windows system (e.g. Windows 10, 64-bit):

    (coming soon ... )


# Running the wrapper scripts in this repository
* In Linux system (e.g. Ubuntu 16.04):

(1) download the repository and save it in the Solstice directory, e.g.

    ~/Solstice-0.8.2-GNU-Linux64

    or clone the develop branch in the Solstice directory
    
```bash
git clone https://github.com/anustg/solstice-scripts.git
git checkout develop
```
    

(2) Setup the python source scripts
```bash
cd ~/Solstice-0.8.2-GNU-Linux64
sudo python setup.py install
```

(3) Initialisation

    3.1 initialise the Solstice directory in the '.\runSolstice\run.py' Python script

```bash
cd ~/Solstice-0.8.2-GNU-Linux64/runSolstice
gedit run.py
```
    Set the 'solstice_dir' = 'your Solstice directory'

    3.2 Define the simulation case
```bash
gedit set_case.py
```
    The set_case.py is the script that initialises the Solstice directory, casefolder, detailed parameters of the case: (I) the sun, (II) the field, (III) the target  

    -- the example case is a solar tower system, including a heliostat field and a billboard receiver       

    -- the example parameters are the case C1.1 in Wang et al., 2019 (in progress), for verification purposes

    -- users are welcome to define their own test cases by changing these parameters


(4) Run
```bash
python run.py
```

(5) Visualise the results in Paraview, e.g.
```bash
~/ParaView-5.6.0-MPI-Linux-64bit/bin/paraview 
```



# Reference
* **Solstice**: https://www.meso-star.com/projects/solstice/solstice.html
* Wang, Y., Potter, D., Asselineau, C.-A., Corsi, C., Wagner, M., Caliot, C., Piaud, B., Blanco, M., Kim, J.-S., Pye, J., 2019. Verification of Optical Modeling on Sunshape and Surface Slope Error. Solar Energy  – in progress



