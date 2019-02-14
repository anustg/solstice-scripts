# solstice-scripts
Scripts and wrappers for Solstice ray-tracing software

# example
* If need installation of Solstice:
get version 0.8.2
'''bash
wget 'https://www.meso-star.com/projects/solstice/downloads/Solstice-0.8.2-GNU-Linux64.tar.gz'
tar xzf ~/Solstice-0.8.2-GNU-Linux64.tar.gz
rm ~/Solstice-0.8.2-GNU-Linux64.tar.gz
'''

* get Solstice ready to run
```bash
source ~/Solstice-0.8.2-GNU-Linux64/etc/solstice.profile
```
* run solstice by the master script
```bash
python source/SOLSTICE_master.py 
```
* visualise the results in Paraview
```bash
~/ParaView-5.6.0-RC2-Qt5-MPI-Linux-64bit/bin/paraview 
```

**Please note:**
* The example case is a solar tower system, including a heliostat field and a billboard receiver
* The master script is ./source/SOLSTICE_master.py
* Two demonstrations are made in the master script: 
    -(1) single sun postion
    -(2) annual performance simulation

  Users can switch the cases in the master script
* The detail parameters are specified in the file:
    ./source/set_parameter.csv

  The current parameters are set the same as the case C1.1 in (Wang et al., 2019 - in progress)

  Users can specify their own test case by changing the parameters.

# Reference
* **Solstice**: https://www.meso-star.com/projects/solstice/solstice.html
* Wang, Y., Potter, D., Asselineau, C.-A., Corsi, C., Wagner, M., Caliot, C., Piaud, B., Blanco, M., Kim, J.-S., Pye, J., 2019. Verification of Optical Modeling on Sunshape and Surface Slope Error. Solar Energy  – in progress



