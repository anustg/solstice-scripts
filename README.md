# solstice-scripts
Scripts and wrappers for Solstice ray-tracing software

# example

* get Solstice ready
```bash
source ~/Solstice-0.8.1-GNU-Linux64/etc/solstice.profile
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
    (1) single sun postion
    (2) annual performance simulation
  Users can switch the cases in the master script
* The detail parameter inputs are specified in file:
    ./source/set_parameter.csv
  The current parameters are set the same as the case C1.1 in (Wang et al., 2019 - in progress).
  Users can specify their own test case by changing the parameters.



