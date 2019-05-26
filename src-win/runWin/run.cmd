 

:: Initialise
:: ==========================================
::
SET solstice_dir=E:\Solstice-0.8.2-Win64
::
:: =========================================
::
::
::
python set_case.py
::
SET case_dir= < %solstice_dir%/casedir.input
SET /P azimuth= < %case_dir%/azimuth.input
SET /P elevation= < %case_dir%/elevation.input
SET /P num_rays= < %case_dir%/rays.input
::
:: Run Solstice
solstice -D%azimuth%,%elevation% -v -n %num_rays% -R %case_dir%/input-rcv.yaml -fo %case_dir%/simul %case_dir%/input.yaml
::
:: Postprocessing
:: make the pp scripts
::mingw32-make %solstice_dir%/src/pp_c
::
solstice -D%azimuth%,%elevation% -g format=obj:split=geometry -fo %case_dir%/geom %case_dir%/input.yaml
solstice -D%azimuth%,%elevation% -q -n 100 -R %case_dir%/input-rcv.yaml -p default %case_dir%/input.yaml > %case_dir%/solpaths
::
:: postprocessing in C (provided by Cyril Caliot)
:: Read "simul" results and produce a text file with the raw results
%solstice_dir%/src/pp_c/solppraw.exe %case_dir%/simul

:: Read "simul" results and produce receiver files (.vtk) of incoming and/or absorbed solar flux per-primitive
%solstice_dir%/src/pp_c/solmaps.exe %case_dir%/simul

:: Read "geom" and "simul" file results and produce primaries and receivers files (.vtk), and .obj geometry files
%solstice_dir%/src/pp_c/solpp.exe %case_dir%/geom %case_dir%/simul

::Read "solpaths" file and produce readable file (.vtk) by paraview to visualize the ray paths

%solstice_dir%/src/pp_c/solpath.exe %case_dir%/solpaths


move *vtk %case_dir%
move *obj %case_dir%
move *txt %case_dir%
del *.input
cd %case_dir%
del *.input

PAUSE
