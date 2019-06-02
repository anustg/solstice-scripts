ECHO OFF

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
SET /P case_dir= < %solstice_dir%/casedir.input
SET /P azimuth= < %case_dir%/azimuth.input
SET /P elevation= < %case_dir%/elevation.input
SET /P num_rays= < %case_dir%/rays.input
SET /P rho_mirror= <%case_dir%/mirror.input
::
:: Run Solstice
solstice -D%azimuth%,%elevation% -v -n %num_rays% -R %case_dir%/input-rcv.yaml -fo %case_dir%/simul %case_dir%/input.yaml
::
:: Postprocessing
:: make the pp scripts
::mingw32-make %solstice_dir%/src-win/pp_c
::
solstice -D%azimuth%,%elevation% -g format=obj:split=geometry -fo %case_dir%/geom %case_dir%/input.yaml
solstice -D%azimuth%,%elevation% -q -n 100 -R %case_dir%/input-rcv.yaml -p default %case_dir%/input.yaml > %case_dir%/solpaths
::
:: postprocessing in C (provided by Cyril Caliot)
:: Read "simul" results and produce a text file with the raw results
%solstice_dir%/src-win/ppWin/solppraw.exe %case_dir%/simul

:: Read "simul" results and produce receiver files (.vtk) of incoming and/or absorbed solar flux per-primitive
%solstice_dir%/src-win/ppWin/solmaps.exe %case_dir%/simul

:: Read "geom" and "simul" file results and produce primaries and receivers files (.vtk), and .obj geometry files
%solstice_dir%/src-win/ppWin/solpp.exe %case_dir%/geom %case_dir%/simul

::Read "solpaths" file and produce readable file (.vtk) by paraview to visualize the ray paths

%solstice_dir%/src-win/ppWin/solpath.exe %case_dir%/solpaths


:: postprocessing
SET rawfile=%case_dir%/simul
python %solstice_dir%/src-win/srcPy/get_raw.py %rawfile% %case_dir% %rho_mirror%


move *vtk %case_dir% >nul
move *obj %case_dir% >nul
move *txt %case_dir% >nul
cd %case_dir% 
del *.input >nul
cd %solstice_dir%
del *.input >nul
PAUSE
