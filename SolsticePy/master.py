import os
import numpy as N
import platform
from SolsticePy.get_raw import *

def run(azimuth, elevation, num_rays, rho_mirror, case_dir='.', version='0.9.0'):

    system=platform.system() # 'Linux' or 'windows'

    if system=='Linux':
        solstice_dir='~/Solstice-%s-GNU-Linux64'%version
    else:
        from SolsticePy.solsticething import find_solstice_root
        solstice_dir=find_solstice_root()

    os.system("solstice -D%s,%s -v -n %s -R %s/input-rcv.yaml -fo %s/simul %s/input.yaml"%(azimuth, elevation, num_rays, case_dir, case_dir, case_dir))
    #

    # preparing post processing
    os.system('solstice -D%s,%s -g format=obj:split=geometry -fo %s/geom %s/input.yaml'%(azimuth, elevation, case_dir, case_dir))
    os.system('solstice -D%s,%s -q -n 100 -R %s/input-rcv.yaml -p default %s/input.yaml > %s/solpaths'%(azimuth, elevation, case_dir, case_dir, case_dir))
    #
    # postprocessing in C (provided by Cyril Caliot)
    #Read "simul" results and produce a text file with the raw results
    if system=='Linux':

        os.system('gcc %s/ppLinux/solppraw.c -o %s/ppLinux/solppraw'%(solstice_dir, solstice_dir))
        os.system('%s/ppLinux/solppraw %s/simul'%(solstice_dir, case_dir))
        #Read "simul" results and produce receiver files (.vtk) of incoming and/or absorbed solar flux per-primitive
        os.system('gcc %s/ppLinux/solmaps.c -o %s/ppLinux/solmaps'%(solstice_dir, solstice_dir))
        os.system('%s/ppLinux/solmaps %s/simul'%(solstice_dir, case_dir))

        #Read "geom" and "simul" file results and produce primaries and receivers files (.vtk), and .obj geometry files
        os.system('gcc %s/ppLinux/solpp.c -o %s/ppLinux/solpp'%(solstice_dir, solstice_dir))
        os.system('%s/ppLinux/solpp %s/geom %s/simul'%(solstice_dir, case_dir, case_dir))

        #Read "solpaths" file and produce readable file (.vtk) by paraview to visualize the ray paths
        os.system('gcc %s/ppLinux/solpaths.c -o %s/ppLinux/solpath'%(solstice_dir, solstice_dir))
        os.system('%s/ppLinux/solpath %s/solpaths'%(solstice_dir, case_dir))
	
        os.system('mv *vtk %s'%case_dir)
        os.system('mv *obj %s'%case_dir)
        os.system('mv *txt %s'%case_dir)

    else:
        os.system('%s/ppWin/solppraw.exe %s/simul'%(solstice_dir, case_dir))

        #Read "simul" results and produce receiver files (.vtk) of incoming and/or absorbed solar flux per-primitive
        os.system('%s/ppWin/solmaps.exe %s/simul'%(solstice_dir, case_dir))

        #Read "geom" and "simul" file results and produce primaries and receivers files (.vtk), and .obj geometry files
        os.system('%s/ppWin/solpp.exe %s/geom %s/simul'%(solstice_dir, case_dir, case_dir))

        #Read "solpaths" file and produce readable file (.vtk) by paraview to visualize the ray paths
        os.system('%s/ppWin/solpath.exe %s/solpaths'%(solstice_dir, case_dir))

	
        os.system('move *vtk %s >nul'%case_dir)
        os.system('move *obj %s >nul'%case_dir)
        os.system('move *txt %s >nul'%case_dir)

    rawfile=case_dir+'/simul'
    proces_raw_results(rawfile, case_dir,rho_mirror)


