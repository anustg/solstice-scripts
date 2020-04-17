import os
import numpy as N
import platform
import os, sys, subprocess, glob, datetime
import colorama
colorama.init()

from .process_raw import *
from .find_solstice import *

def yellow(text):
    return colorama.Fore.YELLOW + colorama.Style.BRIGHT + text + colorama.Style.RESET_ALL

def green(text):
    return colorama.Fore.GREEN + colorama.Style.BRIGHT + text + colorama.Style.RESET_ALL


class Master:

    def __init__(self, casedir='.'):
        if platform.system()=='Linux':
            self.root="~/Solstice-0.9.0-GNU-Linux64"
        else:
            self.root=find_solstice_root()
        self.casedir=casedir
        if not os.path.exists(self.casedir):
	        os.makedirs(self.casedir)
	        assert os.path.isdir(casedir)
        sys.stderr.write("Case directory is '%s'\n" % (yellow(self.casedir),))

    def SPROG(self, n):
        return os.path.join(self.root,'bin',n)

    def in_case(self, fn):
        return os.path.join(self.casedir,fn)


    def run_prog(self, prog, args, output_file=None):
        args1 = [str(a) for a in args]
        sys.stderr.write("Running '%s' with args: %s\n" % (prog," ".join(args1)))
        if output_file is not None:
		    # any error will cause an exception (and we capture the output to a file)
            if platform.system()=='Linux':
                os.system(" ".join(args1))
            else:
                res = subprocess.check_output([prog]+args1)
                with open(output_file,'w') as f:
                    f.write(res.decode('ascii'))
        else:
            if platform.system()=='Linux':
                os.system(" ".join(args1))
            else:
		        # any error will cause an exception...
                subprocess.check_call([prog]+args1)


    def run(self, azimuth, elevation, num_rays, rho_mirror, dni):

        SOLSTICE=self.SPROG("solstice")

        YAML_IN = self.in_case('input.yaml')
        RECV_IN = self.in_case('input-rcv.yaml')

        # main raytrace
        self.run_prog(SOLSTICE,['-D%s,%s'%(azimuth,elevation),'-v','-n',num_rays,'-R',RECV_IN,'-fo',self.in_case('simul'),YAML_IN])
        # post processing
        self.run_prog(SOLSTICE,['-D%s,%s'%(azimuth,elevation),'-g','format=obj:split=geometry','-fo',self.in_case('geom'),YAML_IN])
        self.run_prog(SOLSTICE,['-D%s,%s'%(azimuth,elevation),'-q','-n','100','-R',RECV_IN,'-p','default',YAML_IN], output_file=self.in_case('solpaths'))

        if platform.system()=='Linux':

            os.system('gcc %s/postproc/solppraw.c -o %s/postproc/solppraw'%(self.root, self.root))
            os.system('%s/postproc/solppraw %s'%(self.root, self.in_case('simul')))
            #Read "simul" results and produce receiver files (.vtk) of incoming and/or absorbed solar flux per-primitive
            os.system('gcc %s/postproc/solmaps.c -o %s/postproc/solmaps'%(self.root, self.root))
            os.system('%s/postproc/solmaps %s'%(self.root, self.in_case('simul')))

            #Read "geom" and "simul" file results and produce primaries and receivers files (.vtk), and .obj geometry files
            os.system('gcc %s/postproc/solpp.c -o %s/postproc/solpp'%(self.root, self.root))
            os.system('%s/postproc/solpp %s %s'%(self.root, self.in_case('geom'), self.in_case('simul')))

            #Read "solpaths" file and produce readable file (.vtk) by paraview to visualize the ray paths
            os.system('gcc %s/postproc/solpaths.c -o %s/postproc/solpath'%(self.root, self.root))
            os.system('%s/postproc/solpath %s'%(self.root, self.in_case('solpaths')))
	
            os.system('mv *vtk %s'%self.casedir)
            os.system('mv *obj %s'%self.casedir)
            os.system('mv *txt %s'%self.casedir)

        else:

            # Read "simul" results and produce a text file with the raw results
            run_prog(SPROG('solppraw'),[self.in_case('simul')])

            # Read "simul" results and produce receiver files (.vtk) of incoming and/or absorbed solar flux per-primitive
            run_prog(SPROG('solmaps'),[self.in_case('simul')])

            # Read "geom" and "simul" file results and produce primaries and receivers files (.vtk), and .obj geometry files
            run_prog(SPROG('solpp'),[self.in_case('geom'),self.in_case('simul')])

            # Read "solpaths" file and produce readable file (.vtk) by paraview to visualize the ray paths
            run_prog(SPROG('solpaths'),[self.in_case('solpaths')])

	
            os.system('move *vtk %s >nul'%self.casedir)
            os.system('move *obj %s >nul'%self.casedir)
            os.system('move *txt %s >nul'%self.casedir)

        eta=process_raw_results(self.in_case('simul'), self.casedir,rho_mirror, dni)
        sys.stderr.write('\n' + yellow("Total efficiency: %s\n"%(repr(eta),)))
        sys.stderr.write(green("Completed successfully.\n"))
