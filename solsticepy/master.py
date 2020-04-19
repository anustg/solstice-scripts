import os
import numpy as N
import platform
import os, sys, subprocess, glob, datetime
import colorama
colorama.init()

from .process_raw import *
from .find_solstice import *
from .cal_sun import *

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
                os.system(prog+" "+" ".join(args1)+" > "+output_file)
            else:
                res = subprocess.check_output([prog]+args1)
                with open(output_file,'w') as f:
                    f.write(res.decode('ascii'))
        else:
            if platform.system()=='Linux':
                os.system(prog+" "+" ".join(args1))
            else:
		        # any error will cause an exception...
                subprocess.check_call([prog]+args1)


    def run(self, azimuth, elevation, num_rays, rho_mirror,dni):

        SOLSTICE=self.SPROG("solstice")

        YAML_IN = self.in_case('input.yaml')
        RECV_IN = self.in_case('input-rcv.yaml')

        # main raytrace
        self.run_prog(SOLSTICE,['-D%s,%s'%(azimuth,elevation),'-v','-n',num_rays,'-R',RECV_IN,'-fo',self.in_case('simul'),YAML_IN])
        # post processing
        self.run_prog(SOLSTICE,['-D%s,%s'%(azimuth,elevation),'-g','format=obj:split=geometry','-fo',self.in_case('geom'),YAML_IN])

        self.run_prog(SOLSTICE,['-D%s,%s'%(azimuth,elevation),'-q','-n','100','-R',RECV_IN,'-p','default',YAML_IN], output_file=self.in_case('solpaths'))

        # Read "simul" results and produce a text file with the raw results
        self.run_prog(self.SPROG('solppraw'),[self.in_case('simul')])
        # Read "simul" results and produce receiver files (.vtk) of incoming and/or absorbed solar flux per-primitive
        self.run_prog(self.SPROG('solmaps'),[self.in_case('simul')])

        # Read "geom" and "simul" file results and produce primaries and receivers files (.vtk), and .obj geometry files
        self.run_prog(self.SPROG('solpp'),[self.in_case('geom'),self.in_case('simul')])

        # Read "solpaths" file and produce readable file (.vtk) by paraview to visualize the ray paths
        self.run_prog(self.SPROG('solpaths'),[self.in_case('solpaths')])


        if platform.system()=="Linux":
            os.system('mv *vtk '+self.casedir)
            os.system('mv *obj '+self.casedir)
            os.system('mv *txt '+self.casedir)

        else:
            os.system('move *vtk %s >nul'%self.casedir)
            os.system('move *obj %s >nul'%self.casedir)
            os.system('move *txt %s >nul'%self.casedir)

        eta, performance_hst=process_raw_results(self.in_case('simul'), self.casedir,rho_mirror,dni)
        sys.stderr.write('\n' + yellow("Total efficiency: %s\n"%(repr(eta),)))
        sys.stderr.write(green("Completed successfully.\n"))



    def run_annual(self, nd, nh, latitude, num_rays, num_hst,rho_mirror,dni,gen_vtk=False):

        SOLSTICE=self.SPROG("solstice")
        YAML_IN = self.in_case('input.yaml')
        RECV_IN = self.in_case('input-rcv.yaml')

        sun=SunPosition()
        AZI, ZENITH,table,case_list=sun.annual_angles(latitude, casefolder=self.casedir, nd=nd, nh=nh)
        case_list=case_list[1:]
        SOLSTICE_AZI, SOLSTICE_ELE=sun.convert_convention('solstice', AZI, ZENITH)

        # performance of individual heliostat is recorded
        # TODO note, DNI is not varied in the simulation, 
        # i.e. performance is not dni-weighted
        ANNUAL=np.zeros((num_hst, 9))    
        run=N.r_[0]

        for i in range(len(case_list)):     
            c=int(case_list[i,0].astype(float))
            if c not in run:
                azimuth=SOLSTICE_AZI[c-1]
                elevation=SOLSTICE_ELE[c-1]
                #if np.sin(elevation*np.pi/180.)>=1.e-5:
                #	dni=1618.*np.exp(-0.606/(np.sin(elevation*np.pi/180.)**0.491))
                #else:
                #	dni=0.

                print('')
                print('sun position:', (c))
                print('azimuth:',  azimuth, ', elevation:',elevation)

                onesunfolder=self.casedir+'/sunpos_%s'%(c)
                if not os.path.exists(onesunfolder):
                    os.makedirs(onesunfolder) 
                # run solstice
                if elevation<1.:
                    efficiency_total=ufloat(0,0)
                    performance_hst=np.zeros((num_hst, 9))  
                else:

                    # main raytrace
                    os.system("%s -D%s,%s -v -n %s -R %s -fo %s %s"%(SOLSTICE, azimuth, elevation, num_rays, RECV_IN, self.in_case('simul'), YAML_IN))

                    if gen_vtk:

                        # post processing
                        self.run_prog(SOLSTICE,['-D%s,%s'%(azimuth,elevation),'-g','format=obj:split=geometry','-fo',self.in_case('geom'),YAML_IN])

                        self.run_prog(SOLSTICE,['-D%s,%s'%(azimuth,elevation),'-q','-n','100','-R',RECV_IN,'-p','default',YAML_IN], output_file=self.in_case('solpaths'))
                        # Read "simul" results and produce a text file with the raw results
                        self.run_prog(self.SPROG('solppraw'),[self.in_case('simul')])
                        # Read "simul" results and produce receiver files (.vtk) of incoming and/or absorbed solar flux per-primitive
                        self.run_prog(self.SPROG('solmaps'),[self.in_case('simul')])

                        # Read "geom" and "simul" file results and produce primaries and receivers files (.vtk), and .obj geometry files
                        self.run_prog(self.SPROG('solpp'),[self.in_case('geom'),self.in_case('simul')])

                        # Read "solpaths" file and produce readable file (.vtk) by paraview to visualize the ray paths
                        self.run_prog(self.SPROG('solpaths'),[self.in_case('solpaths')])

                    efficiency_total, performance_hst=process_raw_results(self.in_case('simul'), self.casedir,rho_mirror,dni)
                    sys.stderr.write(yellow("Total efficiency: %s\n"%(repr(efficiency_total),)))
           
                    if platform.system()=="Linux":
                        os.system('mv %s/simul %s'%(self.casedir,onesunfolder))

                        if gen_vtk:
                            os.system('mv *vtk '+onesunfolder)
                            os.system('mv *obj '+onesunfolder)
                            os.system('mv *txt '+onesunfolder)
                            os.system('mv %s/geom %s'%(self.casedir,onesunfolder))
                            os.system('mv %s/solpaths %s'%(self.casedir,onesunfolder))


                    else:
                        os.system('move %s/simul %s >nul'%(self.casedir,onesunfolder))
                        if gen_vtk:
                            os.system('move *vtk %s >nul'%onesunfolder)
                            os.system('move *obj %s >nul'%onesunfolder)
                            os.system('move *txt %s >nul'%onesunfolder)
                            os.system('move %s/geom %s >nul'%(self.casedir,onesunfolder))
                            os.system('move %s/solpaths %s >nul'%(self.casedir,onesunfolder))

                ANNUAL+=performance_hst
            else:
                ANNUAL+=performance_hst


            for a in range(len(table[3:])):
                for b in range(len(table[0,3:])):
                    val=re.findall(r'\d+',    table[a+3,b+3])
                    if val==[]:
                        table[a+3,b+3]=0
                    else:
                        if c==float(val[0]):
                            table[a+3,b+3]=efficiency_total.nominal_value

            run=np.append(run,c)   
        
        annual_title=N.array(['Q_solar','Q_cosine', 'Q_shade', 'Q_hst_abs', 'Q_block', 'Q_atm', 'Q_spil', 'Q_refl', 'Q_rcv_abs']) 
        ANNUAL=N.vstack((annual_title, ANNUAL))
        np.savetxt(self.casedir+'/lookup_table.csv', table, fmt='%s', delimiter=',')
        np.savetxt(self.casedir+'/heliostats_annual_performance.csv', ANNUAL, fmt='%s', delimiter=',')
        sys.stderr.write("\n"+green("Lookup table saved.\n"))
        sys.stderr.write(green("Completed successfully.\n"+"\n"))


