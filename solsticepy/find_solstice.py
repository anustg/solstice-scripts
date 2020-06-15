# This script finds the location of Solstice from the Windows Registry, assuming
# it was installed with ANU's NSIS-based installer. It works with Python 2.7 and
# Python 3.8, with both 32 and 64-bit versions, on both Windows 7 and Windows 10
# -- John Pye, April 2020

import sys, platform, os

def find_solstice_root(version_required=None,verbose=0):
	"""Locate the place where Solstice files are installed on Windows

	``Arguments``

	  * version_required (None or str): if not None, then enforce this specified version of Solstice (eg '0.9.0')
	  * verbose (bool): if True, output file locations to stderr

	``Return``

	  * dirn (str): the directory where Solstice is installed in the Windows system
	"""

	if verbose: sys.stderr.write("Python is running from %s\n"%(sys.executable,))
	if not platform.system()=="Windows":
		raise RuntimeError("This function is only for Windows")
	# otherwise...
	if sys.version_info[0] < 3:
		if verbose: sys.stderr.write("Python 2, ")
		import _winreg as winreg
	else:
		if verbose: sys.stderr.write("Python 3, ")
		import winreg
	if sys.maxsize > 2**32:
		if verbose: sys.stderr.write("64-bit\n")
		# 64-bit Python
		key = winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE,r"Software\Meso-Star\Solstice",0,winreg.KEY_READ)
	else:
		if verbose: sys.stderr.write("32-bit\n")
		key = winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE,r"Software\Meso-Star\Solstice",0,winreg.KEY_WOW64_64KEY+winreg.KEY_READ)
	dirn, typd = winreg.QueryValueEx(key, 'root')
	ver, typv = winreg.QueryValueEx(key, 'version')
	if version_required:
		assert ver == version_required
	return dirn

def find_prog(name,version_required=None):
	"""Find the path to any required Solstice executable program

	``Arguments``

	  * name (str): stem-name (eg 'solpp') of the program (eg 'solpp.exe') required
	  * version_required (None or str): if not None, then enforce this specified version of Solstice (e.g. '0.9.0')

	``Return``

	  * path (str): path of the required Solstice executable program
	"""
	
	if platform.system()=="Windows":
		path = os.path.join(find_solstice_root(version_required),"bin","%s.exe"%(name,));
		if not os.path.exists(path):
			raise RuntimeError("Program '%s' was not found at path '%s'"%(name,path))
		return path
	else:
		# assume all solstice programs are on the PATH
		import subprocess
		DEVNULL = open(os.devnull, 'wb')
		rc = subprocess.call(['which',name],stdout=DEVNULL)
		DEVNULL.close()
		if rc:
			raise RuntimeError("Program '%s' was not found in the PATH" %(name))
		return name

if __name__=="__main__":
	dirn = find_solstice_root('0.9.0',verbose=1)
	sys.stderr.write("Solstice is installed in %s\n\n" %(dirn,)) # works both python2+3 :o)

	import subprocess, os

	spath = os.path.join(dirn,"bin","solstice.exe")

	# check that installed version of solstice is 0.9.0:
	ret = subprocess.check_output([spath,"--version"])
	assert ret.decode('ascii').strip() == "Solstice 0.9.0"

	# output the help text from solstice:
	subprocess.check_call([spath,"-h"])

# you should see the 'solstice -help' output shown in your console.

# checked with Python 2.7 32-bit on Windows 7 
# checked with Python 3.8.2 32-bit on Windows 10
# checked with Python 3.8.2 64-bit on Windows 10

