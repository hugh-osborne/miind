import subprocess 
import submit


dir=sys.argv[1]

build = os.path.join(directories.miind_root(),'build')
path  = os.path.join(build,'jobs',dir)

with submit.cd(build):
        subprocess.call(["ls","-l"])
        subprocess.call(['make'])

with submit.cd(build):
	f=open('joblist')
	lines=f.readlines()

for line in lines:
	name=line.strip()
	subprocess.call(['qsub','./sub.sh',name])
