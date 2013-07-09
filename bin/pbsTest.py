#!/usr/bin/env python

import sys, os, re, glob
import argparse

import subprocess


class Pbs(object):

    def __init__(self, name, queue, mail=None, nodes=None, ppn=None, depend=None,
                 out=None, err=None, workdir=None):
        self.name    = name
        self.queue   = queue
        self.nodes   = nodes
        self.ppn     = ppn
        self.mail    = mail
        self.depend  = depend
        self.workdir = workdir or os.getenv("PWD")
        self.out     = out     or os.path.join(self.workdir, name + ".o$PBS_JOBID")
        self.err     = err     or os.path.join(self.workdir, name + ".e$PBS_JOBID")
        self.cmds = []

    def addCmd(self, cmd, noop=False):
        if noop:
            cmd = "echo \""+cmd+"\""
        self.cmds.append(cmd)

    def script(self):
        s  = "#!/usr/bin/env bash\n"
        s += "#PBS -q %s\n" % (self.queue)
        s += "#PBS -N %s\n" % (self.name)
        if self.nodes and self.ppn:
            s += "#PBS -l nodes=%s:ppn=%s\n" % (str(self.nodes), str(self.ppn))
        s += "#PBS -l walltime=%d\n" % (1800)
        if self.depend:
            s += "#PBS -W depend=%s\n" % (self.depend)
        s += "#PBS -o %s\n" % (self.out)
        s += "#PBS -e %s\n" % (self.err)
        if self.mail:
            s += "#PBS -M %s\n" % (self.mail)

        for c in self.cmds:
            s += c + "\n"
        return s


    def write(self, filename):
        s = self.script()
        fp = open(filename, 'w')
        fp.write(s)
        fp.close()




#############################################################################################
#
# 
#
#############################################################################################
        
def main(db, visit, noop=False, nCcd=10, queue='batch', nodes=None, ppn=None, camera="suprimecam-mit",
         mail=None, rmlog=False):

    ###############################
    # init some variables
    
    if not os.getenv("WWW_ROOT"):
        print "Must set WWW_ROOT"
        sys.exit()
    wwwroot  = os.getenv("WWW_ROOT")
    wwwrerun = os.getenv("WWW_RERUN") or db

    # PBS logs
    cwd = os.getenv("PWD")
    logdir = os.path.join(cwd, "log")
    if not os.path.exists(logdir):
        os.mkdir(logdir)
    if rmlog:
        oldlogs = glob.glob(os.path.join(logdir, "*"))
        print "Removing old logs (%d files)" % (len(oldlogs))
        for log in oldlogs:
            os.unlink(log)

            
    # cmds for PBS script
    envFile = os.path.join(cwd, "envExports.sh")
    envSrc  = "source " + envFile
    if not os.path.exists(envFile):
        os.system("export > %s" % (envFile))
    expWwwRoot  = "export WWW_ROOT="  + wwwroot
    expWwwRerun = "export WWW_RERUN=" + wwwrerun

    # run newQa.py if the www directory is missing
    www = os.path.join(os.getenv("WWW_ROOT"), wwwrerun)
    if not os.path.exists(www):
        cmd =  "newQa.py -c green -p hsc -F %s" % (wwwrerun)
        print cmd
        qaProc = subprocess.call(cmd, shell=True)

        
    ################################
    # make the PBS jobs
    scat = Pbs("scatter", queue, nodes=nodes, ppn=ppn, mail=mail, workdir=logdir)
    gath = Pbs("gather", queue,  mail=mail, workdir=logdir)

    for s in ["cd $PBS_O_WORKDIR", envSrc, expWwwRoot, expWwwRerun]:
        scat.addCmd(s)
        gath.addCmd(s)

        
    #########
    # scatter
    scatCmd = "pipeQa.py --noWwwCache -C %s -b ccd -v %s -d db -c $PBS_ARRAYID -S none %s" % (camera,
                                                                                              visit, db)
    scat.addCmd(scatCmd, noop=noop)

    scatFile = "qsub-scat.sh"
    scat.write(scatFile)
    scatCmd = ["qsub", "-V", "-t",  "0-%d" % (nCcd-1),  scatFile]
    scatProc = subprocess.Popen(scatCmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    scatId = scatProc.communicate()[0].strip()
    print scatCmd, scatId

    
    #########
    # gather
    gath.depend = "afteranyarray:"+scatId
    gathCmd = "pipeQa.py -C %s -b ccd -v %s -d db -S summOnly %s" % (camera, visit, db)
    gath.addCmd(gathCmd, noop=noop)

    gathFile = "qsub-gath.sh"
    gath.write(gathFile)
    gathCmd = ["qsub", "-V", "%s" % (gathFile)]
    gathProc = subprocess.Popen(gathCmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    gathId = gathProc.communicate()[0].strip()
    print gathCmd, gathId
    

    
######################################################################
#
######################################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", action='version', version="0.0")
    parser.add_argument("db", type=str, help="Database name.")
    parser.add_argument("visit", type=int, help="Visit to process.")
    parser.add_argument("-c", "--nCcd", type=int, default=10, help="Number of CCDs to process")
    parser.add_argument("-C", "--camera", default="suprimecam-mit", help="camera name",
                        choices=("suprimecam-mit", "hsc", "suprimecam"))
    parser.add_argument("-M", "--mail", default="steven.bickerton@gmail.com",
                        help="email address for PBS messages")
    parser.add_argument("-n", "--noop", action='store_true', default=False,  help="Don't actually run.")
    parser.add_argument("--nodes", default=None, type=int, help="Number of nodes")
    parser.add_argument("--ppn", default=None, type=int, help="Processes per node")
    parser.add_argument("-Q", "--queue", type=str, default="batch", help="Name of PBS queue")
    parser.add_argument("-r", "--rmlog", action='store_true', default=False, help="Remove old logs.")
    parser.add_argument("-s", "--monitor", action='store_true', default=False, help="Spawn a popup monitor to watch qstat")
    args = parser.parse_args()

    main(args.db, args.visit,
         noop=args.noop, nCcd=args.nCcd, queue=args.queue,
         nodes=args.nodes, ppn=args.ppn,
         camera=args.camera, mail=args.mail, rmlog=args.rmlog)

    if args.monitor:
        import commands
        stat, pupdate = commands.getstatusoutput("which popupdate")
        if not re.search("popupdate", pupdate):
            print "You don't appear to have popupdate installed on your system.  Please ask Steve if you'd like it (it's a short python script)"
            sys.exit()
        os.system("popupdate 'qstat -t' -t 1.0 &")
        
