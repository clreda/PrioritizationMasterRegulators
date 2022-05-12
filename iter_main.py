#coding: utf-8

import subprocess as sb
from glob import glob

from params import args, root_folder, name
niterations = args["niterations"]

while (True):
    fname_ls = glob(root_folder+name+"/solutions-*.zip")
    if (len(fname_ls)!=niterations):
        sb.call("python3 main.py",shell=True)
        fname_ls = glob(root_folder+name+"/solutions-*.zip")
        for fname in fname_ls:
            sz = sb.check_output("file "+fname,shell=True).decode("utf-8").split("\n")[0]
            if ("empty" in sz):
                sb.call("rm -f "+fname,shell=True)
    else:
        break
