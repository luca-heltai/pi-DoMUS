#!/usr/bin/env python

import sys
import shutil
import os
import glob
from sympy import *
from sympy.parsing.sympy_parser import parse_expr

def main(args=None):
    """The main routine."""

    pre_curl = ["-((y**2)/2-y)", 
                "-t*((y**2)/2-y)", 
                "-y**2*(y/3 - 1/2)", 
                "-t*y**2*(y/3 - 1/2)",
                "(1/pi)*cos(pi*y)", 
                "sin(2*pi*t)*(1/pi)*cos(pi*y)",
                "(1/(2*pi))*sin(2*pi*x)*sin(2*pi*y)", 
                "(1/(2*pi))*sin(2*pi*t)*sin(2*pi*x)*sin(2*pi*y)"]   
    testnames = ["patch_test",
                 "time_dep_patch_test",
                 "rotating_patch_test",
                 "time_dep_rotating_patch_test",
                 "rotating_sine_test",
                 "time_dep_rotating_sine_test",
                 "sine_cosine_test",
                 "time_dep_sine_cosine_test"]
    press = "0"

    if args is None:
        args = sys.argv[1:]
    if len(args) == 0:
        err = "you need to specify your test:\nALL: perform all tests\n"
        for it in range(0,len(testnames)):
            err += str(it) + "  : " + testnames[it] +"\n"
        sys.exit(err)
    else:
        i = args[0]
        if i == "ALL":
            start = 0
            end = len(pre_curl)
        else:
            start = int(i)
            end = int(i)+1

    x, y, z, t = var('x,y,z,t')
    local_dict = {"x": x, "y": y, "z": z, "t": t}
    for i in range(start,end):

        pre_curl[i] = parse_expr(pre_curl[i], local_dict=local_dict)
        p = parse_expr(press, local_dict=local_dict)
        pt = p.diff(t)

        d = [0*p,0*p]
        dt = [di.diff(t) for di in d]

        v = [-pre_curl[i].diff(y), pre_curl[i].diff(x)]
        vt = [vi.diff(t) for vi in v]

        print("\nPerforming "+ testnames[i] +
              "\nUsing pattern d, d, u, u, p in \n"+
              "             (" + str(d[0]) + ", " + str(d[1]) + ", " + str(v[0]) + ", " + str(v[1]) + ", " + str(p)+" )\n")

        f_navier = [0, 0, v[0].diff(t) + (v[0])*v[0].diff(x) + (v[1])*v[0].diff(y), v[1].diff(t) + (v[0])*v[1].diff(x) + (v[1])*v[1].diff(y)]
        f_stokes = [0, 0, -v[0].diff(x, 2) - v[0].diff(y, 2) + p.diff(x), -v[1].diff(x, 2) - v[1].diff(y, 2) + p.diff(y)]

        f = []
        for j in range(len(f_navier)):
            f += [f_navier[j] + f_stokes[j]] 

        v_str = ("subsection Dirichlet boundary conditions\n"   
                +"  set IDs and component masks = 0=d;u\n"   
                +"  set IDs and expressions     = 0=" 
                + str(d[0]) + "; " + str(d[1]) + "; " + str(v[0]) + "; " + str(v[1]) + "; " + str(p) + "\n" 
                +"  set Known component names   = d,d,u,u,p\n"
                +"  set Used constants          = \nend\n" +
                "subsection Time derivative of Dirichlet boundary conditions\n"   
                +"  set IDs and component masks = 0=d;u\n"   
                +"  set IDs and expressions     = 0=" 
                + str(dt[0]) + "; " + str(dt[1]) + "; " + str(vt[0]) + "; " + str(vt[1]) + "; " + str(pt) + "\n"
                +"  set Known component names   = d,d,u,u,p\n"
                +"  set Used constants          = \nend\n" +
                "subsection Exact solution\n"
                +"  set Function expression = " 
                + str(d[0]) + "; " + str(d[1]) + "; " + str(v[0]) + "; " + str(v[1]) + "; " + str(p) + "\nend\n" +
                "subsection Forcing terms\n"   
                +"  set IDs and component masks = 0=ALL\n"   
                +"  set IDs and expressions     = 0=" 
                + str(f[0]) + "; " + str(f[1]) + "; " + str(f[2]) + "; " + str(f[3]) + "; " + "0\n"
                +"  set Known component names   = d,d,u,u,p\n  set Used constants          = \nend\n" +
                "subsection Initial solution\n"   
                +"  set Function expression ="   
                + str(d[0].subs(t,0)) + "; " + str(d[1].subs(t,0)) + "; " 
                + str(v[0].subs(t,0)) + "; " + str(v[1].subs(t,0)) + "; " + str(p.subs(t,0)) + "\nend\n" +
                "subsection Initial solution_dot\n"   
                +"  set Function expression ="   
                + str(dt[0].subs(t,0)) + "; " + str(dt[1].subs(t,0)) + "; " 
                + str(vt[0].subs(t,0)) + "; " + str(vt[1].subs(t,0)) + "; " + str(pt.subs(t,0)) +  
                "\nend")

        shutil.copy2('./template.prm', 'ALE_'+ testnames[i]+ '.prm' )

        prm = open('ALE_'+ testnames[i]+ '.prm','a+')
        prm.write(v_str.replace("**", "^"))
        prm.close()

        os.system("sed -i 's/xxx/" + testnames[i] +"/g' ALE_"+ testnames[i]+ '.prm' )
        os.system("sed -i 's/yyy/" + testnames[i] +"/g' ALE_"+ testnames[i]+ '.prm' )

        os.system("mpirun -np 4 ../build/heart --prm=ALE_"+ testnames[i] +".prm")
        os.system("mv err_"+testnames[i]+".txt err"+str(i)+"_"+testnames[i]+".txt")
    files = glob.glob( 'err*' )
    files.sort()
    with open( 'result.txt', 'w' ) as result:
        for file in files:
            result.write("\nerrortable: "+file[5:]+"\n")
            for line in open( file, 'r' ):
                result.write( line )

    os.system("rm -rf err*")
if __name__ == "__main__":
    main()
