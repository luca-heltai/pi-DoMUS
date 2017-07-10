#!/usr/bin/env python

import sys
from sympy import *

from sympy.parsing.sympy_parser import parse_expr

def main(args=None):
    """The main routine."""
    if args is None:
        args = sys.argv[1:]

    pre_curl = "sin(2*pi*x)*sin(2*pi*y)"
    p = "0"

    if len(args) > 0:
        pre_curl = args[0]
    else:
        pre_curl = "t*sin(2*pi*x)*sin(2*pi*y)"

    if len(args) > 1:
        p = args[1]
    else:
        p = "0"

    print("\n\nUsing v = curl(" + pre_curl + " e_z), p = " +p +"\n\n")

    x, y, z, t = var('x,y,z,t')
    local_dict = {"x": x, "y": y, "z": z, "t": t}

    pre_curl = parse_expr(pre_curl, local_dict=local_dict)
    p = parse_expr(p, local_dict=local_dict)

    v = [0, 0, -pre_curl.diff(y), pre_curl.diff(x)]
    vt = [vi.diff(t) for vi in v]

    f_navier = [0, 0, v[2].diff(t) + (v[2])*v[2].diff(x) + (v[3])*v[2].diff(y), v[3].diff(t) + (v[2])*v[3].diff(x) + (v[3])*v[3].diff(y)]
    f_stokes = [0, 0, -v[2].diff(x, 2) - v[2].diff(y, 2) + p.diff(x), -v[3].diff(x, 2) - v[3].diff(y, 2) + p.diff(y)]

    f = []
    for i in range(len(f_navier)):
        f += [f_navier[0] + f_stokes[i]] 
    v_str = "subsection Exact solution\n  set Function expression = " + str(v[0]) + "; " + str(v[1]) + "; " + str(v[2]) + "; " + str(v[3]) + "; " + str(p) + "\nend"
    print(v_str.replace("**", "^"))

    v_str = "subsection Forcing terms\n  set Function expression = " + str(f[0]) + "; " + str(f[1]) + "; " + str(f[2]) + "; " + str(f[3]) + "; 0\nend"
    print(v_str.replace("**", "^"))

if __name__ == "__main__":
    main()

