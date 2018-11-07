import scipy as sp
import scipy.interpolate
import numpy as np

from units import unit

with open("fit-aoadata/output_filt_cl.csv") as cl:
    cl_aoadata = cl.read()
with open("fit-aoadata/output_filt_cd.csv") as cd:
    cd_aoadata = cd.read()

cl_aoa = []
cd_aoa = []
cl_aoadata = cl_aoadata.split("\n")[:-1]
cd_aoadata = cd_aoadata.split("\n")[:-1]
for n in range(len(cl_aoadata)):
    cl_aoa.append([float(i) for i in cl_aoadata[n].split(',')])
for n in range(len(cd_aoadata)):
    cd_aoa.append([float(i) for i in cd_aoadata[n].split(',')])

cl_interp = sp.interpolate.interp1d(x=[i[0] for i in cl_aoa], y=[i[1] for i in cl_aoa], kind="linear")
cd_interp = sp.interpolate.interp1d(x=[i[0] for i in cd_aoa], y=[i[1] for i in cd_aoa], kind="linear")

def getCl(aoa: float):
    if type(aoa) == type(1 * unit.degree) or type(aoa) == type(1*unit.radian):
        aoa = aoa.m
    aoa = abs(aoa)
    return cl_interp(aoa)
def getCd(aoa: float):
    if type(aoa) == type(1 * unit.degree) or type(aoa) == type(1*unit.radian):
        aoa = aoa.m
    aoa = abs(aoa)
    return cd_interp(aoa)
