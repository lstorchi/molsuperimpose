import re
import numpy
import math

#####################################################################

def read_ncxyz (filename, trans = True):

  filep = open(filename, "r")
  
  filep.readline()
  filep.readline()
  
  xlist = []
  ylist = []
  zlist = []
  atoms = []
  
  for line in filep:
    p = re.compile(r'\s+')
    line = p.sub(' ', line)
    line = line.lstrip()
    line = line.rstrip()
  
    plist =  line.split(" ")
  
    if (len(plist) == 4):
     atomname = plist[0]
     x = plist[1]
     y = plist[2]
     z = plist[3]
  
     xlist.append(float(x))
     ylist.append(float(y))
     zlist.append(float(z))
     atoms.append(atomname)
  
  filep.close()
  
  if trans:
    xc = numpy.mean(xlist)
    yc = numpy.mean(ylist)
    zc = numpy.mean(zlist)
  
    for i in range(len(xlist)):
      xlist[i] = xlist[i] - xc
      ylist[i] = ylist[i] - yc
      zlist[i] = zlist[i] - zc
  
  return xlist, ylist, zlist, atoms

#####################################################################

def write_ncxyz (filename, xl, yl, zl, al):

  filep = open(filename, "w")
  
  filep.write("%6d\n"%(len(al)))
  filep.write("\n")
  for i in range(0, len(al)):
      filep.write("%3s %10.5f %10.5f %10.5f\n"%(al[i], xl[i], \
              yl[i], zl[i]))
  
#####################################################################
