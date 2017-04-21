import sys
import math
import numpy
import pybel

sys.path.append("./modules")
import kabsch_minima
import xyzutil

filename1 = ""
filename2 = ""
verbose = True
dumpalsoobmol = False

if (len(sys.argv)) == 3: 
  filename1 = sys.argv[1]
  filename2 = sys.argv[2]
else:
  print "usage :", sys.argv[0] , " target.xyz filename.xyz"
  exit(1)

xlist1, ylist1, zlist1, atoms1 = xyzutil.read_ncxyz (filename1)
xlist2, ylist2, zlist2, atoms2 = xyzutil.read_ncxyz (filename2)

if len(atoms1) == len(atoms2):
  
  dist = []
  mol1list = numpy.zeros((len(atoms1), 3))
  mol2list = numpy.zeros((len(atoms2), 3))

  for i in range(0, len(atoms2)):
      x2 = (xlist1[i] - xlist2[i])**2
      y2 = (ylist1[i] - ylist2[i])**2
      z2 = (zlist1[i] - zlist2[i])**2

      dist.append(math.sqrt(x2 + y2 + z2))

      mol1list[i, 0] = xlist1[i]
      mol1list[i, 1] = ylist1[i]
      mol1list[i, 2] = zlist1[i]

      mol2list[i, 0] = xlist2[i]
      mol2list[i, 1] = ylist2[i]
      mol2list[i, 2] = zlist2[i]

  rmsd = kabsch_minima.rmsd(mol1list, mol2list)

  count = sum(d > rmsd for d in dist)
 
  print count
  print ""
  for i in range(len(dist)):
      if dist[i] > rmsd:
          print atoms1[i], xlist1[i], ylist1[i], zlist1[i],  dist[i]

else:
  print "Wrong dim"
