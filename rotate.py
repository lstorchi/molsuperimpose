import sys
import pybel

filename = ""
verbose = False

if (len(sys.argv)) == 2: 
  filename = sys.argv[1]
else:
  print "usage :", sys.argv[0] , " filename.xyz"
  exit(1)

matrix = pybel.ob.matrix3x3()
matrix.RotAboutAxisByAngle(pybel.ob.vector3(1, 0, 0), 90)

if verbose:
  for i in range(3):
      for j in range(3):
          line = "%10.5f "%(matrix.Get(i,j))
          sys.stdout.write(line)
      sys.stdout.write("\n")

imatrix = matrix.inverse()

if verbose:
  print ""
  for i in range(3):
      for j in range(3):
          line = "%10.5f "%(imatrix.Get(i,j))
          sys.stdout.write(line)
      sys.stdout.write("\n")

rotarray = pybel.ob.doubleArray(9)
matrix.GetArray(rotarray)

mol = pybel.readfile("xyz", filename).next()
mol.OBMol.Rotate(rotarray)
mol.OBMol.Translate(pybel.ob.vector3(1.0, 10.0, 3.0));
print mol.write("xyz")
