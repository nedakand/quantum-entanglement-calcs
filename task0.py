import numpy as np

import math as m



# returns photon wave function linearly polarized at angle theta

def photon(theta):
  a = theta/180. * m.pi
  return np.matrix([[m.cos(a)],[m.sin(a)]])



# returns operator (2x2 matrix) for linear polarizer at angle theta

# use numpy matrixes instead of arrays - then multiplication is by matrix rules, not by element

def pol(theta):
   a = theta/180.*m.pi
   return np.matrix([[m.cos(a)**2,m.cos(a)*m.sin(a)],

[m.cos(a)*m.sin(a),m.sin(a)**2]])



# calculates integral <phi1 | operator | phi2> (expectation value)

def e(phi1,operator,phi2):
   return (phi1.transpose()*(operator * phi2))[0,0]



# calculates expectation value for product of two channel coincidence counts with

# polarizer at angle a in channel 1 and b in channel 2

# returns an array [e1, e2, e3]: (not matrix but numpy array, i.e. mult and div will work by element)

# e1 = expectation value for non-entangled correlated pairs of VV and HH

# e2 = expectation value for entangled symmetric state, VV+HH

# e3 = expectation value for entangled antisymmetric state, VV-HH

def C(a,b):
   v = photon(0) # vertically polarized photon
   h = photon(90) # horizontallly polarised photon
   # correlated but non-entangled case
   ecorr = (e(v,pol(a),v)*e(v,pol(b),v) + e(h,pol(a),h)*e(h,pol(b),h) )/2
   # entangled cross term
   entang = (e(h,pol(a),v)*e(h,pol(b),v) + e(v,pol(a),h)*e(v,pol(b),h))/2
   return np.array([ecorr, ecorr + entang, ecorr - entang])



# Returns expectation value for 4 measurements at angles alpha and betta

# in two channels (plus 90 degrees)

# returns array [E1, E2, E3] for non-entangled (E1) and entangled + case (E2) and – (E3)

# (note that C(a,b) are arrays of size 3)

def E(a,b):
   return (C(a,b) - C(a,b+90) - C(a+90,b) + C(a+90,b+90)) / (C(a,b) + C(a,b+90) + C(a+90,b) + C(a+90,b+90))



# prints array with given precision, sets small values to zeros

# just to avoid printing ~1e-20 numbers, these are just zeros

def print2(string,a,crlf = True):
   print(string,end = '')
   for i in range(len(a)):
      if i!=0: print(', ',end='')
      print(f'{a[i]:.3f}',end = '')
   if crlf: print()
   return


# MAIN PROGRAM

#calculate some c(a,b)

print ('All pairs below are for [non-entangled, entangled +, entangled -] case')

for angles in [[90,90],[0,0],[90,0],[0,90],[45,45],[45,-45]]:
   print2(f'Signal at {angles[0]:4d},{angles[1]:4d} = ',C(angles[0],angles[1]))
   print(" *** DOLP for some angles")

for angle in range(0,91,15):
   par = C(angle,angle)
   a = angle + 90
   if a>90: a=a-180
   per = C(angle,a)
   print(f'DOLP {angle:3d}:')
   print2(' parallel = ',par)
   print2(' perpendi = ',per)
   print(f' DOLP = {(par-per)/(par+per)}')



print(' Bell Inequality:')

a = [[0,22.5],[45,22.5],[0,67.5],[45,67.5]]

S = np.zeros(3)

for i in range(len(a)):
   ee = E(a[i][0],a[i][1])
   print2(f' E{i} = ',ee)
   if i == 2: S = S - ee
   else: S = S + ee
   print2(' S = ',S)
   #print2(' S entangled symm = ', )

