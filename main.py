import sys
import numpy as np

if len(sys.argv) == 2:
	r = float(sys.argv[1])
else:
	print "This program requires", 1, "argument."
	print len(sys.argv) - 1, "given."
	sys.exit("Argument mismatch")


Atoms = ["Pb", "S"]

def coordination_num(atom, crys):

	c = 0
	dist = 0
	x, y, z = atom[1], atom[2], atom[3]
	for at in crys:
		if at != atom:
			x1, y1, z1 = at[1], at[2], at[3]
			dist = ( (x-x1)**2 + (y-y1)**2 + (z-z1)**2 ) ** 0.5
			if dist <= b+0.005:
				c += 1
#	print c
	if c == 0: print x, y, z, atom, dist
	return c







#This function removes all atoms with coordination n
def remove_lessCoordinated(crys, n):
	
	atomstoberemoved = []

	new_crys = crys

	for atom in crys:
		if coordination_num(atom, new_crys) <= n:
#			if coordination_num(atom,new_crys) == 0:
#				print atom[0], atom[1], atom[2], atom[3]
			atomstoberemoved.append(atom)

	for atom in atomstoberemoved:
		new_crys.remove(atom)

	
	return new_crys		




def write2file(crystal, filename):
	f = open(filename, 'w')

	for atom in crystal:
		x, y, z= str(atom[1]), str(atom[2]), str(atom[3])
		f.write(atom[0] + '    ' + x + '    ' + y + '    ' + z + '    ' + '\n')

	f.close()




def create_raw_crystals(crystal):

	main_crystal = []
	shifted_crystal = []
	for atom in crystal:
		x, y, z= atom[1] * b, atom[2] * b, atom[3]* b
		at = [atom[0], x, y, z]
		#shift the main co-ordinate by half the bond length
		x1, y1, z1 = x + b / 2, y, z
		at_shifted = [atom[0], x1, y1, z1]

		if ( (x**2 + y**2 + z**2)**0.5 ) <= r:
			main_crystal.append(at)
			
		if ( (x1**2 + y1**2 + z1**2)**0.5 ) <= r:
			shifted_crystal.append(at_shifted)

	return main_crystal, shifted_crystal







l = 5.9362
b = l / 2

crystal = []

#After the end of this nested loop, we have a cubic crystal of length r
for z in range( int(-r/2-1), int(r/2+1) ):
	temp1, temp2 = Atoms[0], Atoms[1]
	Atoms = [temp2, temp1]
	for x in range( int(-r/2-1), int(r/2+1) ):
		temp1, temp2 = Atoms[0], Atoms[1]
		Atoms = [temp2, temp1]
		for y in range( int(-r/2-1), int(r/2+1) ):
			crystal.append([Atoms[y%2], x, y, z])





#create spherical main crystal and a crystal with shifted origin
main_crystal, shifted_crystal = create_raw_crystals(crystal)

#creating data files

location = './Data/XYZ/'+ str(r) + 'A/'
import os
import shutil
shutil.rmtree(location)
os.mkdir(location)
write2file(main_crystal, location + 'PbS' + str(r) + '.xyz')
write2file(shifted_crystal, location + 'PbS_shifted' + str(r) + '.xyz')

##Remove atoms which have co-ordination number 3 or less

#first removing less co-ordinated atoms from main crystal
main_crystal_removed = remove_lessCoordinated(main_crystal, 3)
print 'end'
#now removing less co-ordinated atoms from the shifted crystal
shifted_crystal_removed = remove_lessCoordinated(shifted_crystal, 3)



write2file(main_crystal_removed, location + 'PbS_removed' + str(r) + '.xyz')
write2file(shifted_crystal_removed, location + 'PbS_shifted_removed' + str(r) + '.xyz')







from addNofAtoms import addNofAtoms

addNofAtoms(location + 'PbS' + str(r) + '.xyz')
addNofAtoms(location + 'PbS_shifted' +str(r) + '.xyz')
addNofAtoms(location + 'PbS_shifted_removed' + str(r) + '.xyz')
addNofAtoms(location + 'PbS_removed' + str(r) + '.xyz')











