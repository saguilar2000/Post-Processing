from tracemalloc import start
import h5py
import numpy as np
import os
import sys
import struct
import matplotlib as mpl
mpl.use('Agg')

# Constants
BATE2SEC = 471300*365*24*60*60
UNITS_TIME = 471300*3.154e+7
GAMMA = 7./5.
M_H = 1.67e-24 # g
M_U = 2.38
K_B = 1.38069e-16 # cm**2 g s**-2 K**-1
CONVERSION = ((GAMMA-1)*M_H*M_U)/K_B
TIME_CONVERSION_FACTOR = 471300 * 365 * 24 * 60 * 60
DENSITY_CONVERSION_FACTOR = 1.989e+33 / (0.1 * 3.086e+18) ** 3
MASS_CONVERSION_FACTOR = 1.989e+33
COORDINATES_CONVERSION_FACTOR = 0.1 * 3.086e+18
INTERNAL_ENERGY_CONVERSION_FACTOR = (0.1 * 3.086e+18) ** 2 / ((471300 * 3.154e+7) ** 2)*CONVERSION

def copy_old(oldname,fnew):
	fold=h5py.File(oldname,'r')
	for gname in fold.keys():
		gnew=fnew.create_group(gname)
		if(gname=='Header'):
			for attr in fold[gname].attrs.keys(): 
				gnew.attrs.create(attr,fold[gname].attrs[attr])
		else:
			for dset in fold[gname].keys():
				if(dset!='KromeSpecies'): gnew.create_dataset(dset,data=fold[gname][dset])
	fold.close()

if(len(sys.argv)<=1):
	print('I need the snapshot path. Bye')
	quit()

snap_list=os.listdir(sys.argv[1])
snap_list=[s for s in snap_list if 'snapshot_' in s]
snap_list.sort()
#snap_list=snap_list[:20]
snap_list.remove('snapshot_0000.hdf5')

pfrac=0.001
if(len(sys.argv)>2):
       pfac=int(sys.argv[2])

#First prepare the density file
density_bin=open('./BinaryFiles/rho.bin','wb')
mass_bin=open('./BinaryFiles/mass.bin','wb')
coordinates_bin=open('./BinaryFiles/pos.bin','wb')
tgas_bin=open('./BinaryFiles/tgas.bin','wb')
tdust_bin=open('./BinaryFiles/tdust.bin','wb')
av_bin=open('./BinaryFiles/av.bin','wb')
for i, snap in enumerate(snap_list):
	# if "1343" not in snap:
	# 	continue
	with h5py.File(sys.argv[1]+snap,'r') as f:
		time=f['Header'].attrs['Time'].item()*TIME_CONVERSION_FACTOR
		id=f['/PartType0/ParticleIDs'][:]
		rho=f['/PartType0/Density'][:]*DENSITY_CONVERSION_FACTOR # g cm^-3
		mass=f['/PartType0/Masses'][:]*MASS_CONVERSION_FACTOR # g
		pos=f['/PartType0/Coordinates'][:]*COORDINATES_CONVERSION_FACTOR # cm
		tgas=f['/PartType0/InternalEnergy'][:]*INTERNAL_ENERGY_CONVERSION_FACTOR # K
		tdust=f['/PartType0/DustTemperature'][:] # K | No Conversion Needed
		av=-np.log(f['/PartType0/Extinction'][:]) # mag
		
		# cond_g=tgas<10. # I DON'T KNOW IF THIS IS A GOOD IDEA
		# tgas[cond_g]=10.

		id=id.ravel()
		sorted_id=np.argsort(id)

		dt=struct.pack('iifi',8,len(rho),time,8) # fixed to get time in seconds
		density_bin.write(dt)
		density_bin.write(struct.pack('i',len(rho)*4))
		dt=struct.pack('f'*len(rho),*rho[sorted_id])
		density_bin.write(dt)
		density_bin.write(struct.pack('i',len(rho)*4))

		dt=struct.pack('iifi',8,len(mass),time,8)
		mass_bin.write(dt)
		mass_bin.write(struct.pack('i',len(mass)*4))
		dt=struct.pack('f'*len(mass),*mass[sorted_id])
		mass_bin.write(dt)
		mass_bin.write(struct.pack('i',len(mass)*4))

		dt=struct.pack('iifi',8,len(pos),time,8)
		coordinates_bin.write(dt)
		coordinates_bin.write(struct.pack('i',len(pos[:,0])*4))
		dt=struct.pack('f'*len(pos[:,0]),*pos[sorted_id,0])
		coordinates_bin.write(dt)
		coordinates_bin.write(struct.pack('i',len(pos[:,0])*4))

		coordinates_bin.write(struct.pack('i',len(pos[:,1])*4))
		dt=struct.pack('f'*len(pos[:,1]),*pos[sorted_id,1])
		coordinates_bin.write(dt)
		coordinates_bin.write(struct.pack('i',len(pos[:,1])*4))

		coordinates_bin.write(struct.pack('i',len(pos[:,2])*4))
		dt=struct.pack('f'*len(pos[:,2]),*pos[sorted_id,2])
		coordinates_bin.write(dt)
		coordinates_bin.write(struct.pack('i',len(pos[:,2])*4))

		dt=struct.pack('iifi',8,len(tgas),time,8)
		tgas_bin.write(dt)
		tgas_bin.write(struct.pack('i',len(tgas)*4))
		dt=struct.pack('f'*len(tgas),*tgas[sorted_id])
		tgas_bin.write(dt)
		tgas_bin.write(struct.pack('i',len(tgas)*4))

		dt=struct.pack('iifi',8,len(tdust),time,8)
		tdust_bin.write(dt)
		tdust_bin.write(struct.pack('i',len(tdust)*4))
		dt=struct.pack('f'*len(tdust),*tdust[sorted_id])
		tdust_bin.write(dt)
		tdust_bin.write(struct.pack('i',len(tdust)*4))

		dt=struct.pack('iifi',8,len(av),time,8)
		av_bin.write(dt)
		av_bin.write(struct.pack('i',len(av)*4))
		dt=struct.pack('f'*len(av),*av[sorted_id])
		av_bin.write(dt)
		av_bin.write(struct.pack('i',len(av)*4))
		print(f"Done {snap}")
		
		del rho,mass,pos,tgas,tdust,av


density_bin.close()
mass_bin.close()
coordinates_bin.close()
tgas_bin.close()
tdust_bin.close()
av_bin.close()
print('Done Binary Files')

#Now run krome on the extracted density
print('call krome post-processing from python')
os.system('./evolve')
print('done krome')
#Finally get back the species and reconstruct a new snapshot file

out=open('./BinaryFiles/xkrome.bin','rb')
dummy,nsp,dummy=struct.unpack('iii',out.read(12))
print('Found ',nsp,' species')
#Now read the list of valid particles
dummy,npt,Time,dummy=struct.unpack('iifi',out.read(16))
for i, snap in enumerate(snap_list):
#	if(i>10): break
	new_snap=snap.split('.')[0]+'_kr.hdf5'
	print('Creating',new_snap)
	f=h5py.File('./Snapshots/'+new_snap,'w')
	#Copy all data but krome species
	copy_old(sys.argv[1]+snap,f)
	id=f['/PartType0/ParticleIDs'][:]
	id=id.ravel()
	sorted_id=np.argsort(id)

	#Now read the binary file
	Time=0.
	dummy=0
	dummy2=0
	#dummy data
	dt=out.read(16)
	if(len(dt)<16):
		break
	#Npart, Time
	dummy,npt,Time,dummy2=struct.unpack('iifi',dt)
	#read species
	dummy,=struct.unpack('i',out.read(4))
	dt=out.read(nsp*npt*4)
	species=struct.unpack('f'*nsp*npt,dt)
	dummy2,=struct.unpack('i',out.read(4))
	species=np.reshape(species,(nsp,npt)).T
	species_copy=np.zeros(shape=(npt,nsp))
	species_copy[sorted_id,:]=species
	valid_copy=np.zeros(npt,dtype='int')
	f['/PartType0'].create_dataset('KromeSpecies',data=species_copy)
	f['/PartType0'].create_dataset('KromeEvolved',data=valid_copy)
	f.close()
	del species_copy,species,dt
out.close()

print('Done!')
