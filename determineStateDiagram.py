
import sys
import numpy as np

# This script determines the boundaries between the different region in the state diagram. 
# The lines of the resulting data file contain tauD as the first value and then the values of nuP that mark the boundaries.

# coefficients
tauR = 1.0
Dtr = 1.0
gamma = 1.0
visc = 1.0
nuV = 1.0
elastMod = 1.0
nuU = 10.0
kappa = 0.0
ellOff = 0.0

tauDStart = 0.01
tauDEnd = 100.0
numtauD = 1000
tauDList = np.logspace(np.log10(tauDStart), np.log10(tauDEnd),num=numtauD, endpoint=True)

# tolerances	
grRateTol = 0.0000001
rotSolTol = 0.00000001
tol = 0.001

# sample wavenumbers 		
kZero = 0.0		
kMax = 5.0
numk = 100
kList = np.linspace(kZero,kMax,num=numk,endpoint=False)

# data file to save data
dataFileName = 'dataStateDiagram_gamma_%05.2f'%(gamma)+'_visc_%05.2f'%(visc)+'_frict_%05.2f'%(nuV)+'_elastMod_%05.2f'%(elastMod)+'_rest_%05.2f'%(nuU)+'_kappa_%05.2f'%(kappa)+'.dat'
dataFileName = dataFileName.replace('.','d',6)
dataFile = open(dataFileName,"w")
dataFile.close()
dataFile = open(dataFileName,"a",buffering=9999)	

# identity matrix
identMat = 1j*np.zeros((2,2))
identMat[0,0] = 1.0
identMat[1,1] = 1.0

# move through values of tauD
for tauD in tauDList:
	
	print('tauD = ' + str(tauD))
	
	
	# nuP above which stationary solution starts to exist
	nuPPol = 3.0*(nuV + nuU*tauD)/(2.0*gamma*tauR)
	dataFile.write('%e'%(tauD)+' %e'%(nuPPol))
	
	
	# find nuP where rotational solutions start to exist
	nuPLow = 0.1
	nuPHigh = 200.0	
	while (nuPHigh - nuPLow) > tol:
		nuP = 0.5*(nuPHigh + nuPLow) + 0.0*1j
		A = np.sqrt(16.0*gamma**2*nuP**2*tauD**2*tauR**2 + 8.0*gamma*nuP*tauD**2*tauR*(2.0*gamma*nuP*tauR - 3.0*nuV) - 96.0*gamma*nuP*tauD*tauR**2*(nuU*tauD + nuV) + tauD**2*(2.0*gamma*nuP*tauR - 3.0*nuV)**2 + 12.0*tauD*tauR*(nuU*tauD + nuV)*(2.0*gamma*nuP*tauR - 3.0*nuV) + 36.0*tauR**2*(nuU*tauD + nuV)**2)
		P0Plus = np.sqrt(2.0)*np.sqrt(3.0*nuV*tauD + 12.0*nuV*tauR - 6.0*tauD*tauR*(gamma*nuP - 2.0*nuU) + A)*np.sqrt(1.0/(-3.0*nuV*tauD + 6.0*nuV*tauR + 6.0*tauD*tauR*(gamma*nuP + nuU) - A))
		u0Plus = np.sqrt(6.0)*np.sqrt(nuP)*np.sqrt(3.0*nuV*tauD + 12.0*nuV*tauR - 6.0*tauD*tauR*(gamma*nuP - 2.0*nuU) + A)/(6.0*np.sqrt(gamma)*np.sqrt(nuU)*np.sqrt(tauR)*np.sqrt(nuU*tauD + nuV))
		phi0Plus = np.arccos(np.sqrt(3.0*nuU*tauD + 3.0*nuV)/(6.0*np.sqrt(gamma)*np.sqrt(nuP)*np.sqrt(nuU)*tauD*np.sqrt(tauR)*np.sqrt(1.0/(-3.0*nuV*tauD + 6.0*nuV*tauR + 6.0*tauD*tauR*(gamma*nuP + nuU) - A))))
		omega0Plus = np.sqrt(nuU*tauD + nuV)*np.sqrt(12.0*gamma*nuP*nuU*tauD**2*tauR - (nuU*tauD + nuV)*(-3.0*nuV*tauD + 6.0*nuV*tauR + 6.0*tauD*tauR*(gamma*nuP + nuU) - A))*np.sqrt(1.0/(-3.0*nuV*tauD + 6.0*nuV*tauR + 6.0*tauD*tauR*(gamma*nuP + nuU) - A))/(nuV*tauD)
		
		# check whether solution is real-valued
		if abs(np.imag(P0Plus)) < rotSolTol and abs(np.imag(u0Plus)) < rotSolTol and abs(np.imag(phi0Plus)) < rotSolTol and abs(np.imag(omega0Plus)) < rotSolTol and abs(np.real(P0Plus)) > rotSolTol and abs(np.real(u0Plus)) > rotSolTol and abs(np.real(phi0Plus)) > rotSolTol and abs(np.real(omega0Plus)) > rotSolTol:
			nuPHigh = np.real(nuP)			
		else:	
			nuPLow = np.real(nuP)	
	
	nuPOscStart = np.real(0.5*(nuPHigh + nuPLow))
	dataFile.write(' %e'%(nuPOscStart))
	
	
	# find nuP below which hysteresis might emerge due to subcritical transition
	nuPLow = nuPOscStart - tol
	nuPHigh = 100.0
	while (nuPHigh - nuPLow) > tol:
		nuP = 0.5*(nuPHigh + nuPLow) + 0.0*1j
		A = -np.sqrt(16.0*gamma**2*nuP**2*tauD**2*tauR**2 + 8.0*gamma*nuP*tauD**2*tauR*(2.0*gamma*nuP*tauR - 3.0*nuV) - 96.0*gamma*nuP*tauD*tauR**2*(nuU*tauD + nuV) + tauD**2*(2.0*gamma*nuP*tauR - 3.0*nuV)**2 + 12.0*tauD*tauR*(nuU*tauD + nuV)*(2.0*gamma*nuP*tauR - 3.0*nuV) + 36.0*tauR**2*(nuU*tauD + nuV)**2)
		P0Minus = np.sqrt(2.0)*np.sqrt(3.0*nuV*tauD + 12.0*nuV*tauR - 6.0*tauD*tauR*(gamma*nuP - 2.0*nuU) + A)*np.sqrt(1.0/(-3.0*nuV*tauD + 6.0*nuV*tauR + 6.0*tauD*tauR*(gamma*nuP + nuU) - A))
		u0Minus = np.sqrt(6.0)*np.sqrt(nuP)*np.sqrt(3.0*nuV*tauD + 12.0*nuV*tauR - 6.0*tauD*tauR*(gamma*nuP - 2.0*nuU) + A)/(6.0*np.sqrt(gamma)*np.sqrt(nuU)*np.sqrt(tauR)*np.sqrt(nuU*tauD + nuV))
		phi0Minus = np.arccos(np.sqrt(3.0*nuU*tauD + 3.0*nuV)/(6.0*np.sqrt(gamma)*np.sqrt(nuP)*np.sqrt(nuU)*tauD*np.sqrt(tauR)*np.sqrt(1.0/(-3.0*nuV*tauD + 6.0*nuV*tauR + 6.0*tauD*tauR*(gamma*nuP + nuU) - A))))
		omega0Minus = np.sqrt(nuU*tauD + nuV)*np.sqrt(12.0*gamma*nuP*nuU*tauD**2*tauR - (nuU*tauD + nuV)*(-3.0*nuV*tauD + 6.0*nuV*tauR + 6.0*tauD*tauR*(gamma*nuP + nuU) - A))*np.sqrt(1.0/(-3.0*nuV*tauD + 6.0*nuV*tauR + 6.0*tauD*tauR*(gamma*nuP + nuU) - A))/(nuV*tauD)
		
		# check whether solution is real-valued
		if abs(np.imag(P0Minus)) < rotSolTol and abs(np.imag(u0Minus)) < rotSolTol and abs(np.imag(phi0Minus)) < rotSolTol and abs(np.imag(omega0Minus)) < rotSolTol and abs(np.real(P0Minus)) > rotSolTol and abs(np.real(u0Minus)) > rotSolTol and abs(np.real(phi0Minus)) > rotSolTol and abs(np.real(omega0Minus)) > rotSolTol:
			nuPLow = np.real(nuP)	
		else:	
			nuPHigh = np.real(nuP)	
			
	nuPOscHyst = np.real(0.5*(nuPHigh + nuPLow))
	
	dataFile.write(' %e'%(nuPOscHyst))	
		
	
	# determine the onset of finite-wavelength instability
	# only outside the region of rotational states (excluding the area of hysteresis)	
	if nuPPol > nuPOscHyst:
		dataFile.write(' nan nan')
		
	else:
		nuPLow = nuPPol + tol
		nuPHigh = nuPOscHyst
		
		while (nuPHigh - nuPLow) > tol:
						
			nuP = 0.5*(nuPHigh + nuPLow) + 0.0*1j
			
			P0 = np.sqrt(1.0 - 3.0*(nuV + nuU*tauD)/(2.0*gamma*tauR*nuP))
			v0 = nuP*P0/(nuV + nuU*tauD)
			u0 = v0*tauD
			
			eigenvalsReZeroMax = 0.0
			eigenvalsReMax = 0.0
			eigenvalsImMaxAtMaxGR = 0.0
			eigenvalsImMaxAtAnyPosGR = 0.0
			
			# perform linear stability analysis for sample values of the wavenumber
			for k in kList:
			
				kx = k
				ky = 0.0

				kSq = kx*kx + ky*ky
				if kSq == 0.0:
					kSq = 0.000000001
	
				MPP = 1j*np.zeros((2,2))
				MPP[0,0] = - v0*1j*kx - Dtr*kSq - 1/tauR - 4*gamma*P0*v0/3
				MPP[1,1] = - v0*1j*kx - Dtr*kSq - 1/tauR - gamma*P0*v0

				MPv = 1j*np.zeros((2,2))
				MPv[0,0] = kappa*1j*kx*P0 + 2*gamma/3 - 2*gamma*P0*P0/3
				MPv[1,0] = - 1j*ky*P0/2 + kappa*1j*ky*P0/2
				MPv[1,1] = 1j*kx*P0/2 + kappa*1j*kx*P0/2 + 2*gamma/3 + gamma*P0*P0/3

				MvP = 1j*np.zeros((2,2))
				MvP[0,0] = nuP
				MvP[1,1] = nuP

				Muv = 1j*np.zeros((2,2))
				Muv[0,0] = 1
				Muv[1,1] = 1 

				proj = 1j*np.zeros((2,2))
				proj[0,0] = 1 - kx*kx/kSq
				proj[0,1] = - kx*ky/kSq
				proj[1,0] = - kx*ky/kSq
				proj[1,1] = 1 - ky*ky/kSq

				identMat = 1j*np.zeros((2,2))
				identMat[0,0] = 1
				identMat[1,1] = 1

				M = 1j*np.zeros((4,4))

				MPartPP = MPP + MPv@proj@MvP/(visc*kSq + nuV)

				MPartPu = - (elastMod*kSq + nuU)*MPv/(visc*kSq + nuV)

				MPartuP = nuP*Muv@proj/(visc*kSq + nuV)

				MPartuu = - (v0*1j*kx + 1/tauD)*identMat - (elastMod*kSq + nuU)*proj@Muv/(visc*kSq + nuV)

				M[0,0] = MPartPP[0,0]
				M[0,1] = MPartPP[0,1]
				M[1,0] = MPartPP[1,0]
				M[1,1] = MPartPP[1,1]
					
				M[0,2] = MPartPu[0,0]
				M[0,3] = MPartPu[0,1]
				M[1,2] = MPartPu[1,0]
				M[1,3] = MPartPu[1,1]

				M[2,0] = MPartuP[0,0]
				M[2,1] = MPartuP[0,1]
				M[3,0] = MPartuP[1,0]
				M[3,1] = MPartuP[1,1]

				M[2,2] = MPartuu[0,0]
				M[2,3] = MPartuu[0,1]
				M[3,2] = MPartuu[1,0]
				M[3,3] = MPartuu[1,1]
				
				# calculate eigenvalues and -vectors
				eigenvalues, eigenvectors = np.linalg.eig(M)
				
				# largest growth rate at zero wavenumber
				if k == kZero:
					for i in range(0,4):
						if eigenvalsReZeroMax < np.real(eigenvalues[i]):
							eigenvalsReZeroMax = np.real(eigenvalues[i])
				# largest growth rate for any wavenumber
				for i in range(0,4):
					if eigenvalsReMax < np.real(eigenvalues[i]):
						eigenvalsReMax = np.real(eigenvalues[i])
			
			if eigenvalsReMax > grRateTol and eigenvalsReMax > eigenvalsReZeroMax:
				nuPHigh = nuP
			else:
				nuPLow = nuP
	
		nuPFiniteWavelength = 0.5*(nuPHigh + nuPLow)
		
		# outside the region of rotational states
		if nuPFiniteWavelength < nuPOscStart: 
			dataFile.write(' %e'%(nuPFiniteWavelength) + ' %e'%(nuPOscStart))
		# inside the region of rotational states
		else:
			dataFile.write(' nan %e'%(nuPFiniteWavelength))
	
	
	# nuP where stationary polar state starts to exist within the region of rotational states
	if nuPPol < nuPOscStart:
		dataFile.write(' %e'%(nuPOscStart))
	else:
		dataFile.write(' %e'%(nuPPol))
	
	dataFile.write('\n')
		
dataFile.close()

