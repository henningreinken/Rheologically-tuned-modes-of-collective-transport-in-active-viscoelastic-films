
import os
import numpy as np
import sys
import time

# coefficients
tauR = 1.0
Dtr = 1.0
gam = 1.0
visc = 1.0
nuV = 1.0
elastMod = 1.0
nuU = 10.0
kappa = 0.0
tauD = 0.52062

N = 64
tEnd = 100.0

dt = 0.001

runNum = 0

nStep = 0

numSteps = round(tEnd/dt) 

dtPrint = 1.0
dnPrint = round(dtPrint/dt)

# get current directory
thisPath = os.getcwd()

# where to save the results of the search
dataFileName = thisPath + '/stripesHyst_tauD_%07.5f'%(tauD)+'.dat'
dataFileName = dataFileName.replace('.','d',1)
dataFile = open(dataFileName,"w")
dataFile.close()
dataFile = open(dataFileName,"a")



# sample wavenumbers to determine finite-wavelength instability
kZero = 0.0		
kMax = 5.0
numk = 100
kList = np.linspace(kZero,kMax,num=numk,endpoint=False)

# tolerances
tol = 0.0001
grRateTol = 0.0000001
rotSolTol = 0.00000001

# nuP above which the stationary polar solution exists 
nuPPol = 3.0*(nuV + nuU*tauD)/(2.0*gam*tauR)
	
# determine nuP above which rotational solutions exist	
nuPLow = 0.1
nuPHigh = 200.0	
while (nuPHigh - nuPLow) > tol:
	
	nuP = 0.5*(nuPHigh + nuPLow) + 0.0*1j
	
	A = np.sqrt(16.0*gam**2*nuP**2*tauD**2*tauR**2 + 8.0*gam*nuP*tauD**2*tauR*(2.0*gam*nuP*tauR - 3.0*nuV) - 96.0*gam*nuP*tauD*tauR**2*(nuU*tauD + nuV) + tauD**2*(2.0*gam*nuP*tauR - 3.0*nuV)**2 + 12.0*tauD*tauR*(nuU*tauD + nuV)*(2.0*gam*nuP*tauR - 3.0*nuV) + 36.0*tauR**2*(nuU*tauD + nuV)**2)
	P0Plus = np.sqrt(2.0)*np.sqrt(3.0*nuV*tauD + 12.0*nuV*tauR - 6.0*tauD*tauR*(gam*nuP - 2.0*nuU) + A)*np.sqrt(1.0/(-3.0*nuV*tauD + 6.0*nuV*tauR + 6.0*tauD*tauR*(gam*nuP + nuU) - A))
	u0Plus = np.sqrt(6.0)*np.sqrt(nuP)*np.sqrt(3.0*nuV*tauD + 12.0*nuV*tauR - 6.0*tauD*tauR*(gam*nuP - 2.0*nuU) + A)/(6.0*np.sqrt(gam)*np.sqrt(nuU)*np.sqrt(tauR)*np.sqrt(nuU*tauD + nuV))
	phi0Plus = np.arccos(np.sqrt(3.0*nuU*tauD + 3.0*nuV)/(6.0*np.sqrt(gam)*np.sqrt(nuP)*np.sqrt(nuU)*tauD*np.sqrt(tauR)*np.sqrt(1.0/(-3.0*nuV*tauD + 6.0*nuV*tauR + 6.0*tauD*tauR*(gam*nuP + nuU) - A))))
	omega0Plus = np.sqrt(nuU*tauD + nuV)*np.sqrt(12.0*gam*nuP*nuU*tauD**2*tauR - (nuU*tauD + nuV)*(-3.0*nuV*tauD + 6.0*nuV*tauR + 6.0*tauD*tauR*(gam*nuP + nuU) - A))*np.sqrt(1.0/(-3.0*nuV*tauD + 6.0*nuV*tauR + 6.0*tauD*tauR*(gam*nuP + nuU) - A))/(nuV*tauD)
		
	if abs(np.imag(P0Plus)) < rotSolTol and abs(np.imag(u0Plus)) < rotSolTol and abs(np.imag(phi0Plus)) < rotSolTol and abs(np.imag(omega0Plus)) < rotSolTol and abs(np.real(P0Plus)) > rotSolTol and abs(np.real(u0Plus)) > rotSolTol and abs(np.real(phi0Plus)) > rotSolTol and abs(np.real(omega0Plus)) > rotSolTol:
		nuPHigh = np.real(nuP)
	else:	
		nuPLow = np.real(nuP)

nuPOscStart = np.real(0.5*(nuPHigh + nuPLow))

print('nuP above which rotational solutions exist: ' + str(nuPOscStart))


# determine nuP above which rotational solutions exist and uniform stationary states are unstable (upper boundary of hysteresis region)	
nuPLow = nuPOscStart - tol
nuPHigh = 100.0
while (nuPHigh - nuPLow) > tol:
	
	nuP = 0.5*(nuPHigh + nuPLow) + 0.0*1j
	
	A = -np.sqrt(16.0*gam**2*nuP**2*tauD**2*tauR**2 + 8.0*gam*nuP*tauD**2*tauR*(2.0*gam*nuP*tauR - 3.0*nuV) - 96.0*gam*nuP*tauD*tauR**2*(nuU*tauD + nuV) + tauD**2*(2.0*gam*nuP*tauR - 3.0*nuV)**2 + 12.0*tauD*tauR*(nuU*tauD + nuV)*(2.0*gam*nuP*tauR - 3.0*nuV) + 36.0*tauR**2*(nuU*tauD + nuV)**2)
	P0Minus = np.sqrt(2.0)*np.sqrt(3.0*nuV*tauD + 12.0*nuV*tauR - 6.0*tauD*tauR*(gam*nuP - 2.0*nuU) + A)*np.sqrt(1.0/(-3.0*nuV*tauD + 6.0*nuV*tauR + 6.0*tauD*tauR*(gam*nuP + nuU) - A))
	u0Minus = np.sqrt(6.0)*np.sqrt(nuP)*np.sqrt(3.0*nuV*tauD + 12.0*nuV*tauR - 6.0*tauD*tauR*(gam*nuP - 2.0*nuU) + A)/(6.0*np.sqrt(gam)*np.sqrt(nuU)*np.sqrt(tauR)*np.sqrt(nuU*tauD + nuV))
	phi0Minus = np.arccos(np.sqrt(3.0*nuU*tauD + 3.0*nuV)/(6.0*np.sqrt(gam)*np.sqrt(nuP)*np.sqrt(nuU)*tauD*np.sqrt(tauR)*np.sqrt(1.0/(-3.0*nuV*tauD + 6.0*nuV*tauR + 6.0*tauD*tauR*(gam*nuP + nuU) - A))))
	omega0Minus = np.sqrt(nuU*tauD + nuV)*np.sqrt(12.0*gam*nuP*nuU*tauD**2*tauR - (nuU*tauD + nuV)*(-3.0*nuV*tauD + 6.0*nuV*tauR + 6.0*tauD*tauR*(gam*nuP + nuU) - A))*np.sqrt(1.0/(-3.0*nuV*tauD + 6.0*nuV*tauR + 6.0*tauD*tauR*(gam*nuP + nuU) - A))/(nuV*tauD)
			
	if abs(np.imag(P0Minus)) < rotSolTol and abs(np.imag(u0Minus)) < rotSolTol and abs(np.imag(phi0Minus)) < rotSolTol and abs(np.imag(omega0Minus)) < rotSolTol and abs(np.real(P0Minus)) > rotSolTol and abs(np.real(u0Minus)) > rotSolTol and abs(np.real(phi0Minus)) > rotSolTol and abs(np.real(omega0Minus)) > rotSolTol:
		nuPLow = np.real(nuP)		
	else:	
		nuPHigh = np.real(nuP)
		
nuPOscHyst = np.real(0.5*(nuPHigh + nuPLow))

print('nuP below which hysteresis (uniform states): ' + str(nuPOscStart))


# determine nuP where finite wavelength instability sets in
nuPLow = nuPPol + tol
nuPHigh = nuPOscHyst
while (nuPHigh - nuPLow) > tol:
	
	nuP = 0.5*(nuPHigh + nuPLow)

	P0 = np.sqrt(1.0 - 3.0*(nuV + nuU*tauD)/(2.0*gam*tauR*nuP))
	v0 = nuP*P0/(nuV + nuU*tauD)
	u0 = v0*tauD
	
	eigenvalsReZeroMax = 0.0
	eigenvalsReMax = 0.0
	eigenvalsImMaxAtMaxGR = 0.0
	eigenvalsImMaxAtAnyPosGR = 0.0
	
	for k in kList:
	
		kx = k
		ky = 0.0

		kSq = kx*kx + ky*ky
		if kSq == 0.0:
			kSq = 0.000000001

		MPP = 1j*np.zeros((2,2))
		MPP[0,0] = - v0*1j*kx - Dtr*kSq - 1/tauR - 4*gam*P0*v0/3
		MPP[1,1] = - v0*1j*kx - Dtr*kSq - 1/tauR - gam*P0*v0

		MPv = 1j*np.zeros((2,2))
		MPv[0,0] = kappa*1j*kx*P0 + 2*gam/3 - 2*gam*P0*P0/3
		MPv[1,0] = - 1j*ky*P0/2 + kappa*1j*ky*P0/2
		MPv[1,1] = 1j*kx*P0/2 + kappa*1j*kx*P0/2 + 2*gam/3 + gam*P0*P0/3

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
		
		if k == kZero:
			for i in range(0,4):
				if eigenvalsReZeroMax < np.real(eigenvalues[i]):
					eigenvalsReZeroMax = np.real(eigenvalues[i])
		for i in range(0,4):
			if eigenvalsReMax < np.real(eigenvalues[i]):
				eigenvalsReMax = np.real(eigenvalues[i])
				eigenvalsImMaxAtMaxGR = np.imag(eigenvalues[i])
		
		if np.max(np.real(eigenvalues)) > grRateTol:
			for i in range(0,4):
				if abs(eigenvalsImMaxAtAnyPosGR) < abs(np.imag(eigenvalues[i])):
					eigenvalsImMaxAtAnyPosGR = np.imag(eigenvalues[i])
	
	if eigenvalsReMax > grRateTol and eigenvalsReMax > eigenvalsReZeroMax:
		nuPHigh = nuP
	else:
		nuPLow = nuP

nuPFiniteWavelength = 0.5*(nuPHigh + nuPLow)


# determine nuP where long wavelengths start to grow faster than finite wavelengths
nuPLow = nuPFiniteWavelength + tol
nuPHigh = 60.0

while (nuPHigh - nuPLow) > tol:
	
	nuP = 0.5*(nuPHigh + nuPLow)

	P0 = np.sqrt(1.0 - 3.0*(nuV + nuU*tauD)/(2.0*gam*tauR*nuP))
	v0 = nuP*P0/(nuV + nuU*tauD)
	u0 = v0*tauD
	
	eigenvalsReZeroMax = 0.0
	eigenvalsReMax = 0.0
	eigenvalsImMaxAtMaxGR = 0.0
	eigenvalsImMaxAtAnyPosGR = 0.0
	
	for k in kList:
	
		kx = k
		ky = 0.0

		kSq = kx*kx + ky*ky
		if kSq == 0.0:
			kSq = 0.000000001

		MPP = 1j*np.zeros((2,2))
		MPP[0,0] = - v0*1j*kx - Dtr*kSq - 1/tauR - 4*gam*P0*v0/3
		MPP[1,1] = - v0*1j*kx - Dtr*kSq - 1/tauR - gam*P0*v0

		MPv = 1j*np.zeros((2,2))
		MPv[0,0] = kappa*1j*kx*P0 + 2*gam/3 - 2*gam*P0*P0/3
		MPv[1,0] = - 1j*ky*P0/2 + kappa*1j*ky*P0/2
		MPv[1,1] = 1j*kx*P0/2 + kappa*1j*kx*P0/2 + 2*gam/3 + gam*P0*P0/3

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
		
		if k == kZero:
			for i in range(0,4):
				if eigenvalsReZeroMax < np.real(eigenvalues[i]):
					eigenvalsReZeroMax = np.real(eigenvalues[i])
		for i in range(0,4):
			if eigenvalsReMax < np.real(eigenvalues[i]):
				eigenvalsReMax = np.real(eigenvalues[i])
				eigenvalsImMaxAtMaxGR = np.imag(eigenvalues[i])
		
		if np.max(np.real(eigenvalues)) > grRateTol:
			for i in range(0,4):
				if abs(eigenvalsImMaxAtAnyPosGR) < abs(np.imag(eigenvalues[i])):
					eigenvalsImMaxAtAnyPosGR = np.imag(eigenvalues[i])
	
	if eigenvalsReZeroMax > grRateTol and abs(eigenvalsReMax - eigenvalsReZeroMax) < grRateTol:
		nuPHigh = nuP
	else:
		nuPLow = nuP

nuPLongWavelength = 0.5*(nuPHigh + nuPLow)


# calculates velocity field
def calcVelo(Px,Py,ux,uy,nuP):
	
	PxF = np.fft.fft2(Px)
	PyF = np.fft.fft2(Py)
	uxF = np.fft.fft2(ux)
	uyF = np.fft.fft2(uy)

	# dealiasing
	PxF = PxF*dealiasingMask
	PyF = PyF*dealiasingMask
	uxF = uxF*dealiasingMask
	uyF = uyF*dealiasingMask
	
	# gradients od polar order
	Px_xF = kx*PxF
	Px_yF = ky*PxF
	Py_xF = kx*PyF
	Py_yF = ky*PyF
	Px_x = np.real(np.fft.ifft2(Px_xF))
	Px_y = np.real(np.fft.ifft2(Px_yF))
	Py_x = np.real(np.fft.ifft2(Py_xF))
	Py_y = np.real(np.fft.ifft2(Py_yF))
	
	PDivF = Px_xF + Py_yF
	PDiv_xF = kx*PDivF
	PDiv_yF = ky*PDivF
	
	# gradients of displacement
	ux_xF = kx*uxF
	ux_yF = ky*uxF
	uy_xF = kx*uyF
	uy_yF = ky*uyF
	ux_x = np.real(np.fft.ifft2(ux_xF))
	ux_y = np.real(np.fft.ifft2(ux_yF))
	uy_x = np.real(np.fft.ifft2(uy_xF))
	uy_y = np.real(np.fft.ifft2(uy_yF))
	
	uDivF = ux_xF + uy_yF
	uDiv_xF = kx*uDivF
	uDiv_yF = ky*uDivF
		
	# velocity field
	vxF = (elastMod*kLaplacian*uxF - nuU*uxF + nuP*PxF)/LinvF
	vyF = (elastMod*kLaplacian*uyF - nuU*uyF + nuP*PyF)/LinvF
			
	# pressure correction
	vDivF = kx*vxF + ky*vyF
	pvF = vDivF/kLaplacianPoisson
	vxF = vxF - kx*pvF
	vyF = vyF - ky*pvF
	
	# back to real space
	vx = np.real(np.fft.ifft2(vxF))
	vy = np.real(np.fft.ifft2(vyF))

	return vx, vy, vxF, vyF

# right-hand side, needed for calculation of Runge-Kutta steps
def rhs(dt,Px,Py,ux,uy,nuP):
	
	PxF = np.fft.fft2(Px)
	PyF = np.fft.fft2(Py)
	uxF = np.fft.fft2(ux)
	uyF = np.fft.fft2(uy)

	# dealiasing
	PxF = PxF*dealiasingMask
	PyF = PyF*dealiasingMask
	uxF = uxF*dealiasingMask
	uyF = uyF*dealiasingMask
	
	# gradients od polar order
	Px_xF = kx*PxF
	Px_yF = ky*PxF
	Py_xF = kx*PyF
	Py_yF = ky*PyF
	Px_x = np.real(np.fft.ifft2(Px_xF))
	Px_y = np.real(np.fft.ifft2(Px_yF))
	Py_x = np.real(np.fft.ifft2(Py_xF))
	Py_y = np.real(np.fft.ifft2(Py_yF))
	
	PDivF = Px_xF + Py_yF
	PDiv_xF = kx*PDivF
	PDiv_yF = ky*PDivF
	
	# gradients of displacement
	ux_xF = kx*uxF
	ux_yF = ky*uxF
	uy_xF = kx*uyF
	uy_yF = ky*uyF
	ux_x = np.real(np.fft.ifft2(ux_xF))
	ux_y = np.real(np.fft.ifft2(ux_yF))
	uy_x = np.real(np.fft.ifft2(uy_xF))
	uy_y = np.real(np.fft.ifft2(uy_yF))
	
	uDivF = ux_xF + uy_yF
	uDiv_xF = kx*uDivF
	uDiv_yF = ky*uDivF
	
	# velocity field
	vxF = (elastMod*kLaplacian*uxF - nuU*uxF + nuP*PxF)/LinvF
	vyF = (elastMod*kLaplacian*uyF - nuU*uyF + nuP*PyF)/LinvF

	# pressure correction
	vDivF = kx*vxF + ky*vyF
	pvF = vDivF/kLaplacianPoisson
	vxF = vxF - kx*pvF
	vyF = vyF - ky*pvF
	
	# back to real space
	vx = np.real(np.fft.ifft2(vxF))
	vy = np.real(np.fft.ifft2(vyF))
	
	# velocity gradients
	vx_xF = kx*vxF
	vx_yF = ky*vxF
	vy_xF = kx*vyF
	vy_yF = ky*vyF
	vx_x = np.real(np.fft.ifft2(vx_xF))
	vx_y = np.real(np.fft.ifft2(vx_yF))
	vy_x = np.real(np.fft.ifft2(vy_xF))
	vy_y = np.real(np.fft.ifft2(vy_yF))
	
	# convective derivative polar order
	convecPx = - vx*Px_x - vy*Px_y + 0.5*(vx_y - vy_x)*Py + kappa*vx_x*Px + 0.5*kappa*(vy_x + vx_y)*Py
	convecPy = - vx*Py_x - vy*Py_y + 0.5*(vy_x - vx_y)*Px + kappa*vy_y*Py + 0.5*kappa*(vx_y + vy_x)*Px
	convecPxF = np.fft.fft2(convecPx)
	convecPyF = np.fft.fft2(convecPy)
	
	# nonlinear velocity alignment terms
	velAlignPx = gam*(2.0*vx - 2.0*vx*Px*Px + vx*Py*Py - 3.0*Px*Py*vy)/3.0
	velAlignPy = gam*(2.0*vy - 2.0*vy*Py*Py + vy*Px*Px - 3.0*Py*Px*vx)/3.0
	velAlignPxF = np.fft.fft2(velAlignPx)
	velAlignPyF = np.fft.fft2(velAlignPy)
	
	# polar order time propagation via operator splitting (linear and nonlinear)
	PxNewF = np.exp(dt*LinPF)*(PxF + dt*(convecPxF + velAlignPxF))
	PyNewF = np.exp(dt*LinPF)*(PyF + dt*(convecPyF + velAlignPyF))
		
	# convection term displacement
	convecux = - vx*ux_x - vy*ux_y
	convecuy = - vx*uy_x - vy*uy_y
	convecuxF = np.fft.fft2(convecux)
	convecuyF = np.fft.fft2(convecuy)
	
	# correction term
	rhsAddTerm = vx_x*ux_x + vy_x*ux_y + vx_y*uy_x + vy_y*uy_y
	rhsAddTermF = np.fft.fft2(rhsAddTerm)
	qF = - rhsAddTermF/kLaplacianPoisson
			
	# displacement time propagation via operator splitting (linear and nonlinear)
	uxNewF = np.exp(dt*LinuF)*(uxF + dt*(convecuxF + vxF - kx*qF))
	uyNewF = np.exp(dt*LinuF)*(uyF + dt*(convecuyF + vyF - ky*qF))
	
	# pressure correction
	uDivF = kx*uxNewF + ky*uyNewF
	qvF = uDivF/kLaplacianPoisson
	uxNewF = uxNewF - kx*qvF
	uyNewF = uyNewF - ky*qvF
	
	# back to real space
	PxNew = np.real(np.fft.ifft2(PxNewF))
	PyNew = np.real(np.fft.ifft2(PyNewF))
	uxNew = np.real(np.fft.ifft2(uxNewF))
	uyNew = np.real(np.fft.ifft2(uyNewF))

	RHSPx = PxNew - Px
	RHSPy = PyNew - Py
	RHSux = uxNew - ux
	RHSuy = uyNew - uy

	return RHSPx, RHSPy, RHSux, RHSuy



# start the actual search for nuP where stripe patterns become unstable against local rotations

nuPLow = nuPFiniteWavelength + tol
nuPHigh = nuPLongWavelength - tol

print()
print('nuPLow = ' + str(nuPLow) )
print('nuPHigh = ' + str(nuPHigh) )
print()

while (nuPHigh - nuPLow) > tol:

	nuP = 0.5*(nuPHigh + nuPLow)
	
	# first determine the stationary polar solution
	P0 = np.sqrt(1.0 - 3.0*(nuV + nuU*tauD)/(2.0*gam*tauR*nuP))
	v0 = nuP*P0/(nuV + nuU*tauD)
	u0 = v0*tauD

	eigenvalsReZeroMax = 0.0
	eigenvalsReMax = 0.0
	eigenvalsImMaxAtMaxGR = 0.0
	eigenvalsImMaxAtAnyPosGR = 0.0


	# determine the wavenumber of the fastest-growing mode
	kZeroLSA = 0.0		
	kMaxLSA = 5.0
	numkLSA = 100
	kListLSA = np.linspace(kZeroLSA,kMaxLSA,num=numkLSA,endpoint=False)

	kLSAMostUnstable = 0.0

	for kLSA in kListLSA:

		kxLSA = kLSA
		kyLSA = 0.0

		kSqLSA = kxLSA*kxLSA + kyLSA*kyLSA
		if kSqLSA == 0.0:
			kSqLSA = 0.000000001

		MPP = 1j*np.zeros((2,2))
		MPP[0,0] = - v0*1j*kxLSA - Dtr*kSqLSA - 1/tauR - 4*gam*P0*v0/3
		MPP[1,1] = - v0*1j*kxLSA - Dtr*kSqLSA - 1/tauR - gam*P0*v0

		MPv = 1j*np.zeros((2,2))
		MPv[0,0] = kappa*1j*kxLSA*P0 + 2*gam/3 - 2*gam*P0*P0/3
		MPv[1,0] = - 1j*kyLSA*P0/2 + kappa*1j*kyLSA*P0/2
		MPv[1,1] = 1j*kxLSA*P0/2 + kappa*1j*kxLSA*P0/2 + 2*gam/3 + gam*P0*P0/3

		MvP = 1j*np.zeros((2,2))
		MvP[0,0] = nuP
		MvP[1,1] = nuP

		Muv = 1j*np.zeros((2,2))
		Muv[0,0] = 1
		Muv[1,1] = 1 

		proj = 1j*np.zeros((2,2))
		proj[0,0] = 1 - kxLSA*kxLSA/kSqLSA
		proj[0,1] = - kxLSA*kyLSA/kSqLSA
		proj[1,0] = - kxLSA*kyLSA/kSqLSA
		proj[1,1] = 1 - kyLSA*kyLSA/kSqLSA

		identMat = 1j*np.zeros((2,2))
		identMat[0,0] = 1
		identMat[1,1] = 1

		M = 1j*np.zeros((4,4))

		MPartPP = MPP + MPv@proj@MvP/(visc*kSqLSA + nuV)

		MPartPu = - (elastMod*kSqLSA + nuU)*MPv/(visc*kSqLSA + nuV)

		MPartuP = nuP*Muv@proj/(visc*kSqLSA + nuV)

		MPartuu = - (v0*1j*kxLSA + 1/tauD)*identMat - (elastMod*kSqLSA + nuU)*proj@Muv/(visc*kSqLSA + nuV)

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

		if kLSA == kZeroLSA:
			for i in range(0,4):
				if eigenvalsReZeroMax < np.real(eigenvalues[i]):
					eigenvalsReZeroMax = np.real(eigenvalues[i])
		for i in range(0,4):
			if eigenvalsReMax < np.real(eigenvalues[i]):
				eigenvalsReMax = np.real(eigenvalues[i])
				kLSAMostUnstable = kLSA
	
	# set the box size to four times the wavelength of the fastest-growing mode
	if kLSAMostUnstable > 0.01:
		L = 8.0*np.pi/kLSAMostUnstable
	else: 
		L = 100.0
		
	dx = L/N

	stripesAmpP = 0.02
	stripesAmpu = 0.01

	
	dx = L/N

	# to check if parameters are forwarded correctly
	print('kappa = ' + str(kappa), flush=True)
	print('gam_a = ' + str(gam), flush=True)
	print('visc = ' + str(visc), flush=True)
	print('elastMod = ' + str(elastMod), flush=True)
	print('nu_v = ' + str(nuV), flush=True)
	print('nu_u = ' + str(nuU), flush=True)
	print('nu_P = ' + str(nuP), flush=True)
	print('tauD = ' + str(tauD), flush=True)
	print('grid points = ' + str(N) + 'x' + str(N), flush=True)
	print('box size = ' + str(L), flush=True)

	# definition of the wavevector
	k = 1j*np.fft.fftfreq(N,d=L/2.0/np.pi/N)
	kx, ky = np.meshgrid(k, k, sparse=False, indexing='ij')

	# for plotting: numbering meshgrid
	nx, ny = np.meshgrid(np.linspace(0,N,num=N,endpoint=False),np.linspace(0,N,num=N,endpoint=False),indexing='ij')

	# mask for dealiaising
	dealiasingFrac = 2.0/6.0
	dealiasingMask = np.ones((N,N))
	dealiasingMask[:,round(N/2)-round(N/2*dealiasingFrac):round(N/2)+round(N/2*dealiasingFrac)] = 0
	dealiasingMask[round(N/2)-round(N/2*dealiasingFrac):round(N/2)+round(N/2*dealiasingFrac),:] = 0

	# Laplacian and biharmonic derivative in Fourier space
	kLaplacian = kx*kx + ky*ky
	kLaplacianPoisson = kLaplacian
	# kLaplacianPoisson[0,0] = -1.0
	kLaplacianPoisson = 1j*np.zeros((N,N))
	kLaplacianPoisson[:,:] = kLaplacian[:,:]
	kLaplacianPoisson[0,0] = 0.000000001 # to avoid divide-by-zero
	kBiharmonic = kLaplacian*kLaplacian

	# random initial values
	initVarP = 0.01
	initVaru = 0.0
	Px = -initVarP + 2.0*initVarP*np.random.rand(N,N)
	Py = -initVarP + 2.0*initVarP*np.random.rand(N,N)
	ux = -initVaru + 2.0*initVaru*np.random.rand(N,N)
	uy = -initVaru + 2.0*initVaru*np.random.rand(N,N)

	# scaling such that mean values are zero
	Px = Px - Px.mean() + 0.0
	Py = Py - Py.mean() + 0.0
	ux = ux - ux.mean() + 0.0
	uy = uy - uy.mean() + 0.0

	# add stripe-like modulation of small amplitude
	for i in range(0,N):
		for j in range(0,N):
			Px[i,j] = Px[i,j] + P0
			Py[i,j] = Py[i,j] + stripesAmpP*np.sin(2.0*np.pi*4.0*float(i)/N)
			ux[i,j] = ux[i,j] + u0
			uy[i,j] = uy[i,j] + stripesAmpu*np.sin(2.0*np.pi*4.0*float(i)/N)


	# to estimate remaining time
	timeSinceLastPrint = time.time()

	# linear operators
	LinPF = - 1.0/tauR + Dtr*kLaplacian
	LinvF = nuV - visc*kLaplacianPoisson
	LinuF = - 1.0/tauD
	
	# local angle of polar order parameter
	angleP = np.arctan2(Py,Px)
	anglePOld = angleP	
	angularStepPSum = np.zeros((N,N))
	angularStepPSumMean = 0.0
	
	nStep = 0
		
	while nStep < numSteps and angularStepPSumMean < 2.0*np.pi :
		
		t = round(nStep*dt,6)
		
		# calculate the local angle and determine how it has changed
		angleP = np.arctan2(Py,Px)
		angularStepP = angleP - anglePOld
		for i in range(0,N):
			for j in range(0,N):
				if angularStepP[i,j] < -0.5*np.pi:
					angularStepP[i,j] += np.pi
				if angularStepP[i,j] > 0.5*np.pi:
					angularStepP[i,j] -= np.pi
		anglePOld = angleP
		angularStepPSum = angularStepPSum + angularStepP
		
		# mean local angle change is used to identify rotational states
		angularStepPSumMean = np.mean(abs(angularStepPSum))
		
		# show progress
		if nStep%dnPrint == 0:
			vx, vy, vxF, vyF = calcVelo(Px,Py,ux,uy,nuP)
			vDivF = kx*vxF + ky*vyF		
			vDiv = np.real(np.fft.ifft2(vDivF))	
			
			uxF = np.fft.fft2(ux)
			uyF = np.fft.fft2(uy)
			uDivF = kx*uxF + ky*uyF		
			uDiv = np.real(np.fft.ifft2(uDivF))	
			
			print('t = ' + str(t), flush=True)
			print('time remaining: ' + str((time.time() - timeSinceLastPrint)/60/60*(tEnd-t)/dtPrint) + ' hours', flush=True)
			print('Px = ' + str(np.mean(Px)), flush=True)
			print('Py = ' + str(np.mean(Py)), flush=True)
			print('PRMS = ' + str(np.sqrt(np.mean(Px*Px + Py*Py))), flush=True)
			print('vx = ' + str(np.mean(vx)), flush=True)
			print('vy = ' + str(np.mean(vy)), flush=True)
			print('vRMS = ' + str(np.sqrt(np.mean(vx*vx + vy*vy))), flush=True)
			print('ux = ' + str(np.mean(ux)), flush=True)
			print('uy = ' + str(np.mean(uy)), flush=True)
			print('uRMS = ' + str(np.sqrt(np.mean(ux*ux + uy*uy))), flush=True)
			print('local rotations in P : ' + str(np.mean(abs(angularStepPSum))/2.0/np.pi))
			timeSinceLastPrint = time.time()
			

		# Runge-Kutta intermediate steps
		k1Px, k1Py, k1ux, k1uy = rhs(dt,Px,Py,ux,uy,nuP)
		k2Px, k2Py, k2ux, k2uy = rhs(dt,Px + 0.5*k1Px,Py + 0.5*k1Py,ux + 0.5*k1ux,uy + 0.5*k1uy,nuP)
		k3Px, k3Py, k3ux, k3uy = rhs(dt,Px + 0.5*k2Px,Py + 0.5*k2Py,ux + 0.5*k2ux,uy + 0.5*k2uy,nuP)
		k4Px, k4Py, k4ux, k4uy = rhs(dt,Px + k3Px,Py + k3Py,ux + k3ux,uy + k3uy,nuP)

		# calculate new velocity field at the next step
		Px = Px + (k1Px + 2.0*k2Px + 2.0*k3Px + k4Px)/6.0
		Py = Py + (k1Py + 2.0*k2Py + 2.0*k3Py + k4Py)/6.0
		ux = ux + (k1ux + 2.0*k2ux + 2.0*k3ux + k4ux)/6.0
		uy = uy + (k1uy + 2.0*k2uy + 2.0*k3uy + k4uy)/6.0
		
		nStep += 1
		
	
	print('local rotations in P : ' + str(np.mean(abs(angularStepPSum))/2.0/np.pi))

	# check if the polar order parameter has on average rotated for at least one full rotation
	if angularStepPSumMean > 2.0*np.pi - tol:
		nuPHigh = nuP
	else:
		nuPLow = nuP
		
	print()
	print('nuPLow = ' + str(nuPLow) )
	print('nuPHigh = ' + str(nuPHigh) )
	print()
		
	dataFile.write('%9.5f'%(nuPLow) + ' %9.5f'%(nuPHigh) + ' %9.5f'%(0.5*(nuPHigh + nuPLow)) + '\n')
	dataFile.flush()
	
dataFile.close()
