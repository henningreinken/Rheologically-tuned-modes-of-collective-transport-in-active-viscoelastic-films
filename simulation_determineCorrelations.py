
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
tauD = 100.0
nuP = 20.0

N = 128
L = 100.0
tEnd = 20.0

dt = 0.001
dx = L/N

runNum = 0

nStep = 0

numSteps = round(tEnd/dt) 

dtPrint = 1.0
dnPrint = round(dtPrint/dt)

dtMean = 0.01
dnMean = round(dtMean/dt)

tStartMean = 1.0

NMean = int(round((tEnd - tStartMean)/dtMean))
nMean = 0

PMeanData = np.zeros((NMean,4))
vMeanData = np.zeros((NMean,4))
uMeanData = np.zeros((NMean,4))
defectsMeanData = np.zeros((NMean,5))

NCorr = int(N/2)
phaseDiffPvCorr = np.zeros(NCorr)

# to check if parameters are forwarded correctly
print('gamma_a = ' + str(gam), flush=True)
print('visc = ' + str(visc), flush=True)
print('elastMod = ' + str(elastMod), flush=True)
print('kappa = ' + str(kappa), flush=True)
print('nu_v = ' + str(nuV), flush=True)
print('nu_u = ' + str(nuU), flush=True)
print('nu_P = ' + str(nuP), flush=True)
print('tauD = ' + str(tauD), flush=True)
print('grid points = ' + str(N) + 'x' + str(N), flush=True)
print('box size = ' + str(L), flush=True)
print('run ' + str(runNum), flush=True)
print()

# get current directory
thisPath = os.getcwd()
# where to save everything
folderName = thisPath

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

# initial values: data from data files
if runNum != 0:
	
	Px = np.genfromtxt(folderName+'/dataPxInit%04d.dat'%(runNum))
	Py = np.genfromtxt(folderName+'/dataPyInit%04d.dat'%(runNum))
	vx = np.genfromtxt(folderName+'/datavxInit%04d.dat'%(runNum))
	vy = np.genfromtxt(folderName+'/datavyInit%04d.dat'%(runNum))
	ux = np.genfromtxt(folderName+'/datauxInit%04d.dat'%(runNum))
	uy = np.genfromtxt(folderName+'/datauyInit%04d.dat'%(runNum))

# counts the number of times field is saved
numbering = 0

# to estimate remaining time
timeSinceLastPrint = time.time()

# linear operators
LinPF = - 1.0/tauR + Dtr*kLaplacian
LinvF = nuV - visc*kLaplacianPoisson
LinuF = - 1.0/tauD

def calcVelo(Px,Py,ux,uy):
	
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
def rhs(dt,Px,Py,ux,uy):
	
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



for nStep in range(0,numSteps):
	
	t = round(nStep*dt,6)
	
	if nStep%dnMean == 0 and t>= tStartMean:
		
		vx, vy, vxF, vyF = calcVelo(Px,Py,ux,uy)

		PxMean = np.mean(Px)
		PyMean = np.mean(Py)
		PRMS = np.sqrt(np.mean(Px*Px + Py*Py))

		vxMean = np.mean(vx)
		vyMean = np.mean(vy)
		vRMS = np.sqrt(np.mean(vx*vx + vy*vy))

		uxMean = np.mean(ux)
		uyMean = np.mean(uy)
		uRMS = np.sqrt(np.mean(ux*ux + uy*uy))

		PMeanData[nMean,0] = dt*nStep
		PMeanData[nMean,1] = PxMean
		PMeanData[nMean,2] = PyMean
		PMeanData[nMean,3] = PRMS

		vMeanData[nMean,0] = dt*nStep
		vMeanData[nMean,1] = vxMean
		vMeanData[nMean,2] = vyMean
		vMeanData[nMean,3] = vRMS

		uMeanData[nMean,0] = dt*nStep
		uMeanData[nMean,1] = uxMean
		uMeanData[nMean,2] = uyMean
		uMeanData[nMean,3] = uRMS
		
		angleP = np.arctan2(Py,Px)
		anglev = np.arctan2(vy,vx)
		angleDiffPv = angleP - anglev
		angleDiffPv[angleDiffPv>np.pi] = angleDiffPv[angleDiffPv>np.pi] - 2.0*np.pi
		angleDiffPv[angleDiffPv<-np.pi] = angleDiffPv[angleDiffPv<-np.pi] + 2.0*np.pi
		
		for k in range(0,NCorr):
			phaseDiffPvCorr[k] = phaseDiffPvCorr[k] + 0.25*np.mean(angleDiffPv*(np.roll(angleDiffPv,k,axis=0) + np.roll(angleDiffPv,-k,axis=0) + np.roll(angleDiffPv,k,axis=1) + np.roll(angleDiffPv,-k,axis=1)))
		
		angleP[angleP>np.pi] = angleP[angleP>np.pi] - 2.0*np.pi
		angleP[angleP<-np.pi] = angleP[angleP<-np.pi] + 2.0*np.pi
		
		# topological charge of grid points
		topoChargeP = np.zeros((N,N))
		topoChargev = np.zeros((N,N))
		
		# determine topological charges
		for i in range(0,N):
			for j in range(0,N):
				
				# periodic boundary conditions
				iP = i + 1
				if iP == N: iP = 0
				jP = j + 1
				if jP == N: jP = 0
				iM = i - 1
				if iM == -1: iM = N-1
				jM = j - 1
				if jM == -1: jM = N-1
				
				topoChargeP[i,j] += np.fmod(angleP[iP,jP] - angleP[iP,j ],np.pi)/np.pi
				topoChargeP[i,j] += np.fmod(angleP[i ,jP] - angleP[iP,jP],np.pi)/np.pi
				topoChargeP[i,j] += np.fmod(angleP[iM,jP] - angleP[i ,jP],np.pi)/np.pi
				topoChargeP[i,j] += np.fmod(angleP[iM,j ] - angleP[iM,jP],np.pi)/np.pi
				topoChargeP[i,j] += np.fmod(angleP[iM,jM] - angleP[iM,j ],np.pi)/np.pi
				topoChargeP[i,j] += np.fmod(angleP[i ,jM] - angleP[iM,jM],np.pi)/np.pi
				topoChargeP[i,j] += np.fmod(angleP[iP,jM] - angleP[i ,jM],np.pi)/np.pi
				topoChargeP[i,j] += np.fmod(angleP[iP,j ] - angleP[iP,jM],np.pi)/np.pi
			
				topoChargev[i,j] += np.fmod(anglev[iP,jP] - anglev[iP,j ],np.pi)/np.pi
				topoChargev[i,j] += np.fmod(anglev[i ,jP] - anglev[iP,jP],np.pi)/np.pi
				topoChargev[i,j] += np.fmod(anglev[iM,jP] - anglev[i ,jP],np.pi)/np.pi
				topoChargev[i,j] += np.fmod(anglev[iM,j ] - anglev[iM,jP],np.pi)/np.pi
				topoChargev[i,j] += np.fmod(anglev[iM,jM] - anglev[iM,j ],np.pi)/np.pi
				topoChargev[i,j] += np.fmod(anglev[i ,jM] - anglev[iM,jM],np.pi)/np.pi
				topoChargev[i,j] += np.fmod(anglev[iP,jM] - anglev[i ,jM],np.pi)/np.pi
				topoChargev[i,j] += np.fmod(anglev[iP,j ] - anglev[iP,jM],np.pi)/np.pi
				
		# to save the locations of defects
		defectPosLocP = np.array([]).reshape(0,2)
		defectNegLocP = np.array([]).reshape(0,2)
		defectPosLocv = np.array([]).reshape(0,2)
		defectNegLocv = np.array([]).reshape(0,2)
		
		# determine locations of defects
		for i in range(0,N):
			for j in range(0,N):
				
				# periodic boundary conditions
				iP = i + 1
				if iP == N: iP = 0	
				jP = j + 1
				if jP == N: jP = 0
				iM = i - 1
				if iM == -1: iM = N-1
				jM = j - 1
				if jM == -1: jM = N-1
				
				if topoChargeP[i,j] > 0.9 and topoChargeP[iP,j] > 0.9 and topoChargeP[iP,jP] > 0.9 and topoChargeP[i,jP] > 0.9:
					defectPosLocP = np.vstack((defectPosLocP,np.array([i*dx+0.5*dx,j*dx+0.5*dx])))
					
				if topoChargeP[i,j] < -0.9 and topoChargeP[iP,j] < -0.9 and topoChargeP[iP,jP] < -0.9 and topoChargeP[i,jP] < -0.9:
					defectNegLocP = np.vstack((defectNegLocP,np.array([i*dx+0.5*dx,j*dx+0.5*dx])))
			
				if topoChargev[i,j] > 0.9 and topoChargev[iP,j] > 0.9 and topoChargev[iP,jP] > 0.9 and topoChargev[i,jP] > 0.9:
					defectPosLocv = np.vstack((defectPosLocv,np.array([i*dx+0.5*dx,j*dx+0.5*dx])))
					
				if topoChargev[i,j] < -0.9 and topoChargev[iP,j] < -0.9 and topoChargev[iP,jP] < -0.9 and topoChargev[i,jP] < -0.9:
					defectNegLocv = np.vstack((defectNegLocv,np.array([i*dx+0.5*dx,j*dx+0.5*dx])))
		
		# determine number of defects
		numPlus1DefectsP  = np.shape(defectPosLocP)[0]
		numMinus1DefectsP = np.shape(defectNegLocP)[0]
		numPlus1Defectsv  = np.shape(defectPosLocv)[0]
		numMinus1Defectsv = np.shape(defectNegLocv)[0]
		
		# save number of defects
		defectsMeanData[nMean,0] = dt*nStep
		defectsMeanData[nMean,1] = numPlus1DefectsP	
		defectsMeanData[nMean,2] = numMinus1DefectsP	
		defectsMeanData[nMean,3] = numPlus1Defectsv	
		defectsMeanData[nMean,4] = numMinus1Defectsv
		
		nMean += 1
	
	# show progress
	if nStep%dnPrint == 0:
		vx, vy, vxF, vyF = calcVelo(Px,Py,ux,uy)
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
		timeSinceLastPrint = time.time()

	# Runge-Kutta intermediate steps
	k1Px, k1Py, k1ux, k1uy = rhs(dt,Px,Py,ux,uy)
	k2Px, k2Py, k2ux, k2uy = rhs(dt,Px + 0.5*k1Px,Py + 0.5*k1Py,ux + 0.5*k1ux,uy + 0.5*k1uy)
	k3Px, k3Py, k3ux, k3uy = rhs(dt,Px + 0.5*k2Px,Py + 0.5*k2Py,ux + 0.5*k2ux,uy + 0.5*k2uy)
	k4Px, k4Py, k4ux, k4uy = rhs(dt,Px + k3Px,Py + k3Py,ux + k3ux,uy + k3uy)

	# calculate new fields at the next step
	Px = Px + (k1Px + 2.0*k2Px + 2.0*k3Px + k4Px)/6.0
	Py = Py + (k1Py + 2.0*k2Py + 2.0*k3Py + k4Py)/6.0
	ux = ux + (k1ux + 2.0*k2ux + 2.0*k3ux + k4ux)/6.0
	uy = uy + (k1uy + 2.0*k2uy + 2.0*k3uy + k4uy)/6.0
	
# save mean values in data files
np.savetxt(folderName+'/dataPMean.dat',np.asarray(PMeanData))
np.savetxt(folderName+'/datavMean.dat',np.asarray(vMeanData))
np.savetxt(folderName+'/datauMean.dat',np.asarray(uMeanData))
np.savetxt(folderName+'/dataDefectsMean%04d.dat'%(runNum),np.asarray(defectsMeanData))

# save spatial correlation of phase shift between P and v
phaseDiffPvCorr/float(nMean)
np.savetxt(folderName+'/dataPhaseDiffCorrPv%04d.dat'%(runNum),phaseDiffPvCorr)

vx, vy, vxF, vyF = calcVelo(Px,Py,ux,uy)
	
# save last time step for the field for continuation of calculation
np.savetxt(folderName+'/dataPxInit%04d.dat'%(runNum+1),Px)
np.savetxt(folderName+'/dataPyInit%04d.dat'%(runNum+1),Py)
np.savetxt(folderName+'/datavxInit%04d.dat'%(runNum+1),vx)
np.savetxt(folderName+'/datavyInit%04d.dat'%(runNum+1),vy)
np.savetxt(folderName+'/datauxInit%04d.dat'%(runNum+1),ux)
np.savetxt(folderName+'/datauyInit%04d.dat'%(runNum+1),uy)

