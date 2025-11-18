

import numpy as np


def pfss(br0,nr,rss):
	r"""
	Adapted from pfsspy, https://github.com/dstansby/pfsspy/
	pfsspy was originally by David Stansby, with the main solver from Anthony
	Yeates; this version largely removes the counterproductive OOP wrappers
	put in by David.
	pfsspy was released under GPLv3, this code under Apache-2.0

	Computes PFSS model coronal magnetic field.

	Extrapolates a 3D PFSS using an eigenfunction method in :math:`r,s,p`
	coordinates, on the dumfric grid
	(equally spaced in :math:`\rho = \ln(r/r_{sun})`,
	:math:`s= \cos(\theta)`, and :math:`p=\phi`).

	Inputs
	-----
	br0 - 2D numpy array of the radial magnetic field of shape (s,p)
	nr - integer for size of new rho grid
	rss - source surface radius (units of solar radius, 
		i.e. min(rho) = 0 => r_{inner} = r_{sun})

	Notes
	-----
	In order to avoid numerical issues, the monopole term (which should be zero
	for a physical magnetic field anyway) is explicitly excluded from the
	solution.

	The output should have zero current to machine precision,
	when computed with the DuMFriC staggered discretization.
	"""


	HAS_NUMBA = False
	try:
		import numba
		HAS_NUMBA = True
	except Exception:
		pass


	def _eigh(A):
		return np.linalg.eigh(A)

	def _compute_r_term(m, k, ns, Q, brt, lam, ffm, nr, ffp, psi, psir):
		for l in range(ns):
			# Ignore the l=0 and m=0 term; for a globally divergence free field
			# this term is zero anyway, but numerically it may be small which
			# causes numerical issues when solving for c, d
			if l == 0 and m == 0:
				continue
			# - sum (c_{lm} + d_{lm}) * lam_{l}
			# lam[l] is small so this blows up
			cdlm = np.dot(Q[:, l], np.asfortranarray(brt[:, m])) / lam[l]
			# - ratio c_{lm}/d_{lm} [numerically safer this way up]
			ratio = (ffm[l]**(nr - 1) - ffm[l]**nr) / (ffp[l]**nr - ffp[l]**(nr - 1))
			dlm = cdlm / (1.0 + ratio)
			clm = ratio * dlm
			psir[:, l] = clm * ffp[l]**k + dlm * ffm[l]**k

		# - compute entry for this m in psit = Sum_l c_{lm}Q_{lm}**j
		psi[:, :, m] = np.dot(psir, Q.T)
		return psi, psir

	def _als_alp(nr, nphi, Fs, psi, Fp, als, alp):
		for j in range(nr + 1):
			for i in range(nphi + 1):
				als[i, :, j] = Fs * (psi[j, :, ((i - 1) % nphi)] - psi[j, :, ((i) % nphi)])
			for i in range(nphi):
				alp[i, 1:-1, j] = Fp[1:-1] * (psi[j, 1:, i] - psi[j, :-1, i])
		return als, alp

	def _A_diag(A, ns, Vg, Uc, mu, m):
		for j in range(ns):
			A[j, j] = Vg[j] + Vg[j + 1] + Uc[j] * mu[m]
		return A

	if HAS_NUMBA:
		_eigh = numba.jit(nopython=True)(_eigh)
		_compute_r_term = numba.jit(nopython=True)(_compute_r_term)
		_als_alp = numba.jit(nopython=True)(_als_alp)
		_A_diag = numba.jit(nopython=True)(_A_diag)


	ns = br0.shape[0]
	nphi = br0.shape[1]

	# Coordinates:
	ds = 2.0 / ns
	dp = 2 * np.pi / nphi
	dr = np.log(rss) / nr

	rg = np.linspace(0, np.log(rss), nr + 1)
	rc = np.linspace(0.5 * dr, np.log(rss) - 0.5 * dr, nr)
	sg = np.linspace(-1, 1, ns + 1)
	sc = np.linspace(-1 + 0.5 * ds, 1 - 0.5 * ds, ns)
	pg = np.linspace(0, 2 * np.pi, nphi + 1)
	pc = np.linspace(0.5 * dp, 2 * np.pi - 0.5 * dp, nphi)

	k = np.linspace(0, nr, nr + 1)

	Fp = sg * 0  # Lp/Ls on p-ribs
	Fp[1:-1] = np.sqrt(1 - sg[1:-1]**2) / (np.arcsin(sc[1:]) - np.arcsin(sc[:-1])) * dp
	Vg = Fp / ds / dp
	Fs = (np.arcsin(sg[1:]) - np.arcsin(sg[:-1])) / np.sqrt(1 - sc**2) / dp  # Ls/Lp on s-ribs
	Uc = Fs / ds / dp

	# FFT in phi of photospheric distribution at each latitude:
	brt = np.fft.rfft(br0, axis=1)
	brt = brt.astype(np.complex128)

	# Prepare tridiagonal matrix:
	# - create off-diagonal part of the matrix:
	A = np.zeros((ns, ns))
	for j in range(ns - 1):
		A[j, j + 1] = -Vg[j + 1]
		A[j + 1, j] = A[j, j + 1]
	# - term required for m-dependent part of matrix:
	mu = np.fft.fftfreq(nphi)
	mu = 4 * np.sin(np.pi * mu)**2
	# - initialise:
	psir = np.zeros((nr + 1, ns), dtype='complex')
	psi = np.zeros((nr + 1, ns, nphi), dtype='complex')
	e1 = np.exp(dr)
	fact = np.sinh(dr) * (e1 - 1)

	# Loop over azimuthal modes (positive m):
	for m in range(nphi // 2 + 1):
		# - set diagonal terms of matrix:
		A = _A_diag(A, ns, Vg, Uc, mu, m)

		# - compute eigenvectors Q_{lm} and eigenvalues lam_{lm}:
		#   (note that A is symmetric so use special solver)
		lam, Q = _eigh(A)
		Q = Q.astype(np.complex128)
		# - solve quadratic:
		Flm = 0.5 * (1 + e1 + lam * fact)
		ffp = Flm + np.sqrt(Flm**2 - e1)
		ffm = e1 / ffp

		# - compute radial term for each l (for this m):
		psi, psir = _compute_r_term(m, k, ns, Q, brt, lam, ffm, nr, ffp, psi, psir)

		if (m > 0):
			psi[:, :, nphi - m] = np.conj(psi[:, :, m])

	# Past this point only psi, Fs, Fp are needed
	# Compute psi by inverse fft:
	psi = np.real(np.fft.ifft(psi, axis=2))

	# Hence compute vector potential [note index order, for netcdf]:
	# (note that alr is zero by definition)
	alr = np.zeros((nphi + 1, ns + 1, nr))
	als = np.zeros((nphi + 1, ns, nr + 1))
	alp = np.zeros((nphi, ns + 1, nr + 1))

	als, alp = _als_alp(nr, nphi, Fs, psi, Fp, als, alp)

	rcp2 = np.linspace(-0.5 * dr, np.log(rss) + 0.5 * dr, nr + 2)
	rrc = np.exp(rcp2)
	thc = np.zeros(ns + 2) - 1
	thc[1:-1] = np.arccos(sc)
	# Centre of cells in phi (including ghost cells)
	np.linspace(-0.5 * dp, 2 * np.pi + 0.5 * dp, nphi + 2)

	# Required face normals:
	dnp = np.zeros((ns + 2, 2))
	dns = np.zeros((ns + 1, 2))
	dnr = np.zeros(ns + 2)
	for k in range(2):
		for j in range(1, ns + 1):
			dnp[j, k] = rrc[k] * np.sqrt(1 - sc[j - 1]**2) * dp
		dnp[0, k] = dnp[1, k]
		dnp[-1, k] = dnp[-2, k]
		for j in range(1, ns):
			dns[j, k] = rrc[k] * (np.arcsin(sc[j]) - np.arcsin(sc[j - 1]))
		dns[0, k] = dns[1, k]
		dns[-1, k] = dns[-2, k]
	for j in range(ns + 2):
		dnr[j] = rrc[0] * (np.exp(dr) - 1)
	dnr[0] = -dnr[0]
	dnr[-1] = -dnr[-1]

	# Required area factors:
	Sbr = np.zeros((ns + 2, nr + 1))
	for k in range(nr + 1):
		Sbr[1:-1, k] = np.exp(2 * rg[k]) * ds * dp
		Sbr[0, k] = Sbr[1, k]
		Sbr[-1, k] = Sbr[-2, k]
	Sbs = np.zeros((ns + 1, nr + 2))
	for k in range(nr + 2):
		for j in range(1, ns):
			Sbs[j, k] = 0.5 * np.exp(2 * rcp2[k] - dr) * dp * (np.exp(2 * dr) - 1) * np.sqrt(1 - sg[j]**2)
		Sbs[0, k] = Sbs[1, k]
		Sbs[-1, k] = Sbs[-2, k]
	Sbp = np.zeros((ns + 2, nr + 2))
	for k in range(nr + 2):
		for j in range(1, ns + 1):
			Sbp[j, k] = 0.5 * np.exp(2 * rcp2[k] - dr) * (np.exp(2 * dr) - 1) * (np.arcsin(sg[j]) - np.arcsin(sg[j - 1]))
		Sbp[0, k] = Sbp[1, k]
		Sbp[-1, k] = Sbp[-2, k]

	# Compute br*Sbr, bs*Sbs, bp*Sbp at cell centres by Stokes theorem:
	br = np.zeros((nphi + 2, ns + 2, nr + 1))
	bs = np.zeros((nphi + 2, ns + 1, nr + 2))
	bp = np.zeros((nphi + 1, ns + 2, nr + 2))
	br[1:-1, 1:-1, :] = als[1:, :, :] - als[:-1, :, :] + alp[:, :-1, :] - alp[:, 1:, :]
	bs[1:-1, :, 1:-1] = alp[:, :, 1:] - alp[:, :, :-1]
	bp[:, 1:-1, 1:-1] = als[:, :, :-1] - als[:, :, 1:]

	# Fill ghost values with boundary conditions:
	# - zero-gradient at outer boundary:
	bs[1:-1, :, -1] = 2 * bs[1:-1, :, -2] - bs[1:-1, :, -3]
	bp[:, 1:-1, -1] = 2 * bp[:, 1:-1, -2] - bp[:, 1:-1, -3]
	# - periodic in phi:
	bs[0, :, :] = bs[-2, :, :]
	bs[-1, :, :] = bs[1, :, :]
	br[0, :, :] = br[-2, :, :]
	br[-1, :, :] = br[1, :, :]
	# js = jp = 0 at photosphere:
	for i in range(nphi + 1):
		bp[i, :, 0] = Sbp[:, 0] / dnp[:, 0] * (bp[i, :, 1] * dnp[:, 1] / Sbp[:, 1] + br[i, :, 0] * dnr[:] / Sbr[:, 0] - br[i + 1, :, 0] * dnr[:] / Sbr[:, 0])
	for i in range(nphi + 2):
		bs[i, :, 0] = Sbs[:, 0] / dns[:, 0] * (bs[i, :, 1] * dns[:, 1] / Sbs[:, 1] + br[i, :-1, 0] * dnr[:-1] / Sbr[:-1, 0] - br[i, 1:, 0] * dnr[1:] / Sbr[1:, 0])
	# - polar boundaries as in dumfric:
	for i in range(nphi + 2):
		i1 = (i + nphi // 2) % nphi
		br[i, -1, :] = br[i1, -2, :]
		br[i, 0, :] = br[i1, 1, :]
		bs[i, -1, :] = 0.5 * (bs[i, -2, :] - bs[i1, -2, :])
		bs[i, 0, :] = 0.5 * (bs[i, 1, :] - bs[i1, 1, :])
	for i in range(nphi + 1):
		i1 = (i + nphi // 2) % nphi
		bp[i, -1, :] = -bp[i1, -2, :]
		bp[i, 0, :] = -bp[i1, 1, :]

	brc = br.copy()
	bsc = bs.copy()
	bpc = bp.copy()
	for i in range(nphi + 2):
		brc[i, :, :] = br[i, :, :] / Sbr
		bsc[i, :, :] = bs[i, :, :] / Sbs
	for i in range(nphi + 1):
		bpc[i, :, :] = bp[i, :, :] / Sbp

	brc = brc[1:-1, 1:-1, :]
	bsc = -bsc[1:-1, :, 1:-1]
	bpc = bpc[:, 1:-1, 1:-1]

	brg = br[:-1, :-1, :] + br[1:, :-1, :] + br[1:, 1:, :] + br[:-1, 1:, :]
	bsg = bs[:-1, :, :-1] + bs[1:, :, :-1] + bs[1:, :, 1:] + bs[:-1, :, 1:]
	bpg = bp[:, :-1, :-1] + bp[:, 1:, :-1] + bp[:, 1:, 1:] + bp[:, :-1, 1:]
	for i in range(nphi + 1):
		brg[i, :, :] /= 2 * (Sbr[:-1, :] + Sbr[1:, :])
		bsg[i, :, :] /= 2 * (Sbs[:, :-1] + Sbs[:, 1:])
	for i in range(nphi + 1):
		bpg[i, :, :] /= (Sbp[:-1, :-1] + Sbp[1:, :-1] + Sbp[1:, 1:] + Sbp[:-1, 1:])
	bsg *= -1


	return rg,rc,sg,sc,pg,pc,np.stack((brg, bsg, bpg), axis=-1),(brc, bsc, bpc)


def resample(br0,new_shape):
	"""
	Resize surface magnetic the same way that sunpy handles a magnetogram
	with its own Map.resample method with default (linear) methods

	Inputs
	-----
	br0 - 2D numpy array of the radial magnetic field of shape (s,p)
	new_shape - tuple containing 2 integers (new_s, new_p)

	"""
	from scipy.interpolate import RegularGridInterpolator

	num_s = new_shape[0]
	num_p = new_shape[1]

	#Weird 0.5 pixel shift from sunpy (result of center=True argument by default)
	interp = RegularGridInterpolator((np.linspace(0.5,br0.shape[0]-0.5,num=br0.shape[0]),np.linspace(0.5,br0.shape[1]-0.5,num=br0.shape[1])),br0)
	s_out = (np.arange(num_s) + 0.5) * br0.shape[0]/num_s
	p_out = (np.arange(num_p) + 0.5) * br0.shape[1]/num_p

	s_out,p_out = np.meshgrid(s_out,p_out,indexing='ij')
	Br_out = interp((s_out,p_out))
	return Br_out


