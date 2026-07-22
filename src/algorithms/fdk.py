"""FDK (Feldkamp-Davis-Kress) cone-beam reconstruction algorithm."""
import math

import numpy as np
from scipy.fft import fft, ifft

from ..optional_deps import PI


class OptimizedFDKReconstructor:
    """Optimized FDK reconstruction from ConeBeam_FDK_corrected2.ipynb"""

    def __init__(self, SOD=1000.0, delta_dd=0.1075, Nimage=100):
        self.SOD = SOD
        self.delta_dd = delta_dd
        self.Nimage = Nimage

    def filter_SL(self, N, d):
        """Shepp-Logan filter"""
        fh_SL = np.zeros(N)
        for k1 in range(0, N, 1):
            fh_SL[k1] = -2.0/(PI*PI*d*d*(4*(k1-N/2.0)**2-1))
        return fh_SL

    def nearestPowerOf2(self, N):
        """Find nearest power of 2"""
        a = int(math.log2(N))
        if 2**a == N:
            return N
        return 2**(a + 1)

    def Fun_Weigth_Projection(self, projection_beta, SOD, delta_dd):
        """Weight projection"""
        Nrows, Ncolumns = projection_beta.shape
        dd_column = delta_dd*np.arange(-Ncolumns/2+0.5, (Ncolumns/2+1)-0.5, 1.0)
        dd_row = delta_dd*np.arange(-Nrows/2+0.5, (Nrows/2+1)-0.5, 1.0)
        dd_row2D, dd_column2D = np.meshgrid(dd_row, dd_column, indexing='ij')
        weighted_projection = projection_beta*SOD/np.sqrt(SOD*SOD+np.power(dd_row2D, 2.0)+np.power(dd_column2D, 2.0))
        return weighted_projection

    def optimize_convolution(self, weighted_projection, fh_RL):
        """Optimized convolution"""
        Nrows, Ncolumns = weighted_projection.shape
        Nfft = self.nearestPowerOf2(2 * Ncolumns - 1)
        fh_RL_padded = np.zeros(Nfft)
        fh_RL_padded[:len(fh_RL)] = fh_RL / 2.0

        fh_RL_fft = fft(fh_RL_padded)

        projection_padded = np.zeros((Nrows, Nfft))
        projection_padded[:, :Ncolumns] = weighted_projection

        projection_fft = fft(projection_padded, axis=1)
        convoluted_freq = projection_fft * fh_RL_fft
        convoluted_time = ifft(convoluted_freq, axis=1).real
        # filter_SL builds the kernel with its peak at sample Nfft/2, not at 0, so
        # the circular convolution comes out delayed by that much. Reading from
        # offset 0 returns the wrong half of the signal (wrong shape and sign).
        filtered_projection = convoluted_time[:, Nfft // 2:Nfft // 2 + Ncolumns]

        return filtered_projection

    def Fun_BackProjection(self, filtered_projection, SOD, beta_num, beta_m, delta_dd, Nimage):
        """Backprojection"""
        Nrows, Ncolumns = filtered_projection.shape
        MX, MZ = Nimage, int(Nimage*Nrows/Ncolumns)

        roi = delta_dd*np.array([-Ncolumns/2.0+0.5, Ncolumns/2.0-0.5, -Nrows/2.0+0.5, Nrows/2.0-0.5])
        hx = (roi[1]-roi[0])/(MX-1)
        xrange = roi[0]+hx*np.arange(0, MX)
        hy = (roi[3]-roi[2])/(MZ-1)
        yrange = roi[2]+hy*np.arange(0, MZ)
        XX, YY, ZZ = np.meshgrid(xrange, xrange, yrange, indexing='ij')
        temp_rec = np.zeros((MX, MX, MZ))
        U = (SOD+XX*np.sin(beta_m)-YY*np.cos(beta_m))/SOD
        a = (XX*np.cos(beta_m)+YY*np.sin(beta_m))/U
        xx = np.int32(np.floor(a/delta_dd))
        u1 = a/delta_dd-xx
        b = ZZ/U
        yy = np.int32(np.floor(b/delta_dd))
        u2 = b/delta_dd-yy
        xx = xx+int(Ncolumns/2)
        yy = yy+int(Nrows/2)

        mask = np.where((xx >= 0) & (xx < Ncolumns-1) & (yy >= 0) & (yy < Nrows-1))
        xx = xx[mask]
        yy = yy[mask]
        u1 = u1[mask]
        u2 = u2[mask]

        temp = ((1-u1)*(1-u2)*filtered_projection[yy, xx]+(1-u1)*u2*filtered_projection[yy+1, xx]+
                (1-u2)*u1*filtered_projection[yy, xx+1]+u1*u2*filtered_projection[yy+1, xx+1])
        temp_rec[mask] = temp_rec[mask]+temp/(np.power(U[mask], 2))*2*PI/beta_num

        return temp_rec
