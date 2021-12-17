import numpy as np
import pandas as pd
import astropy.constants as cst
from sherpa.models.model import ArithmeticModel
from sherpa.models.parameter import Parameter
import edibles.src.math as eMath
from scipy.interpolate import interp1d,interp2d

class BoltzDistribution():
    def __init__(self, n=10):
        if n < 1: n=10
        if n > 50: n=50
        self.B = 0.4305883
        self.D = 1.437 / 10 ** 6
        self.H = 1.129 / 10 ** 10
        self.Jlvl = np.arange(n)*2
        JJp1 = self.Jlvl * (self.Jlvl + 1)

        # assuming E(J) = BJ(J+1) - D[J(J+1)]^2 + H[J(J+1)]^3
        self.JdE = self.B * JJp1 - self.D * JJp1 + self.H * JJp1
        self.JdE = self.JdE * 1.2398/10000  # in the unit of EV
        self.k = 8.6173/100000 # in the unit of EV
        self.Jpop = np.ones_like(self.Jlvl)
        self.Jpop[0] = 1

    def calc(self,t_low=30, t_high=100):
        if t_low <= 0: t_low=0.0001
        kT_low = self.k * t_low
        Jpop_low = (2 * self.Jlvl + 1) * np.exp(-(self.JdE / kT_low))
        Jpop_low = Jpop_low / Jpop_low[7]

        if t_high <= 0: t_high = 0.0001
        kT_high = self.k * t_high
        Jpop_high = (2 * self.Jlvl + 1) * np.exp(-(self.JdE / kT_high))
        Jpop_high = Jpop_high / Jpop_high[7]

        Jpop = np.append(Jpop_low[0:8], Jpop_high[8:])
        self.Jpop = Jpop / np.sum(Jpop)

class BoltzModel(ArithmeticModel):
    def __init__(self, T_low=30, T_high=150, N_tot=3, b=1.6, v_off=0, name=None):
        if name is None: name="BoltzModel"
        #b_inst = cst.c.to('km/s').value / 80000 / np.sqrt(2 * np.log(2)) / 2
        linelist = "/Users/haoyufan/C2C3Data/C3LineList.txt"
        self.linelist = pd.read_csv(linelist)

        self.T_low = Parameter(name,'Temperature Low J', T_low, frozen=False, min=0.1)
        self.T_high = Parameter(name, 'Temperature High J', T_high, frozen=False, min=0.1)
        self.N_tot = Parameter(name,'ColumnDensity E12', N_tot, frozen=False, min=0.00001)
        self.b = Parameter(name, 'b_sigma', b, frozen=False, min=1.27, max=5.0)
        self.v_off = Parameter(name, 'VelocityOffset', v_off, frozen=False, min=v_off-10, max=v_off+10)

        ArithmeticModel.__init__(self, name,
                                 (self.T_low, self.T_high, self.N_tot, self.b, self.v_off))

    def calc(self, pars, x, *args, **kwargs):
        T_low, T_high, N_tot, b, v_off = pars
        flux = np.ones_like(x)

        Boltzmann = BoltzDistribution(n=50)
        Boltzmann.calc(t_low=T_low, t_high=T_high)
        Jpop = Boltzmann.Jpop
        Jlvl = Boltzmann.Jlvl

        for i in range(len(self.linelist)):
            data_slice = self.linelist.iloc[i]
            name, wave, fjj = data_slice.Name, data_slice.Wave, data_slice.fjj

            if "_" in name: name = name.split("_")[0]
            wave = (1 + v_off / cst.c.to('km/s').value) * wave
            jlvl = int(name[1:])
            N_J = N_tot * Jpop[Jlvl == jlvl][0] * 1e12  # N_tot in unit of 10E12
            flux_line = eMath.gaussianAbsorption(lam=x, lam_0=wave, b=b, N=N_J, f=fjj)
            flux = flux * flux_line

        return flux

class C2Model(ArithmeticModel):
    def __init__(self, T=30, dens=50, N_tot=10, b=0.5, d=0.0005, v_off=0, name=None, Jmax=16):
        if name is None: name = "C2Model"

        linelist = "/Users/haoyufan/C2C3Data/C2LineList.txt"
        self.linelist = pd.read_csv(linelist)

        popGrid = "/Users/haoyufan/C2C3Data/vDB82_McCall_pops.txt"
        self.popInterp = {}
        popData = np.loadtxt(popGrid, unpack=1)
        densities = np.unique(popData[0])
        temperatures = np.unique(popData[1])
        for k in range(11):
            pops = (popData[k+2] / popData[13]).reshape([len(densities), len(temperatures)])
            self.popInterp[k*2] = interp2d(temperatures, densities, pops, kind='cubic', bounds_error=True)

        if Jmax > 16: Jmax = 16
        if Jmax < 0: Jmax = 2
        self.Jmax = Jmax

        self.T = Parameter(name,'Temperature', T, frozen=False, min=10, max=100)
        self.dense = Parameter(name,'Density', dens, frozen=False, min=10, max=1000)
        self.N_tot = Parameter(name, 'ColumnDensity E13', N_tot, frozen=False, min=0.00001)
        self.b = Parameter(name, 'b_sigma', b, frozen=False, min=0.000001, max=3.0)
        self.d = Parameter(name, 'd_gamma', d, frozen=False, min=0, max=3.0)
        self.v_off = Parameter(name, 'VelocityOffset', v_off, frozen=False, min=v_off-10, max=v_off+10)

        ArithmeticModel.__init__(self, name,
                                 (self.T, self.dense, self.N_tot, self.b, self.d, self.v_off))

    def calc(self, pars, x, *args, **kwargs):
        T, dense, N_tot, b, d, v_off = pars

        """
        # to resample, but this takes a lot of times...     
        deltaGrid = (x[-1] - x[0]) / len(x)
        psfkms = 2.8
        deltakms = np.min([1.665 * b / 3., psfkms/3.])
        deltaGood = (deltakms / cst.c.to("km/s").value) * np.median(x)
        deltaNew = deltaGrid
        i = 1
        while (deltaGood < deltaNew):
            i += 1.
            deltaNew = deltaGrid / i
        x_resample = np.arange(x[0], x[-1], deltaNew)
        """

        x_resample = x

        flux = np.ones_like(x_resample)

        for i in range(len(self.linelist)):
            data_slice = self.linelist.iloc[i]
            J, wave, fjj = data_slice.J, data_slice.Wave, data_slice.fjj
            pop = self.popInterp[J](T, dense)[0]
            N_J = pop * N_tot * 1e13
            lam_0 = wave * (1 + v_off/cst.c.to("km/s").value)
            transmission_line = eMath.voigtOpticalDepthAbsorption(lam=x_resample, lam_0=lam_0, b=b, d=d, N=N_J, f=fjj)
            flux = flux * transmission_line

        modelInterp = interp1d(x_resample, flux, kind='linear', bounds_error=False, fill_value=1)
        flux_out = modelInterp(x)
        return flux_out

    def exportPop(self):
        pop = []
        for J in [0,2,4,6,8,10,12,14,16]:
            pop.append(self.popInterp[J](self.T.val, self.dense.val)[0])
        return np.array(pop)

class IsotopicC2Model(ArithmeticModel):
    def __init__(self, T=30, N_tot=10, b=1.69, v_off=0, name=None):
        if name is None:
            name = "Isotopic_C2Model"

        self.linelist = np.array([8804.892, 8805.597, 8806.655, 8808.065])

        popGrid = "/Users/haoyufan/C2C3Data/13C2_Excitation.csv"
        popData = pd.read_csv(popGrid)
        popData = popData.to_numpy().T
        temperatures = np.unique(popData[0])
        self.popInterp = {}
        for i in range(4):
            self.popInterp[i+1] = interp1d(temperatures, popData[i+1], kind='cubic')

        self.T = Parameter(name,'Temperature', T, frozen=False, min=10, max=100)
        self.N_tot = Parameter(name, 'ColumnDensity E12', N_tot, frozen=False, min=0.00001)
        self.b = Parameter(name, 'b_sigma', b, frozen=False, min=2.8 / 2.355 * np.sqrt(2), max=3.0)
        self.v_off = Parameter(name, 'VelocityOffset', v_off, frozen=False, min=v_off-10, max=v_off+10)

        ArithmeticModel.__init__(self, name, (self.T, self.N_tot, self.b, self.v_off))

    def calc(self, pars, x, *args, **kwargs):
        T, N_tot, b, v_off = pars
        N_tot = N_tot * 10 ** 12
        flux = np.ones_like(x)

        intensities = [float(self.popInterp[i](T)) for i in [1,2,3,4]]
        for wave, portion in zip(self.linelist, intensities):
            N_J = portion * N_tot
            lam_0 = wave * (1 + v_off / cst.c.to("km/s").value)
            transmission_line = eMath.voigtOpticalDepthAbsorption(lam=x, lam_0=lam_0, b=b, d=0, N=N_J, f=0.0014)
            flux = flux * transmission_line
        return flux