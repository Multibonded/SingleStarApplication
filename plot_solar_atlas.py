import ISPy.spec.atlas as atlas
import matplotlib.pyplot as plt
from specutils import analysis
from specutils import Spectrum1D
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
from lmfit.models import GaussianModel, VoigtModel, LinearModel, ConstantModel


lines = [4281.363, 4465.805, 4758.118, 5022.866, 5490.147, 6092.789, 7357.726, 4409.52, 4657.212, 4719.533]
for line in lines:
    width = 0.1#[half spectral range in angstrom, so use like like 1 or 2]
    s = atlas.atlas()
    lambda_atlas,intensity_atlas,c = s.get(line-width,  line+width, cgs=True)

    intensity_atlas = intensity_atlas/max(intensity_atlas)

    spec1d = Spectrum1D(spectral_axis=lambda_atlas * u.AA, flux=intensity_atlas * u.Jy)
    ew = analysis.equivalent_width(spec1d)
    print(line)
    print("EW", ew)
    plt.title("1984 Solar spectra")
    plt.ylabel("Intensity")
    plt.xlabel("Wavelength (AA)")
    plt.figtext(0.2, 0.14, "EW: " + str(np.round(ew, 8)))

    plt.plot(lambda_atlas, intensity_atlas)
    #plt.savefig(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\graphs\spectra\Neckel_atlas\Neckel_line"+str(int(np.round(line)))+".png")





    plt.show()