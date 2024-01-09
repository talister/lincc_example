import os

from astropy.tests.helper import assert_quantity_allclose
from importlib_resources import files
from synphot import units
from synphot.spectrum import BaseUnitlessSpectrum, SourceSpectrum, SpectralElement

from example_package import example_module
from example_package.utils import read_element  # , read_eso_spectra, percentage_difference


def test_greetings() -> None:
    """Verify the output of the `greetings` function"""
    output = example_module.greetings()
    assert output == "Hello from LINCC-Frameworks!"


def test_meaning() -> None:
    """Verify the output of the `meaning` function"""
    output = example_module.meaning()
    assert output == 42


class TestReadElement:
    eso_unit = units.u.photon / units.u.s / units.u.m**2 / units.u.um

    def test_read_element_as_sourcespectrum(self):
        test_fp = files("tests.example_package.data").joinpath("test_radiance.dat")

        print(f"test_fp= {test_fp}")
        element = read_element(test_fp.as_posix(), element_type="spectrum")

        assert isinstance(element, SourceSpectrum)
        assert element.waveset[0].unit == units.u.AA
        assert_quantity_allclose(element.waveset[0], 300 * units.u.nm)
        assert_quantity_allclose(element(element.waveset[0]), 16.4805 * units.PHOTLAM, rtol=1e-3)
