import os
import sys
import warnings

try:
    if sys.version_info >= (3, 9):
        # This exists in 3.8 but a different API
        import importlib.resources as pkg_resources
    else:
        raise ImportError
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources

from astropy import units as u
from astropy.utils.exceptions import AstropyUserWarning
from synphot import SourceSpectrum, SpectralElement, specio, units
from synphot.spectrum import BaseUnitlessSpectrum, Empirical1D


def read_element(filtername_or_filename, element_type="element", wave_units=u.nm, flux_units=units.THROUGHPUT):
    """Generic reader for optical elements, filters and others.
    The passed <filtername_or_filename> is first looked up in the `Conf` mapping
    dictionary for a match; if there is no match, it is assumed to be a filename
    or location specifier. These can be of 3 types:
    1. LCO Imaging Lab scanned optical element CSV files (contain 'LCO_' and '.csv' in the filename)
    2. SVO filter service references (contain 'http://svo')
    3. Local files (either ASCII or FITS)

    Checking on the read wavelengths is performed for local files:
    * if the first wavelength value is <100nm and the user didn't override the
    units through [wave_units], the wavelengths are assumed to be in, and are
    converted to, microns.
    * if the first wavelength value is >3000nm and the user didn't override the
    units through [wave_units], the wavelengths are assumed to be in, and are
    converted to, angstroms.
    """

    element_type = element_type.lower()
    filename = None  # conf.mapping.get(filtername_or_filename, None)
    if filename is None:
        filename = filtername_or_filename
    else:
        filename = filename()
        element_type = "spectral_element"
    if "LCO_" in filename.upper() and ".csv" in filename.lower():
        file_path = pkg_resources.files("example_package.data").joinpath(os.path.expandvars(filename))
        source = "LCO iLab format"
        header, wavelengths, throughput = read_lco_filter_csv(file_path)
    elif "http://svo" in filename.lower():
        source = "SVO filter service"
        header, wavelengths, throughput = specio.read_remote_spec(filename, wave_unit=u.AA, flux_unit=units.THROUGHPUT)
    else:
        source = "local file"
        file_path = os.path.expandvars(filename)
        if not os.path.exists(file_path):
            file_path = str(pkg_resources.files("example_package.data").joinpath(filename))
        warnings.simplefilter("ignore", category=AstropyUserWarning)
        if filename.lower().endswith("fits") or filename.lower().endswith("fit"):
            wave_col = "lam"
            flux_col = "trans"
            if element_type == "radiance":
                flux_col = "flux"
            try:
                header, wavelengths, throughput = specio.read_spec(
                    file_path, wave_col=wave_col, flux_col=flux_col, wave_unit=u.nm, flux_unit=flux_units
                )
            except KeyError:
                # ESO-SM01 format; different column name for transmission and micron vs nm
                header, wavelengths, throughput = specio.read_spec(
                    file_path, wave_col="lam", flux_col="flux", wave_unit=u.micron, flux_unit=flux_units
                )
        else:
            header, wavelengths, throughput = specio.read_ascii_spec(
                file_path, wave_unit=wave_units, flux_unit=flux_units
            )
        if wavelengths[0].value < 100.0 and wave_units == u.nm:
            # Small values seen, Convert to microns
            wavelengths = wavelengths.value * u.micron
        elif wavelengths[0].value > 3000.0 and wave_units == u.nm:
            # Large values seen, Convert to angstroms
            wavelengths = wavelengths.value * u.AA
        if element_type != "spectrum" and element_type != "radiance" and throughput.mean() > 1.5:
            # Test for mean throughput is above 1 to catch case where a throughput
            # fudge may be in the range ~1 to a few e.g. ESO Omegacam optics fudge
            # which goes to 3.2 and averages out to ~1.4
            throughput /= 100.0
            header["notes"] = "Divided by 100.0 to convert from percentage"
    header["source"] = source
    header["filename"] = filename
    if element_type == "spectrum" or element_type == "radiance":
        # SourceSpectrum can't use the default units.THROUGHPUT so we need to
        # change to an assumed units.PHOTLAM (u.photon / (u.cm**2 * u.s * u.AA))
        # if nothing was passed by the user
        if flux_units == units.THROUGHPUT:
            if element_type == "radiance":
                # Default units for ESO skycalc output (minus the arcsec^2)
                throughput = throughput.value * u.photon / u.s / u.m**2 / u.um
            else:
                throughput = throughput.value * units.PHOTLAM
        element = SourceSpectrum(
            Empirical1D, points=wavelengths, lookup_table=throughput, keep_neg=False, meta={"header": header}
        )
    elif element_type == "spectral_element":
        element = SpectralElement(
            Empirical1D, points=wavelengths, lookup_table=throughput, keep_neg=False, meta={"header": header}
        )
    else:
        element = BaseUnitlessSpectrum(
            Empirical1D, points=wavelengths, lookup_table=throughput, keep_neg=False, meta={"header": header}
        )

    return element
