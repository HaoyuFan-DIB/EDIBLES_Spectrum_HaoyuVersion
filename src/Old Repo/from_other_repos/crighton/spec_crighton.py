# Stuff to deal with a qso spectrum.

from itertools import izip
import copy
import os, pdb
from math import sqrt
from pprint import pformat

import numpy as np

import matplotlib.pyplot as pl
from matplotlib.ticker import FormatStrFormatter, MultipleLocator

from utilities import nan2num, between
from convolve import convolve_psf
from io import readtxt, readtabfits
from plot import axvlines, axvfill, puttext

PATH = os.path.abspath(os.path.dirname(__file__))
from constants import c_kms

debug = False
def getwave(hd):
    """ Given a fits header, get the wavelength solution.
    """
    dv=None
    if hd.has_key('CDELT1'):
        dw = hd['CDELT1']
    else:
        dw = hd['CD1_1']
    CRVAL = hd['CRVAL1']
    CRPIX = hd['CRPIX1']
    # wavelength of pixel 1
    wstart = CRVAL + (1 - CRPIX) * dw
    # check if it's log-linear scale (heuristic)
    if CRVAL < 10:
        wstart = 10**wstart
        dv = c_kms * (1. - 1. / 10. ** -dw)
        print 'constant dv = %.3f km/s (assume CRVAL1 in log(Angstroms))' % dv

    npts = hd['NAXIS1']
    return _make_wa_scale(wstart, dw, npts, constantdv=dv)

def find_wa_edges(wa):
    """ Given wavelength bin centres, find the edges of wavelengh
    bins.

    Examples
    --------
    >>> print find_wa_edges([1, 2.1, 3.3, 4.6])
    [ 0.45  1.55  2.7   3.95  5.25]
    """
    wa = np.asarray(wa)
    edges = wa[:-1] + 0.5 * (wa[1:] - wa[:-1])
    edges = [2*wa[0] - edges[0]] + edges.tolist() + [2*wa[-1] - edges[-1]]
    return np.array(edges)

def _make_wa_scale(wstart, dw, npts, constantdv=False, verbose=False):
    """ Generates a wavelength scale from the wstart, dw, and npts
    values.

    >>> wa = _make_wa_scale(40, 1, 5)
    >>> np.allclose(wa, [40., 41., 42., 43., 44.])
    True
    >>> wa = _make_wa_scale(3.5, 1e-3, 5, constantdv=True)
    >>> np.allclose(wa, [3162.278,  3169.568,  3176.874,  3184.198, 3191.537])
    True
    """
    if constantdv:
        if verbose:  print '_make_wa_scale(): Using log-linear scale'
        wa = 10**(wstart + np.arange(npts, dtype=float) * dw)
    else:
        if verbose:  print '_make_wa_scale(): Using linear scale'
        wa = wstart + np.arange(npts, dtype=float) * dw
    return wa

class Spectrum(object):
    """ A class to hold information about a spectrum.
    
    wa:       array of wavelength values (overrides all wavelength keywords)
    fl:       array of flux values
    er:       array of error (same units as flux) values
    co:       array of continuum values
    dw:       Wavelength difference between adjacent pixel centres
    dv:       Velocity difference (km/s) (overrides wdelt keyword)
    fwhm:     instrumental FWHM in km/s
    filename: filename of spectrum

    If enough information is given, the wavelength scale will be
    generated.  Note that there is no error check if you give
    conflicting wavelength scale information in the keywords!  In this
    case certain combinations of keywords take precendence.  See the
    code comments for details.

    Notes for FITS header:
    
    wstart = CRVAL - (CRPIX - 1.0) * CDELT,  dw = CDELT

    Conversion between velocity width and log-linear pixel width:
    
    dv / c_kms = 1 - 10**(-dw)

    Examples
    --------
    >>> sp = Spectrum(wstart=4000, dw=1, npts=500)
    >>> np.allclose(sp.fl, 0)
    True
    >>> sp = Spectrum(wstart=4000, dv=60, npts=500)
    >>> np.allclose(sp.dw, 8.692773e-05)
    True
    >>> sp = Spectrum(wstart=4000, wend=4400, npts=500)
    >>> np.allclose(sp.dw, 0.80160)
    True
    >>> wa = np.linspace(4000, 5000, 500)
    >>> fl = np.ones(len(wa))
    >>> sp = Spectrum(wa=wa, fl=fl)
    >>> np.allclose(sp.dw, 2.004008)
    True
    >>> sp = Spectrum(CRVAL=4000, CRPIX=1, CDELT=1, fl=np.ones(500))
    >>> print np.allclose(sp.dw, 1), np.allclose(sp.wa[0],4000)
    True True
    """
    def __init__(self,
                 dw=None, dv=None, wstart=None, wend=None, npts=None,
                 CRVAL=None, CRPIX=None, CDELT=None,
                 wa=None, fl=None, er=None, co=None,
                 fwhm=None, filename=None):
        """ Create the wavelength scale and initialise attributes."""
        if fl is not None:
            fl = np.asarray(fl)
            fl[np.isinf(fl)] = np.nan
            self.fl = fl
            npts = len(fl)
        if er is not None:
            er = np.asarray(er)
            # replace bad values with NaN
            er[np.isinf(er)|(er<=0.)] = np.nan
            self.er = er
            npts = len(er)
        if co is not None:
            co = np.asarray(co)
            co[np.isinf(co)] = np.nan
            self.co = co
            npts = len(co)

        # Check whether we need to make a wavelength scale.
        makescale = True

        if dv is not None:
            dw = np.log10(1. / (1. - dv / c_kms))

        if None not in (CRVAL, CRPIX, CDELT) :
            wstart = CRVAL - (CRPIX - 1.0) * CDELT
            dw = CDELT
            # check if it's log-linear scale (heuristic)
            if CRVAL < 10:
                wstart = 10**wstart
                dv = c_kms * (1. - 1. / 10. ** -dw)

        if wa is not None:
            wa = np.asarray(wa, float)
            npts = len(wa)
            makescale = False
        elif None not in (wstart, dw, npts):
            if dv is not None:
                wstart = np.log10(wstart)
        elif None not in (wstart, wend, dw):
            if dv is not None:
                wstart, wend = np.log10([wstart,wend])
            # make sure the scale is the same or bigger than the
            # requested wavelength range
            npts = int(np.ceil((wend - wstart) / float(dw)))
        elif None not in (wstart, wend, npts):
            # Make a linear wavelength scale
            dw = (wend - wstart) / (npts - 1.0)
        elif None not in (wend, dw, npts):
            raise ValueError('Please specify wstart instead of wend')
        else:
            raise ValueError('Not enough info to make a wavelength scale!')

        if makescale:
            if debug: print 'making wav scale,', wstart, dw, npts, bool(dv)
            wa = _make_wa_scale(wstart, dw, npts, constantdv=bool(dv))
        else:
            # check whether wavelength scale is linear or log-linear
            # (constant velocity)
            diff = wa[1:] - wa[:-1]
            if np.allclose(diff, diff[0]):
                dw = np.median(diff)
            else:
                diff = np.log10(wa[1:]) - np.log10(wa[:-1])
                if np.allclose(diff, diff[0]):
                    dw = np.median(diff)
                    dv = c_kms * (1. - 1. / 10. ** dw)

        # assign remaining attributes
        if fl is None:
            self.fl = np.zeros(npts)
        if er is None:
            self.er = np.empty(npts) * np.nan  # error (one sig)
        if co is None:
            self.co = np.empty(npts) * np.nan

        self.fwhm = fwhm
        self.dw = dw
        self.dv = dv
        self.filename = filename
        self.wa = wa

    def __getitem__(self, a):
        """
        >>> sp = Spectrum(wstart=4000, dw=1, npts=500)
        >>> sp1 = sp[sp.wa < 4100]
        """
        sp = self.view()
        sp.fl = self.fl[a]
        sp.wa = self.wa[a]
        sp.er = self.er[a]
        sp.co = self.co[a]
        return sp

    def __getslice__(self, i, j):
        """
        >>> sp = Spectrum(wstart=4000, dw=1, npts=500)
        >>> i,j = sp.wa.searchsorted([4050,4100])
        >>> sp1 = sp[i:j]
        """
        sp = self.view()
        sp.fl = self.fl[i:j]
        sp.wa = self.wa[i:j]
        sp.er = self.er[i:j]
        sp.co = self.co[i:j]
        return sp

    def __repr__(self):
        return 'Spectrum(wa, fl, er, co, dw, dv, fwhm, filename)'

    def __str__(self):
        temp = ("%s = %s" % (attr, pformat(self.__dict__[attr])) for attr
                in 'wa fl er co dw dv fwhm filename'.split())
        return 'Spectrum(\n%s)' % ',\n'.join(temp)

    def view(self):
        """ Creates a copy of the spectrum, without creating new copys
        of the wa,fl,er,co arrays.
        """
        return Spectrum(wa=self.wa, fl=self.fl, er=self.er, co=self.co)

    def multiply(self, val):
        """
        >>> sp = Spectrum(wstart=4000, dw=1, npts=500, fl=np.ones(500))
        >>> sp.multiply(2)
        >>> np.allclose(sp.fl, 2 * np.ones(500))
        True
        """
        self.fl *= val
        self.er *= val
        self.co *= val

    def velplot(self, *args, **kwargs):
        """ Method version of spec.velplot. See that function
        for documentation and a description of function arguments.
        """
        if np.all(np.isnan(self.co)):
            raise ValueError('A continuum must be defined!')
        return velplot(self.wa, self.fl, self.er, self.co, *args, **kwargs)

    def plot(self, ax=None, show=True, yperc=0.98, alpha=0.8,
             linewidth=1., linestyle='steps-mid',
             flcolor='blue', cocolor='red'):
        """ Plots spectrum.

        Returns the matplotlib artists that represent the flux, error
        and continuum curves.
        
        >>> sp = read('data/spSpec-52017-0516-139.fit.gz')
        >>> fig = pl.figure()
        >>> lines = sp.plot(show=False)
        """
        f,e,c,w = self.fl, self.er, self.co, self.wa
        return plot(w, f, e, c, ax=ax, show=show, yperc=yperc, alpha=alpha,
                    linewidth=linewidth, linestyle=linestyle,
                    flcolor=flcolor, cocolor=cocolor)

    def stats(self, wa1, wa2, show=False):
        """Calculates statistics (mean, standard deviation (i.e. RMS), mean
        error, etc) of the flux between two wavelength points.

        Returns
        -------
        mean flux, RMS of flux, mean error, SNR:
           SNR = (mean flux / RMS)
        
        >>> wa = np.linspace(10,20,10)
        >>> np.random.seed(77)
        >>> fl = np.random.randn(len(wa)) + 10
        >>> er = np.ones(len(wa))
        >>> sp = Spectrum(wa=wa, fl=fl, er=er)
        >>> data = sp.stats(11, 18)
        >>> np.allclose(data, [9.66873, 0.98909, 1.0, 9.77542])
        True
        """
        i,j = self.wa.searchsorted([wa1, wa2])
        fl = self.fl[i:j]
        er = self.er[i:j]
        good = (er > 0) & ~np.isnan(fl)
        if len(good.nonzero()[0]) == 0:
            print('No good data in this range!')
            return np.nan, np.nan, np.nan, np.nan
        fl = fl[good]
        er = er[good]
        mfl = fl.mean()
        std = fl.std()
        mer = er.mean()
        snr = mfl / std
        if show:
            print 'mean %g, std %g, er %g, snr %g' % (mfl, std, mer, snr)
        return mfl, std, mer, snr

    def rebin(self, **kwargs):
        """ Class method version of spec.rebin() """
        return rebin(self.wa, self.fl, self.er, **kwargs)

    def rebin_simple(self, n):
        """ Class method version of spec.rebin_simple()."""
        return rebin_simple(self.wa, self.fl, self.er, self.co, n)
    
    def write(self, filename, header=None, overwrite=False):
        """ Writes out a spectrum, as ascii for now - wavelength,
        flux, error, continuum.

        Keyword overwrite can be True or False.

        header is a string to be written to the file before the
        spectrum. A special case is header='RESVEL', which means the
        instrumental fwhm in km/s will be written on the first line
        (VPFIT style).
        """
        if os.path.lexists(filename) and not overwrite:
            c = raw_input('File %s exists - overwrite? (y) or n: ' % filename)
            if c != '':
                if c.strip().lower()[0] == 'n':
                    print 'returning without writing anything...'
                    return
        fh = open(filename, 'w')
        if header is not None:
            if header == 'RESVEL':
                if self.fwhm is None:
                    raise ValueError('Instrumental fwhm is not set!')
                fh.write('RESVEL %.2f' % self.fwhm)
            else:
                fh.write(header)
        fl = np.nan_to_num(self.fl)
        er = np.nan_to_num(self.er)
        if np.all(np.isnan(self.co)):
            for w,f,e in izip(self.wa, fl, er):
                fh.write("% .12g % #12.8g % #12.8g\n" % (w,f,e))
        else:
            co = np.nan_to_num(self.co)
            for w,f,e,c in izip(self.wa, fl, er, co):
                fh.write("% .12g % #12.8g % #12.8g % #12.8g\n" % (w,f,e,c))
        fh.close()
        if self.filename is None:
            self.filename = filename

def read(filename, comment='#', debug=False):
    """
    Reads in QSO spectrum given a filename.  Returns a Spectrum class
    object.

    comment = '#' :  String that marks beginning of comment line,
                     only used when reading in ascii files

    Examples
    --------
    >>> sp = read('data/HE0940m1050m.txt.gz')
    >>> sp = read('data/spSpec-52017-0516-139.fit.gz')
    >>> sp = read('data/runA_h1_100.txt.gz')
    >>> np.allclose(sp.fwhm, 6.6)
    True
    """
    if filename.endswith('gz'):
        import gzip
        fh = gzip.open(filename)
    else:
        fh = open(filename,'r')
    test = fh.next()
    fh.close()

    if len(test) < 9 or test[8] != '=': 
    # Then probably not a fits file
        fwhm = None
        skip = 0
        test = test.split()
        try:
            name = test[0].strip()
        except IndexError:
            pass
        else:
            if name.upper() == 'RESVEL':
                fwhm = float(test[1])
                skip = 1
        try:         # uves_popler .dat file            
            wa,fl,er,co = np.loadtxt(filename, usecols=(0,1,2,4), unpack=True,
                                     comments=comment, skiprows=skip)
        except IndexError:
            try:
                wa,fl,er,co = np.loadtxt(filename, usecols=(0,1,2,3),
                                         unpack=True, comments=comment,
                                         skiprows=skip)
            except IndexError:
                try:
                    wa,fl,er = np.loadtxt(filename, usecols=(0,1,2),
                                          unpack=True, comments=comment,
                                          skiprows=skip)
                except IndexError:
                    wa,fl = np.loadtxt(filename, usecols=(0,1),
                                       unpack=True, comments=comment,
                                       skiprows=skip)
                    er = find_err(fl, find_cont(fl))
                co = None
        else:
            # heuristic to check for Jill Bechtold's FOS spectra
            if filename.endswith('.XY'):
                wa,fl,er,co = np.loadtxt(
                    filename, usecols=(0,1,2,3), unpack=True, comments=comment,
                    skiprows=skip)
            else:
                fl *= co
                er *= co

        if wa[0] > wa[-1]:
            wa = wa[::-1];  fl = fl[::-1];
            if er is not None:  er = er[::-1]
            if co is not None:  co = co[::-1]
        sp = Spectrum(wa=wa, fl=fl, er=er, co=co, filename=filename, fwhm=fwhm)
        return sp

    # Otherwise assume fits file
    import pyfits
    f = pyfits.open(filename)
    hd = f[0].header

    # try record array
    try:
        data = f[1].data
    except IndexError:
        pass
    else:
        good = False
        names = data.dtype.names
        if 'wa' in  names and 'fl' in names:
            wa = data.wa
            fl = data.fl
            good = True
        elif 'wavelength' in names and 'flux' in names:
            wa = data.wavelength
            fl = data.flux
            good = True
        if good:
            er = np.ones_like(fl)
            try:
                er = data.er
            except AttributeError:
                pass
            co = np.empty_like(fl) * np.nan
            try:
                co = data.co
            except AttributeError:
                pass
            return Spectrum(wa=wa, fl=fl, er=er, co=co, filename=filename)

    ##################################################################
    #  First generate the wavelength scale.  Look for CTYPE, CRVAL,
    #  CDELT header cards.  Then read in flux values, and look for
    #  cont and error axes/files.
    ##################################################################

    #naxis1 = hd['NAXIS1']     # 1st axis length (no. data points)
    # pixel stepsize
    if hd.has_key('CDELT1'):
        cdelt = hd['CDELT1']
    elif hd.has_key('CD1_1'):
        cdelt = hd['CD1_1']
    else:
        # Read Songaila's spectra
        wa = f[0].data[0]
        fl = f[0].data[1]
        npts = len(fl)
        try:
            er = f[0].data[2]
            if len(er[er > 0] < 0.5 * npts):
                er = f[0].data[3]
        except IndexError:
            i = int(npts * 0.75) 
            er = np.ones(npts) * np.std(fl[i-50:i+50])
        f.close()
        return Spectrum(wa=wa, fl=fl, er=er, filename=filename)

    ##########################################################
    # Check if SDSS spectrum
    ##########################################################
    if hd.has_key('TELESCOP'):
        if hd['TELESCOP'] == 'SDSS 2.5-M':  # then Sloan spectrum
            data = f[0].data
            fl = data[0]
            er = data[2]
            f.close()
            return Spectrum(fl=fl, er=er, filename=filename, CDELT=cdelt,
                            CRVAL=hd['CRVAL1'], CRPIX=hd['CRPIX1'])

    ##########################################################
    # Check if HIRES spectrum
    ##########################################################
    if hd.has_key('INSTRUME'):   #  Check if Keck spectrum
        if hd['INSTRUME'].startswith('HIRES'):
            if debug:  print 'Looks like Makee output format'
            fl = f[0].data       # Flux
            f.close()
            errname = filename[0:filename.rfind('.fits')] + 'e.fits'
            try:
                er = pyfits.getdata(errname)
            except IOError:
                er = np.ones(len(fl))
            return Spectrum(fl=fl, er=er, filename=filename, CDELT=cdelt,
                            CRVAL=hd['CRVAL1'], CRPIX=hd['CRPIX1'])

    ##########################################################
    # Check if UVES_popler output
    ##########################################################
    history = hd.get_history()
    for row in history:
        if 'UVES POst Pipeline Echelle Reduction' in row:
            co = f[0].data[3]
            fl = f[0].data[0] * co #  Flux
            er = f[0].data[1] * co
            f.close()
            return Spectrum(fl=fl, er=er, co=co, filename=filename, CDELT=cdelt,
                            CRVAL=hd['CRVAL1'], CRPIX=hd['CRPIX1'])

    data = f[0].data
    fl = data[0]
    er = data[2]
    f.close()
    if hd.has_key('CRPIX1'):
        crpix = hd['CRPIX1']
    else:
        crpix = 1
    return Spectrum(fl=fl, er=er, filename=filename, CDELT=cdelt,
                    CRVAL=hd['CRVAL1'], CRPIX=crpix)
    #raise Exception('Unknown file format')

def rebin_simple(wa, fl, er, co, n):
    """ Bins up the spectrum by averaging the values of every n
    pixels. Not very accurate, but much faster than rebin().

    TODO: re-write using masked arrays to preserve masked values.
    """ 
    remain = -(len(wa) % n) or None
    wa = wa[:remain].reshape(-1, n)
    fl = fl[:remain].reshape(-1, n)
    er = er[:remain].reshape(-1, n)
    co = co[:remain].reshape(-1, n)
    n = float(n)
    wa = np.nansum(wa, axis=1) / n 
    co = np.nansum(co, axis=1) / n 
    er = np.nansum(er, axis=1) / n / sqrt(n) 
    fl = np.nansum(fl, axis=1) / n 

    return Spectrum(wa=wa, fl=fl, er=er, co=co)

def rebin(wav, fl, er, **kwargs):
    """ Rebins spectrum to a new wavelength scale generated using the
    keyword parameters.

    Returns the rebinned spectrum.

    Accepts the same keywords as Spectrum.__init__() (see that
    docstring for a description of those keywords)

    Note - will probably get the flux and errors for the first and
    last pixel of the rebinned spectrum wrong.

    General pointers about rebinning if you care about errors in the
    rebinned values:

    1. Don't rebin to a smaller bin size.
    2. Be aware when you rebin you introduce correlations between
       neighbouring points and between their errors.
    3. Rebin as few times as possible.

    >>> wa = np.linspace(11,20,10)
    >>> np.random.seed(77)
    >>> fl = np.random.randn(len(wa)) + 10
    >>> er = np.ones(len(wa))
    >>> er[3] = 0
    >>> sp = Spectrum(wa=wa, fl=fl, er=er)
    >>> rsp = sp.rebin(wstart=sp.wa[0], wend=sp.wa[-1], dw=1.5*sp.dw)
    >>> np.allclose(rsp.fl, [10.3119,10.0409,9.94336,9.2459,9.5919,9.9963])
    True
    >>> np.allclose(rsp.er, [0.89443,0.81650,1.4142,0.81650,0.81650,0.81650])
    True
    >>> np.allclose(rsp.wa, [11., 12.5, 14., 15.5, 17., 18.5])
    True
    """
    # Note: 0 suffix indicates the old spectrum, 1 the rebinned spectrum.
    colors= 'brgy'
    debug = kwargs.pop('debug', False)
    
    # Create rebinned spectrum wavelength scale
    sp1 = Spectrum(**kwargs)
    # find pixel edges, used when rebinning
    edges0 = find_wa_edges(wav)
    edges1 = find_wa_edges(sp1.wa)
    if debug:
        pl.clf()
        x0,x1 = edges1[0:2]
        yh, = pl.bar(x0, 0, width=(x1-x0),color='gray',
                    linestyle='dotted',alpha=0.3)
    widths0 = edges0[1:] - edges0[:-1]
    npts0 = len(wav)
    npts1 = len(sp1.wa)
    df = 0.
    de2 = 0.
    npix = 0    # number of old pixels contributing to rebinned pixel,
    j = 0                # index of rebinned array
    i = 0                # index of old array

    # sanity check
    if edges0[-1] < edges1[0] or edges1[-1] < edges0[0]:
        raise ValueError('Wavelength scales do not overlap!')
    
    # find the first contributing old pixel to the rebinned spectrum
    if edges0[i+1] < edges1[0]:
        # Old wa scale extends lower than the rebinned scale. Find the
        # first old pixel that overlaps with rebinned scale.
        while edges0[i+1] < edges1[0]:
            i += 1
        i -= 1
    elif edges0[0] > edges1[j+1]:
        # New rebinned wa scale extends lower than the old scale. Find
        # the first rebinned pixel that overlaps with the old spectrum
        while edges0[0] > edges1[j+1]:
            sp1.fl[j] = np.nan
            sp1.er[j] = np.nan
            j += 1
        j -= 1
    lo0 = edges0[i]      # low edge of contr. (sub-)pixel in old scale
    while True:
        hi0 = edges0[i+1]  # upper edge of contr. (sub-)pixel in old scale
        hi1 = edges1[j+1]  # upper edge of jth pixel in rebinned scale

        if hi0 < hi1:
            if er[i] > 0:
                dpix = (hi0 - lo0) / widths0[i]
                df += fl[i] * dpix
                # We don't square dpix below, since this causes an
                # artificial variation in the rebinned errors depending on
                # how the old wav bins are divided up into the rebinned
                # wav bins.
                #
                # i.e. 0.25**2 + 0.75**2 != 0.5**2 + 0.5**2 != 1**2
                de2 += er[i]**2 * dpix
                npix += dpix
            if debug:
                yh.set_height(df/npix)
                c0 = colors[i % len(colors)]
                pl.bar(lo0, fl[i], width=hi0-lo0, color=c0, alpha=0.3)
                pl.text(lo0, fl[i], 'lo0')
                pl.text(hi0, fl[i], 'hi0')
                pl.text(hi1, fl[i], 'hi1')
                raw_input('enter...')
            lo0 = hi0
            i += 1
            if i == npts0:  break
        else:
            # We have all old pixel flux values that contribute to the
            # new pixel; append the new flux value and move to the
            # next new pixel.
            if er[i] > 0:
                dpix = (hi1 - lo0) / widths0[i]
                df += fl[i] * dpix
                de2 += er[i]**2 * dpix
                npix += dpix
            if debug:
                yh.set_height(df/npix)
                c0 = colors[i % len(colors)]
                pl.bar(lo0,  fl[i], width=hi1-lo0, color=c0, alpha=0.3)
                pl.text(lo0, fl[i], 'lo0')
                pl.text(hi0, fl[i], 'hi0')
                pl.text(hi1, fl[i], 'hi1')
                raw_input('df, de2, npix: %s %s %s   enter...' %
                          (df, de2, npix))
            if npix > 0:
                # find total flux and error, then divide by number of
                # pixels (i.e. conserve flux density).
                sp1.fl[j] = df / npix
                sp1.er[j] = sqrt(de2) / npix
            else:
                sp1.fl[j] = np.nan
                sp1.er[j] = np.nan
            df = 0.
            de2 = 0.
            npix = 0.
            lo0 = hi1
            j += 1
            if j == npts1:  break
            if debug:
                x0,x1 = edges1[j:j+2]
                yh, = pl.bar(x0, 0, width=x1-x0, color='gray',
                       linestyle='dotted', alpha=0.3)
                raw_input('enter...')

    return sp1

def combine(spectra, cliphi=None, cliplo=None, verbose=False):
    """ Combine spectra pixel by pixel, weighting by the inverse variance
    of each pixel.  Clip high sigma values by sigma times clip values
    Returns the combined spectrum.

    If the wavelength scales of the input spectra differ, combine()
    will rebin the spectra to a common linear (not log-linear)
    wavelength scale, with pixel width equal to the largest pixel
    width in the input spectra. If this is not what you want, rebin
    the spectra by hand with rebin() before using combine().

    >>> wa = np.linspace(11,20,10)
    >>> np.random.seed(77)
    >>> fl1 = np.random.randn(len(wa)) + 10
    >>> fl2 = np.random.randn(len(wa)) + 10
    >>> er1 = np.ones(len(wa))
    >>> er1[7] = -1
    >>> er2 = np.ones(len(wa))
    >>> er2[7] = -1
    >>> er2[3] = -1
    >>> sp1 = Spectrum(wa=wa, fl=fl1, er=er1)
    >>> sp2 = Spectrum(wa=wa, fl=fl2, er=er2)
    >>> csp = combine([sp1,sp2])
    >>> np.allclose(csp.er[~np.isnan(csp.er)], [0.70711, 0.70711, 0.70711, \
              1., 0.70711, 0.70711, 0.70711, 0.70711, 0.70711])
    True
    >>> np.allclose(csp.fl[~np.isnan(csp.fl)], [10.2652, 10.9127, 9.1856, \
              10.4078, 9.9137, 9.8614, 9.4947, 10.7472, 9.5573])
    True
    """

    def clip(cliphi, cliplo, s_rebinned):
        # clip the rebinned input spectra

        # find pixels where we can clip: where we have at least three
        # good contributing values.
        goodpix = np.zeros(len(s_rebinned[0].wa))
        for s in s_rebinned:
            goodpix += (s.er > 0).astype(int)
        canclip = goodpix > 2
        # find median values
        medfl = np.median([s.fl[canclip] for s in s_rebinned], axis=0)
        nclipped = 0
        for i,s in enumerate(s_rebinned):
            fl = s.fl[canclip]
            er = s.er[canclip]
            diff = (fl - medfl) / er
            if cliphi is not None:
                badpix = diff > cliphi
                s_rebinned[i].er[canclip][badpix] = np.nan
                nclipped += len(badpix.nonzero()[0])
            if cliplo is not None:
                badpix = diff < -cliplo
                s_rebinned[i].er[canclip][badpix] = np.nan
                nclipped += len(badpix.nonzero()[0])
        if debug: print nclipped, 'pixels clipped across all input spectra'
        return nclipped

    nspectra = len(spectra)
    if verbose:
        print '%s spectra to combine' % nspectra
    if nspectra < 2:
        raise Exception('Need at least 2 spectra to combine.')

    if cliphi is not None and nspectra < 3:  cliphi = None
    if cliplo is not None and nspectra < 3:  cliplo = None

    # Check if wavescales are the same:
    spec0 = spectra[0]
    wa = spec0.wa
    npts = len(wa)
    needrebin = True
    for sp in spectra:
        if len(sp.wa) != npts:
            if verbose: print 'Rebin required'
            break
        if (np.abs(sp.wa - wa) / wa[0]).max() > 1e-8:
            if verbose:
                print (np.abs(sp.wa - wa) / wa[0]).max(), 'Rebin required'
            break
    else:
        needrebin = False
        if verbose:  print 'No rebin required'

    # interpolate over 1 sigma error arrays
    
    if needrebin:
        # Make wavelength scale for combined spectrum.  Only linear for now.
        wstart = min(sp.wa[0] for sp in spectra)
        wend = max(sp.wa[-1] for sp in spectra)
        # Choose largest wavelength bin size of old spectra.
        if verbose:  print 'finding new bin size'
        maxwidth = max((sp.wa[1:] - sp.wa[:-1]).max() for sp in spectra)
        npts = int(np.ceil((wend - wstart) / maxwidth))      # round up
        # rebin spectra to combined wavelength scale
        if verbose:  print 'Rebinning spectra'
        s_rebinned = [s.rebin(wstart=wstart, npts=npts, dw=maxwidth)
                      for s in spectra]
        combined = Spectrum(wstart=wstart, npts=npts, dw=maxwidth)
        if verbose:
            print ('New wavelength scale wstart=%s, wend=%s, npts=%s, dw=%s'
                   % (wstart, combined.wa[-1], npts, maxwidth))
    else:
        combined = Spectrum(wa=spec0.wa)
        s_rebinned = copy.deepcopy(spectra)

    # sigma clipping, if requested
    if cliphi is not None or cliplo is not None:
        clip(cliphi, cliplo, s_rebinned)
        # repeat, clipping to 4 sigma this time
        #npixclipped = clip(4.,4.,s_rebinned)

    # Now add the spectra
    for i in xrange(len(combined.wa)):
        wtot = fltot = ertot = 0.
        npix = 0            # num of old spectrum pixels contributing to new
        for s in s_rebinned:
            # if not a sensible flux value, skip to the next pixel
            if s.er[i] > 0:
                npix += 1
                # Weighted mean (weight by inverse variance)
                variance = s.er[i] ** 2
                w = 1. / variance
                fltot += s.fl[i] * w
                ertot += (s.er[i] * w)**2
                wtot += w
        if npix > 0:
            combined.fl[i] = fltot / wtot
            combined.er[i] = np.sqrt(ertot) / wtot
        else:
            combined.fl[i] = np.nan
            combined.er[i] = np.nan

        #contributing.fl[i] = npix_contrib
    return combined

def cr_reject(flux, error, nsigma=15.0, npix=2, verbose=False):
    """ Given flux and errors, rejects cosmic-ray type or dead
    pixels. These are defined as pixels that are more than
    nsigma*sigma above or below the median of the npixEL pixels on
    either side.

    Returns newflux,newerror where the rejected pixels have been
    replaced by the median value of npix to either side, and the
    error has been set to NaN.

    The default values work ok for S/N~20, Resolution=500 spectra.
    """
    if verbose:  print nsigma,npix
    flux,error = list(flux), list(error)  # make copies
    i1 = npix
    i2 = len(flux) - npix
    for i in range(i1, i2):
        # make a list of flux values used to find the median
        fl = flux[i-npix:i] + flux[i+1:i+1+npix]
        er = error[i-npix:i] + error[i+1:i+1+npix]
        fl = [f for f,e in zip(fl,er) if e > 0]
        er = [e for e in er if e > 0]
        medfl = np.median(fl)
        meder = np.median(er)
        if np.abs((flux[i] - medfl) / meder) > nsigma:
            flux[i] = medfl
            error[i] = np.nan
            if verbose:  print len(fl), len(er)

    return np.array(flux), np.array(error)

def cr_reject2(fl, er, nsig=10.0, fwhm=2, grow=1, debug=True):
    """ interpolate across features that have widths smaller than the
    expected fwhm resolution.

    fwhm: int, resolution fwhm in pixels
    fl: flux array
    er: error array

    returns
    fl,er: interpolated fl, er arrays
    """
    fl, er = (np.array(a, dtype=float) for a in (fl, er))
    # interpolate over bad pixels 
    fl1 = convolve_psf(fl, fwhm)
    ibad = np.where(np.abs(fl1 - fl) > nsig*er)[0]
    if debug: print len(ibad)
    extras1 = np.concatenate([ibad + 1 + i for i in range(grow)]) 
    extras2 = np.concatenate([ibad - 1 - i for i in range(grow)]) 
    ibad = np.union1d(ibad, np.union1d(extras1, extras2))
    ibad = ibad[(ibad > -1) & (ibad < len(fl))]
    igood = np.setdiff1d(np.arange(len(fl1)), ibad)
    fl[ibad] = np.interp(ibad, igood, fl[igood])
    er[ibad] = np.nan
    return fl,er

def scalemult(w0, f0, e0, w1, f1, e1, mask=None):
    """ find the constant to multipy f1 by so its median will match
    f0 where they overlap in wavelength.
    
    Errors are just used to identify bad pixels, if they are all > 0
    they are ignored.

    mask is optional, and of the form [(3500,4500), (4650,4680)] to mask
    wavelengths from 3500 to 4500 Ang, and 4650 to 4680 Ang, etc.
    """
    w0,f0,e0,w1,f1,e1 = map(np.asarray, [w0,f0,e0,w1,f1,e1])
    masked0 = np.zeros(len(w0), bool)
    masked1 = np.zeros(len(w0), bool)
    if mask is not None:
        for wmin,wmax in mask:
            masked0 |= (wmin < w0) & (w0 < wmax)  
            masked1 |= (wmin < w1) & (w1 < wmax)  

    wmin = max(w0.min(), w1.min())
    wmax = min(w0.max(), w1.max())
    
    good0 = (e0 > 0) & ~np.isnan(f0) & ~masked0 & (wmin < w0) & (w0 < wmax)
    good1 = (e1 > 0) & ~np.isnan(f1) & ~masked1 & (wmin < w1) & (w1 < wmax)
    if good0.sum() < 3 or good1.sum() < 3:
        raise ValueError('Too few good pixels to use for scaling')

    med0 = np.median(f0[good0])
    med1 = np.median(f1[good1])

    if not (med0 > 0) or not (med1 > 0):
        raise ValueError('bad medians:', med0, med1)

    return med0 / med1

def scale_overlap(w0, f0, e0, w1, f1, e1):
    """ Scale two spectra to match where they overlap. Assumes
    spectrum 0 covers a lower wavelength range than spectrum 1. If no
    good regions overlap, regions close to the overlap are searched
    for.

    Returns
    -------
    scale_factor : float
        Multiply spectrum 1 by scale_factor to match spectrum 0.
    """

   # find overlapping regions
    dtype = [('wa','f8'),('fl','f8'),('er','f8')]
    sp0 = np.rec.fromarrays([w0,f0,e0], dtype=dtype)
    sp1 = np.rec.fromarrays([w1,f1,e1], dtype=dtype)
    if sp0.wa.max() < sp1.wa.min():
        raise ValueError('No overlap!')
    # find first overlapping good pixel
    i = 0
    while np.isnan(sp1.er[i]):  i += 1
    i0min = sp0.wa.searchsorted(sp1.wa[i])
    i = -1
    while np.isnan(sp0.er[i]):  i -= 1
    i1max = sp1.wa.searchsorted(sp0.wa[i])
    while True:
        s0 = sp0[i0min:]
        s1 = sp1[:i1max]
        good0 = (s0.er > 0) & ~np.isnan(s0.fl)
        good1 = (s1.er > 0) & ~np.isnan(s1.fl)
        if good0.sum() == 0 or good1.sum() == 0:
            i0min = i0min - (len(sp0) - i0min)
            i1max = 2 * i1max
            if i0min < 0 or i1max > len(sp1)-1:
                raise ValueError('No good pixels to use for scaling')
            continue
        m0, m1 = np.median(s0.fl[good0]), np.median(s1.fl[good1])
        if m0 <= 0:
            print 'looping'
            i0min = max(0, i0min - (len(sp0) - i0min))
        elif m1 <= 0:
            print 'looping'
            i1max = min(len(sp1)-1, 2 * i1max)
        else:
            break
        
    if debug:  print m0,m1
    return m0 / m1

# these are from Vanden Berk et al 2001. Note some line IDs in the Lya
# forest are uncertain
REST = {'LyLim'   : 912.0    ,
        'CIIIa'   : 977.020  ,
        'Lyb'     : 1025.72  ,
        'ArI'     : 1065.10  ,
        'FeIII'   : 1117.26  ,
        'CIII*'   : 1175.35  ,
        'Lya'     : 1215.6701,
        'NV'      : 1240.14  ,
        'OI/SiII' : 1303.0   ,
        'CII'     : 1334.5323,
        'SiIV'    : 1400.    ,
        'CIV'     : 1550.    ,
        'CIIIb'   : 1908.73  ,
        'MgII'    : 2798.75  ,
        'OII'     : 3728.48 ,
        'Hg'      : 4341.7   ,
        'Hb'      : 4862.    ,
        'Ha'      : 6564.61  ,
        }
REST = np.rec.fromrecords(REST.items(), names='name,wa')
REST.sort(order='wa')

def plotlines(z, ax, atmos=None, lines=REST, labels=False, ls='dotted',
              color='k', trim=False, **kwargs):
    """ Draw vertical dotted lines showing expected positions of
    absorption and emission lines, given a redshift.

    atmos: list of wavelength pairs or True (None) Regions of
           atmospheric absorption to plot. If True, it uses an
           internal list of regions taken from aaomega spectra.

    lines: If given, it must be a record array with fields 'name' and 'wa'.
    

    Returns the mpl artists representing the lines.
    """
    lines = np.rec.fromrecords([(l['name'], l['wa']) for l in lines],
                           names='name,wa')
    autoscale = ax.get_autoscale_on()
    if autoscale:
        ax.set_autoscale_on(False)
    artists = []
    w0,w1 = ax.get_xlim()
    wa = lines.wa * (z+1)
    if trim:
        c0 = between(wa,w0,w1)
        wa = wa[c0]
        lines = lines[c0]
    artists.append(axvlines(wa, ax=ax, ls=ls, color=color, **kwargs))
    if labels:
        for i in range(3):
            for w,l in zip(wa[i::3], lines[i::3]):
                if not (w0 < w < w1) and trim:
                    continue
                #name = l.name + '%.2f' % l.wa
                name = l.name
                artists.append(puttext(
                    w, 0.7 + i*0.08, name, ax,
                    xcoord='data', alpha=1, fontsize=10,
                    rotation=90, ha='right'))
    if atmos:
        if atmos == True:
            atmos = None
        artists.append(plotatmos(ax, atmos=atmos))
    if autoscale:
        ax.set_autoscale_on(True)

    return artists

def plotatmos(ax, atmos=None, color='y'):
    """ Plot rough areas where atmospheric absorption is expected.
    """
    autoscale = ax.get_autoscale_on()
    if autoscale:
        ax.set_autoscale_on(False)
    artists = []
    if atmos is None:
        atmos = [(5570, 5590),
                 (5885, 5900),
                 (6275, 6325),
                 (7240, 7400),
                 (7580, 7690),
                 (7700, 8070),
                 (8240, 8550),
                 (8750, 9100),
                 (9300, 11000)]

    artists = axvfill(atmos, ax=ax, color=color, alpha=0.15)
    if autoscale:
        ax.set_autoscale_on(True)
    return artists

def plot(w, f=None, e=None, c=None, ax=None, show=True, yperc=0.98, alpha=0.8,
         linewidth=1., linestyle='steps-mid', flcolor='blue', cocolor='red'):
    """ Plots spectrum.

    Returns the matplotlib artists that represent the flux, error
    and continuum curves.

    Can also give a single argument that is a record array with fields
    wa, fl and optionally er, co.

    >>> sp = read('data/spSpec-52017-0516-139.fit.gz')
    >>> lines = plot(sp.wa, sp.fl, sp.er, sp.co, show=False)
    """
    witherr = True
    if f is None and e is None and c is None:
        rec = w
        w = rec.wa
        f = rec.fl
        try:
            e = rec.er
        except AttributeError:
            e = w * np.nan
            witherr = False
        try:
            c = rec.co
        except AttributeError:
            c = w * np.nan
    elif e is None:
        e = w * np.nan
        witherr = False
    elif c is None:
        c = w * np.nan

    if ax is None:
        fig = pl.figure(figsize=(10,5))
        fig.subplots_adjust(left=0.08, right=0.94)
        ax = pl.gca()

    artists = []
    # check we have something to plot...
    
    good = ~np.isnan(f)
    if witherr:
        good &= (e > 0)
    if np.any(good):
        artists.extend(ax.plot(w, f, lw=linewidth, color=flcolor,
                               alpha=alpha, ls=linestyle))
        if witherr:
            artists.extend(ax.plot(w, e, lw=linewidth, color=flcolor,
                                   alpha=alpha, ls='dashed'))
        artists.extend(ax.plot(w, c, lw=linewidth, color=cocolor,
                               alpha=alpha))
        # plotting limits
        f = f[good]
        ymax = 1.5 * np.percentile(f, 100*yperc)
        if witherr:
            e = e[good]
            ymin = min(-0.1 * np.median(f), -1.0 * np.median(e))
        else:
            ymin = -abs(0.1*ymax)

        if ymax < ymin:
            ymax = 2. * abs(ymin)
        ax.axis((w.min(), w.max(), ymin, ymax))

    ax.axhline(y=0., color='k', alpha=alpha)
    if show:
        pl.show()
    return artists


def velplot(wa, fl, er, co, transitions, z=None, vmin=-5000., vmax=5000.,
            fig=None, at=None, ticks=None, show=True,
            b=None, logN=None, resolution=None):
    """ Vertical stacked plots of expected positions of absorption
    lines at the given redshift on a velocity scale.

    Can also plot tickmarks.

    Can also plot transitions (if sample b's are given)

    Input
    -----
    wa,fl,er,co:
        Spectrum wavelength, flux, error and continuum.
    transitions: string of comma-separated values
        Each entry can be either an ion and approximate wavelength
        ('HI 1215, FeII 2600') or only an ion ('CIV, SiIV'). If only
        an ion is given, the two strongest transitions in atom.dat are
        used (useful for plotting doublets).
    z: float (None)
        Redshift.
    at: (None)
        An atom.dat table.
    vmin, vmax:  (-5000, 5000)
        Limits of the velocity plots in km/s.
    ticks: array of floats (None)
        The wavelengths of tickmarks to plot.
    b: sequence of floats (None)
        The velocity widths of generated lines. If b, logN and
        resolution are given, sample absorption lines are generated.
    logN: sequence of floats (None)
        log10 of column densities of generated lines (cm**-2).
    resolution: float (None)
        Resolution of the spectrum if generating lines

    Returns
    -------
    fig
        The matplotlib figure of the velocity plot.

    Examples
    --------
    >>> s = read('data/spSpec-52017-0516-139.fit.gz')
    >>> fig = velplot(s.wa,s.fl,s.er,s.co,'CIV 1548, CIV 1550, HI 1215', z=3.0, show=0)
    """

    from utilities import findtrans,readatom, calc_abs
    from matplotlib.ticker import NullLocator
    findabs = False
    if None not in (b, logN, resolution):
        findabs = True

    if at is None:
        at = readatom()

    if isinstance(transitions, basestring):
        transitions = transitions.split(',')
    
    n = len(transitions)
    zp1 = z + 1
    betamin = vmin / c_kms
    betamax = vmax / c_kms

    if fig is None:
        fig = pl.gcf()
    fig.clf()
    #fig.subplots_adjust(hspace=0.001)
    ax = fig.add_subplot(111)
    ax.set_autoscale_on(0)
    # plot top down, so we need reversed()
    artists_line = []
    artists_tick = []
    for i,tstring in enumerate(reversed(transitions)):
        offset = i*1.5
        ion, trans = findtrans(tstring, atomdat=at)
        watrans = trans['wav']
        obswa = watrans * zp1
        wmin = obswa * (1 + betamin)
        wmax = obswa * (1 + betamax)
        if i == 0:
            ax.set_title('z = %.4f, wa = %.2f' % (z,obswa))

        if ticks is not None:
            cond = (wmin < ticks) & (ticks < wmax)
            if np.any(cond):
                vel = (ticks[cond] / obswa - 1) * c_kms
                a = ax.vlines(vel, 1.1+offset, 1.35+offset, colors='k')
                artists_tick.append(a)

        cond = (wmin < wa) & (wa < wmax)
        #good = ~np.isnan(fl) & (er > 0) & ~np.isnan(co)
        
        if np.any(cond):
            nfl = fl[cond] / co[cond]
            ner = er[cond] / co[cond]
            # generate absorption from the ion
            nco = np.ones(len(nfl))
            if findabs:
                for tr in trans:
                    #print tstring, tr
                    nco *= calc_abs(wa[cond], at[ion], tr, zp1, resolution,
                                    logN=logN[i], b=b[i])
            vel = (wa[cond] / obswa - 1) * c_kms
            artists_line.extend(
                ax.plot(vel, nfl + offset, lw=1, alpha=0.7, ls='steps-mid') )
            artists_line.extend( 
                ax.plot(vel, ner + offset, lw=1, color='orange', alpha=0.7) )
            if findabs:
                artists_line.extend(
                    ax.plot(vel, nco + offset, 'r', alpha=0.7) )
            ax.axhline(offset, color='gray')

        ax.text(0.8 * vmax, 0.5 + offset, tstring)
    ax.axvline(0, color='k', linestyle='dotted')
    ax.set_xlim(vmin, vmax)
    ax.set_ylim(-0.5, len(transitions) * 1.5)
    ax.yaxis.set_major_locator(NullLocator())
    ax.set_xlabel('Velocity offset km s$^{-1}$')
    if show:
        pl.show()

    return fig, artists_line, artists_tick


def writesp(filename, sp, resvel=None, overwrite=False):
    """ Writes out a spectrum, as ascii for now - wavelength, flux,
    error, continuum.

    sp must have attributes wa, fl, er and optionally co.

    Keyword overwrite can be True or False.

    resvel means the instrumental fwhm in km/s will be written on the
    first line (VPFIT style).
    """
    if os.path.lexists(filename) and not overwrite:
        c = raw_input('File %s exists - overwrite? (y) or n: ' % filename)
        if c != '':
            if c.strip().lower()[0] == 'n':
                print 'returning without writing anything...'
                return
    fh = open(filename, 'w')
    if resvel is not None:
        fh.write('RESVEL %.2f' % resvel)

    fl = np.nan_to_num(sp.fl)
    er = np.nan_to_num(sp.er)
    if not hasattr(sp, 'co') or np.all(np.isnan(sp.co)):
        for w,f,e in izip(sp.wa, fl, er):
            fh.write("% .12g % #12.8g % #12.8g\n" % (w,f,e))
    else:
        co = np.nan_to_num(sp.co)
        for w,f,e,c in izip(sp.wa, fl, er, co):
            fh.write("% .12g % #12.8g % #12.8g % #12.8g\n" % (w,f,e,c))
    fh.close()


    
def find_cont(fl, fwhm1=300, fwhm2=200, nchunks=4):
    """ Given the flux, estimate the continuum. fwhm values are
    smoothing lengths.
    """
    # smooth flux, with smoothing length much longer than expected
    # emission line widths.
    fl = nan2num(fl.astype(float), replace='mean')
    co = convolve_psf(fl, fwhm1, edge=10)
    
    npts = len(fl)
    indices = np.arange(npts)
    
    # throw away top and bottom 2% of data that deviates from
    # continuum and re-fit new continuum. Go chunk by chunk so that
    # points are thrown away evenly across the spectrum.
    nfl = fl / co
    step = npts // nchunks  + 1
    ind = range(0, npts, step) + [npts]
    igood = []
    for i0,i1 in zip(ind[:-1], ind[1:]):
        isort = nfl[i0:i1].argsort()
        len_isort = len(isort)
        j0,j1 = int(0.05 * len_isort), int(0.95 * len_isort)
        igood.extend(isort[j0:j1]+i0)

    good = np.in1d(indices, igood)
    sfl = fl.copy()
    #pdb.set_trace()
    sfl[~good] = np.interp(indices[~good], indices[good], sfl[good])

    co = convolve_psf(sfl, fwhm2, edge=10)
    return co
    
def find_err(fl, co, nchunks=10):
    """ Given a continuum and flux array, return a very rough estimate
    of the error.
    """
    rms = []
    midpoints = []
    npts = len(fl)
    step = npts // nchunks  + 1
    indices = range(0, npts, step) + [npts]
    for i,j in zip(indices[:-1], indices[1:]):
        #print i,j
        imid = int(0.5 * (i + j))
        midpoints.append(imid)
        df = fl[i:j] - co[i:j]
        n = len(df)
        # throw away top and bottom 5% of these
        df = np.sort(df)[int(0.05*n):int(0.95*n)]
        rms.append(df.std())
    er = np.interp(np.arange(len(fl)), midpoints, rms)
    return er

def pca_qso_cont(nspec, seed=None, return_weights=False):
    """ Make qso continua using the PCA and weights from N. Suzuki et
    al. 2005 and N. Suzuki 2006.

    Input:   nspec (int), number of spectra to create
    Returns: wavelength (shape N), array of spectra [shape (nspec, N)]

    Memory use might be prohibitive for nspec > ~1e4.
    """
    # read the principle components
    path =  os.path.abspath(os.path.dirname(__file__))
    filename = path + '/PCAcont/Suzuki05/tab3.txt'
    names = 'wa,mu,musig,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10'
    co = readtxt(filename, skip=23, names=names)
    # use only the first 7 eigenvetors
    eig = [co['e%i' % i] for i in range(1,8)]
    # from Suzuki et al 2006.
    csig = np.array([7.563, 3.604, 2.351, 2.148, 1.586,
                     1.479, 1.137]) #, 0.778, 0.735, 0.673])

    # generate weights for each eigenvector
    if seed is not None:  np.random.seed(seed)
    weights = []
    for sig in csig:
        temp = np.random.randn(2*nspec)
        # make sure we don't have any very large deviations from the mean
        temp = temp[np.abs(temp) < 3][:nspec]
        assert len(temp) == nspec
        weights.append(temp * sig)

    # generate nspec continua. loop over pixels
    sp = []
    for i in range(len(co.wa)):
        sp.append(co.mu[i] + np.sum(w*e[i] for w,e in zip(weights, eig)))

    sp = np.transpose(sp)
    if return_weights:
        return co.wa, sp, weights
    else:
        return co.wa, sp

def vac2air_ciddor(vacw):
    """ Convert vacuum wavelengths in Angstroms to air wavelengths.

    This uses the relation from Ciddor 1996, Applied Optics LP,
    vol. 35, Issue 9, p.1566. Only valid for wavelengths > 2000 Ang.
    """
    k0 = 238.0185
    k1 = 1e-8 * 5792105.
    k2 = 57.362
    k3 = 1e-8 * 167917.
    s2 = (1e4 / vacw)**2
    n = 1 + k1/(k0 - s2) + k3/(k2 - s2)
    airw = vacw / n

    return airw

def vac2air_morton(vacw):
    """ Convert vacuum wavelengths in Angstroms to air wavelengths.
    
    This uses the relation from Morton 1991, ApJS, 77, 119. Only valid
    for wavelengths > 2000 Ang.  Use this for compatibility with older
    spectra that may have been corrected using the (older) Morton
    relation.  The Ciddor relation used in vac2air_ciddor() is claimed
    to be more accurate at IR wavelengths.
    """
    temp = (1e4 / vacw) ** 2
    airw = 1. / (1. + 6.4328e-5 + 2.94981e-2/(146 - temp) +
                 2.5540e-4/(41 - temp)) * vacw
    return airw

def air2vac_morton(airw):
    """ Convert air wavelengths in Angstroms to vacuum wavelengths.
    
    Uses linear interpolation of the inverse transformation for
    vac2air_morton. The fractional error (wa - watrue) / watrue
    introduced by interpolation is < 5e-8.

    Only valid for wa > 2000 Angstroms.
    """
    nstep = int((airw[-1] - airw[0]) // 50.) + 1
    testvac = np.linspace(airw[0], airw[-1]+10, nstep)
    testair = vac2air_morton(testvac)
    vacw = np.interp(airw, testair, testvac)
    
    return vacw

def air2vac_ciddor(airw):
    """ Convert air wavelengths in Angstroms to vacuum wavelengths.
    
    Uses linear interpolation of the inverse transformation for
    vac2air_ciddor. The fractional error (wa - watrue) / watrue
    introduced by interpolation is < 5e-8.

    Only valid for wa > 2000 Angstroms.
    """
    nstep = int((airw[-1] - airw[0]) // 50.) + 1
    testvac = np.linspace(airw[0], airw[-1]+10, nstep)
    testair = vac2air_ciddor(testvac)
    vacw = np.interp(airw, testair, testvac)

    return vacw

def qso_template(wa, z):
    """ return a QSO spectrum at redshift z.

    The SDSS composite spectrum as a function of F_lambda is returned
    at each wavelength of wa. wa must be in angstroms.
    """
    r = readtabfits(PATH + '/data/templates/qso/dr1QSOspec.fits')
    return np.interp(wa, r.wa*(1+z), r.fl)

def make_constant_dv_wa_scale(wmin, wmax, dv):
    """ Make a constant velocity width scale given a start and end
    wavelength, and velocity pixel width.
    """
    dlogw = np.log10(1 + dv/c_kms)
    # find the number of points needed.
    npts = int(np.log10(wmax / wmin) / dlogw)
    wa = wmin * 10**(np.arange(npts)*dlogw)
    return wa

def convolve_constant_dv(wa, fl, wa_dv=None, npix=4., vfwhm=None):
    """ Convolve a wavelength array with a gaussian of constant
    velocity width.

    If `vfwhm` is specified an intermediate wavelength array with
    constant velocity pixel width is calculated. Otherwise, both
    `wa_dv` and `npix` must be given -- this is faster because no
    intermediate array needs to be calculated.

    Parameters
    ----------
    fl, wa: arrays of floats, length N
        The array to be convolved and its wavelengths.

    vfwhm: float, optional
        Full width at half maximum in velocity space (km/s) of the
        gaussian kernel with which to convolve `fl`.

    wa_dv: array of floats, optional
        Wavelength array with a constant velocity width (this can be
        generated with make_constant_dv_wa_scale()).

    npix: float, optional (default 4)
        Number of pixels in the `wa_dv` array coresponding to the
        gaussian FWHM. This only needs to be large enough to sample
        the gaussian line spread function adequately.

    Returns
    -------
    fl_out: array of length N
        fl convolved with the gaussian kernel with the specified FWHM.
    
    """
    # interpolate to the log-linear scale, convolve, then
    # interpolate back again.
    # convolve with the gaussian
    if vfwhm is not None:
        wa_dv = make_constant_dv_wa_scale(wa[0], wa[-1], float(vfwhm)/npix)
    fl_dv = np.interp(wa_dv, wa, fl)
    fl_dv_smoothed = convolve_psf(fl_dv, npix)
    fl_out = np.interp(wa, wa_dv, fl_dv_smoothed)
    return fl_out
