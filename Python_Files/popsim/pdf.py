# -*- coding: utf-8 -*-
"""
pdf.py: module for working with probability density functions
"""
# Copyright (C) 2011-2012 Steffen Waldherr waldherr@ist.uni-stuttgart.de
# Time-stamp: <Last change 2013-11-25 11:08:20 by Steffen Waldherr>

import scipy as sp
import numpy as np


from scipy import stats
from math import exp

class ScalarPDF(object):
    """
    probability density function in one variable
    """
    def __init__(self, varax, pdf, normalize=False):
        """
        initialize new ScalarPDF over vector varax with
        probability density values pdf.
        """
        if np.min(np.diff(varax)) <= 0:
            raise ValueError("Entries of varax must be continuously increasing.")
        if np.min(pdf) < 0:
            raise ValueError("Entries of pdf must be non-negative.")
        self.varax = np.array(varax, dtype=np.float64, copy=True)
        self.pdf = np.array(pdf, dtype=np.float64, copy=True)
        if normalize:
            self.normalize()
        else:
            self.cdf = np.copy(self.pdf)
            self.cdf[0] = 0
            for i in xrange(1,len(self.cdf)):
                self.cdf[i] = self.cdf[i-1]+0.5*(self.pdf[i]+self.pdf[i-1])*(self.varax[i]-self.varax[i-1])

    def __mul__(self, scalar):
        """
        return scalar * PDF (no renormalization)
        """
        res = ScalarPDF(self.varax, self.pdf)
        res.pdf *= scalar
        res.cdf *= scalar
        return res

    def __rmul__(self, scalar):
        return self.__mul__(scalar)

    def __add__(self, other):
        """
        return PDF + otherPDF (doing renormalization)
        new varax is computed as "average" of old varaxes
        """
        if other == 0:
            return self
        vmin = min(self.varax[0], other.varax[0])
        vmax = max(self.varax[-1], other.varax[-1])
        num = 0.5*(len(self.varax) + len(other.varax))
        newvarax = np.linspace(vmin, vmax, num)
        pdf = self.evalpdf(newvarax) + other.evalpdf(newvarax)
        return ScalarPDF(newvarax, pdf)

    def __radd__(self, other):
        return self.__add__(other)

    def evalpdf(self, values):
        """
        get PDF at values
        """
        vals = values
        return np.interp(vals, self.varax, self.pdf, left=0, right=0)

    def sample(self,num=1):
        """
        generate a vector of num samples from this pdf.
        use num=0 to sample a scalar
        """
        s = np.random.random((max(num,1),))
        sample = np.interp(s, self.cdf, self.varax)
        if num == 0:
            return sample[0]
        else:
            return sample

    def normalize(self):
        """
        normalize this PDF to integral over density = 1
        """
        self.cdf = np.copy(self.pdf)
        self.cdf[0] = 0
        for i in xrange(1,len(self.cdf)):
            self.cdf[i] = self.cdf[i-1]+0.5*(self.pdf[i]+self.pdf[i-1])*(self.varax[i]-self.varax[i-1])
        scale = self.cdf[-1]
        self.cdf /= scale
        self.pdf /= scale
        return self

    def mean(self):
        """
        compute mean value for this PDF
        """
        return sp.integrate.trapz(self.varax*self.pdf, x=self.varax)

    def std(self):
        """
        compute standard deviation for this PDF
        """
        return np.sqrt(sp.integrate.trapz(self.varax**2 * self.pdf, x=self.varax) - self.mean()**2)

    def median(self):
        """
        compute median value for this PDF
        """
        return np.interp(0.5, self.cdf, self.varax)

    def shift_to_median(self, newmedian):
        """
        shift variables axis such that median becomes newmedian
        """
        oldmed = self.median()
        self.varax += newmedian - oldmed
        self.normalize()

    def scale_to_median(self, newmedian):
        """
        scale variables axis such that median becomes newmedian
        """
        oldmed = self.median()
        self.varax *= newmedian/oldmed
        self.normalize()

    def scale_to_max(self, newmax):
        """
        scale variables axis such that maximum becomes newmax
        """
        self.varax *= newmax/self.varax[-1]
        self.normalize()

class LogScaleScalarPDF(object):
    """
    probability density function in one variable
    """
    def __init__(self, varax, pdf, logscale=False):
        """
        initialize new ScalarPDF over vector varax with
        probability density values pdf.
        """
        if np.min(np.diff(varax)) <= 0:
            raise ValueError("Entries of varax must be continuously increasing.")
        if np.min(pdf) < 0:
            raise ValueError("Entries of pdf must be non-negative.")
        self.varax = np.array(varax, dtype=np.float64, copy=True)
        self.pdf = np.array(pdf, dtype=np.float64, copy=True)
#         self.pdf /= np.sum(self.pdf)
        self.logscale = logscale
        self.cdf = np.copy(self.pdf)
        self.cdf[0] = 0
        for i in xrange(1,len(self.cdf)):
            self.cdf[i] = self.cdf[i-1]+0.5*(self.pdf[i]+self.pdf[i-1])*(self.varax[i]-self.varax[i-1])
        self.cdf /= self.cdf[-1]

    # def evalpdf(self, values, logscale=None):
    #     """
    #     get PDF at values
    #     if logscale is True, values are given in log10 base.
    #     """
    #     if logscale is None:
    #         logscale = self.logscale
    #     vals = values
    #     if logscale and not self.logscale:
    #         vals = 10**values
    #     if not logscale and self.logscale:
    #         vals = np.log10(values)
    #     pdf = np.interp(vals, self.varax, self.pdf, left=0, right=0)
    #     if logscale:

    def sample(self,num=1,logscale=None):
        """
        generate a vector of num samples from this pdf.
        use num=0 to sample a scalar
        """
        if logscale is None:
            logscale = self.logscale
        s = np.random.random((max(num,1),))
        sample = np.interp(s, self.cdf, self.varax)
        if logscale and not self.logscale:
            sample = np.log10(sample)
        if not logscale and self.logscale:
            sample = 10**sample
        if num == 0:
            return sample[0]
        else:
            return sample 

    def logtransform(self):
        """
        return new pdf, transformed to or from log scale
        """
        if self.logscale:
            ax = 10**self.varax
            pdf = self.pdf/ax/np.log(10)
            return LogScaleScalarPDF(ax, pdf, logscale=False)
        else:
            if np.min(self.varax) <= 0:
                raise ValueError("Log transform is impossible for non-positive variables.")
            logax = np.log10(self.varax)
            logpdf = np.log(10)*self.varax*self.pdf
            return LogScaleScalarPDF(logax, logpdf, logscale=True)

def norm2log(pdf):
    """
    transform to logscale
    """
    if np.min(pdf.varax) <= 0:
        raise ValueError("Log transform is impossible for non-positive variables.")
    logax = np.log10(pdf.varax)
    logpdf = np.log(10)*pdf.varax*pdf.pdf
    return ScalarPDF(logax, logpdf)

def log2norm(pdf):
    """
    transform from logscale to standard coordinates
    """
    ax = 10**pdf.varax
    pdftab = pdf.pdf/ax/np.log(10)
    return ScalarPDF(ax, pdftab)

def estimate(scalardata, method="naive", varax=None, h=None, points=200):
    """
    estimate a ScalarPDF from 1d-array scalardata using a density estimator.

    Choices for method:
    - naive: naive density estimator
    - hist: histogram estimator
    - kde: kernel density estimator from scipy.stats.gaussian_kde

    varax is a 1d-array to define the variable axis, if it is None,
      use numpy.linspace(min(scalardata)-h/2, max(scalardata)+h/2, num=200)
    h is the bandwidth, if it is None, use 0.05*(max(scalardata)-min(scalardata))
    """
    scalardata = np.asarray(scalardata,dtype=np.float64).flatten()
    if h is None:
        h = 0.05*(np.max(scalardata)-np.min(scalardata))
    if varax is None:
        varax = np.linspace(np.min(scalardata)-0.5*h, np.max(scalardata)+0.5*h, num=points, endpoint=True)
    if method == "naive":
        return naive_estimator(scalardata, varax, h)
    elif method == "hist":
        binnum = (np.max(scalardata)-np.min(scalardata))/h
        histrange = (np.min(varax), np.max(varax))
        hist, binedges = np.histogram(scalardata, normed=True, bins=binnum, range=histrange, new=True)
        bincenters = [0.5*(b[0]+b[1]) for b in zip(binedges[:-1], binedges[1:])]
        pdf = np.interp(varax, bincenters, hist, left=0, right=0)
        return ScalarPDF(varax, pdf)
    elif method == "kde":
        kde = stats.gaussian_kde(scalardata)
        return ScalarPDF(varax, kde(varax))

def naive_estimator(scalardata, varax=None, h=None):
    """
    get a ScalarPDF from 1d-array scalardata using a naive estimator.
    varax is a 1d-array to define the variable axis, if it is None,
      use numpy.linspace(min(scalardata)-h/2, max(scalardata)+h/2, num=200)
    h is the bandwidth, if it is None, use 0.05*(max(scalardata)-min(scalardata))
    """
    scalardata = np.asarray(scalardata,dtype=np.float64).flatten()
    if h is None:
        h = 0.05*(np.max(scalardata)-np.min(scalardata))
    if varax is None:
        varax = np.linspace(np.min(scalardata)-0.5*h, np.max(scalardata)+0.5*h, num=200, endpoint=True)
    pdf = 0*varax
    n = scalardata.shape[0]
    for i,d in enumerate(scalardata):
        pdf[np.abs(d-varax) <= h] += 0.5
    pdf /= n*h
    return ScalarPDF(varax, pdf)

def gamma(theta, k, numpoints=200, rangetheta=10.0):
    """
    compute a gamma distribution with parameters theta and k
    options:
    numpoints - number of points to compute pdf at
    rangestd - factor of theta up to which varax should range
    """
    varax = np.linspace(0, rangetheta*theta, num=numpoints)
    pdf = cmix.wrap_gsl_ran_gamma_pdf(theta, k, varax)
    return ScalarPDF(varax, pdf)

def lognormal(center, std, numpoints=200, rangestd=4.0):
    """
    compute a log10normal distribution centered at center with log10 standard deviation std
    options:
    numpoints - number of points to compute pdf at
    rangestd - factor of standard deviation over which varax should range
    """
    logcenter = np.log10(center)
    varax = np.linspace(logcenter-rangestd*std, logcenter+rangestd*std, num=numpoints)
    pdf = sp.stats.norm.pdf(varax, logcenter, std)
    return log2norm(ScalarPDF(varax, pdf))

def normal(center, std, numpoints=200, rangestd=4.0, varax=None):
    """
    compute a normal distribution centered at center with standard deviation std
    options:
    numpoints - number of points to compute pdf at
    rangestd - factor of standard deviation over which varax should range
    """
    if varax is None:
        varax = np.linspace(center-rangestd*std, center+rangestd*std, num=numpoints)
    pdf = sp.stats.norm.pdf(varax, center, std)
    return ScalarPDF(varax, pdf)

def fromprobfun(varax, pdf):
    """
    construct ScalarPDF from callable pdf along varax
    """
    pdfvals = np.asarray([pdf(xi) for xi in varax]).flatten()
    return ScalarPDF(varax, pdfvals)

def equal(minimum=0.,maximum=1.,numpoints=200):
    """
    compute equal distribution
    numpoints - number of points to compute pdf at
    """
    varax=np.linspace(minimum, maximum, numpoints)
    pdf=np.ones(numpoints)/(maximum-minimum)
    return ScalarPDF(varax,pdf)



def trapez(minimum=0.,maximum=1., left=0., right=0.,numpoints=200):
    """
    compute trapezoidal distribution
    numpoints - number of points to compute pdf at
    normalizes the distribution
    """
    varax=np.linspace(minimum, maximum, numpoints)
    pdf=np.linspace(left,right,numpoints)/(0.5*(left+right)*(maximum-minimum))
    return ScalarPDF(varax,pdf)


def trapeznew(minimum=0.,maximum=1., left=0., right=0.,numpoints=200):
    """
    compute trapezoidal distribution
    numpoints - number of points to compute pdf at
    normalizes the distribution
    """
    varax=np.linspace(minimum, maximum, numpoints)
    pdf=np.linspace(left,right,numpoints)/(0.5*(left+right)*(maximum-minimum))
    return ScalarPDF(varax,pdf)


'''
uniform distibution: numpy.random.uniform(a = low, b = high) --> including low, excluding high
    --> its pdf is p(x) = 1/(b - a) --> all floats
numpy.randint bzw. numpy.random_integers --> all integers (discrete uniform distribution)
    --> pdf: equal(minimum, maximum, numpoints)
    
p(a) (a = age): p(a) = 2 * gamma * exp(-gamma * a), mit gamma = log(2) / T, T = 25.1h
    --> cumulative density function: cdf is integralover 0-T of pdf 
        cdf = sp.integrate.quad(pdf, 0, T)
'''
# x = sample from uniform distribution
# pdf = P(a)
# x = P(a) --> P-1(x) = P-1(P(a)) = a
def my_pdf(minimum=0., maximum=1., gamma=1., numpoints=300):
    varax = np.linspace(minimum, maximum, numpoints)
    pdf = 2 * gamma * np.exp(-gamma * varax)
    return ScalarPDF(varax, pdf)
    
    