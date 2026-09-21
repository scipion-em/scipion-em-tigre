# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Jose Luis Vilas (jlvilas@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
from enum import Enum
import pyworkflow.utils as pwutils
from pyworkflow import VERSION_3_0
from pyworkflow.object import Set
from pyworkflow.protocol.params import (PointerParam, IntParam, FloatParam, BooleanParam, LabelParam, EnumParam, StringParam, LEVEL_ADVANCED, GPU_LIST)
import pyworkflow.utils.path as path

from pwem.protocols import EMProtocol
import pwem.emlib.metadata as md

from tomo.protocols.protocol_base import ProtTomoBase
from tomo.objects import Tomogram, SetOfTomograms, SetOfTiltSeries
from pyworkflow import BETA
from tigre import Plugin


class ProtTigreForm(ProtTomoBase):

    # Definition of the families
    FAMILIES = ['Exact algorithms', 'Gradient-based algorithms', 'Krylov algorithms', 'Statistical reconstructions', 'Variational methods']
    FAMILY_EXACT       = 0
    FAMILY_GRADIENT    = 1
    FAMILY_KRYLOV      = 2
    FAMILY_STATISTICAL = 3
    FAMILY_VARIATIONAL = 4

    # Exact algorithms
    ALGORITHMS_EXACT = ['WBP', 'FDK']
    WBP  = 0
    FDK  = 1

    #Gradient based algorithms
    ALGORITHMS_GRADIENT = ['SART', 'SIRT', 'OS-SART', 'ASD-POCS', 'OS-ASD-POCS', 'PCSD', 'AW-PCSD', 'AW-ASD-POCS']
    SART       = 0
    SIRT       = 1
    OSSART     = 2
    ASDPOCS    = 3
    OSASDPOCS  = 4
    PCSD       = 5
    AWPCSD     = 6
    AWASDPOCS  = 7

    #Krylov subspace algorithms
    ALGORITHMS_KRYLOV = ['CGLS', 'LSQR', 'LSMR', 'hybrid-LSQR', 'AB-GMRES', 'BA-GMRES', 'IRN-TV-CGLS']
    CGLS           = 0
    LSQR           = 1
    LSMR           = 2
    HYBRIDLSQR     = 3
    ABGMRES        = 4
    BAGMRES        = 5
    IRNTVCGLS      = 6

    # Statistical algorithms
    ALGORITHMS_STATISTICAL = ['MLEM']
    MLEM = 0

    # Variational algorithms
    ALGORITHMS_VARIATIONAL = ['FISTA', 'SART-TV']
    FISTA  = 0
    SARTTV = 1

    exactList     = []
    grandientList = []
    krylovList    = []
    statsList     = []
    varList       = []

    #Filters list
    FILTER_LIST = ['Ram-lak', 'Shepp-Logan', 'cosine', 'Hamming', 'Hann']
    RAMLAK      = 0
    SHEPP_LOGAN = 1
    COSINE      = 2
    HAMMING     = 3
    HANN        = 4

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.Tomograms = None

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):

        self._defineWBPParams(form)
        self._defineFBPParams(form)

        self._defineSARTParams(form)
        self._defineSIRTParams(form)
        self._defineOSSARTParams(form)
        self._defineASDPOCSParams(form)
        self._defineOSASDPOCSParams(form)
        self._definePCSDParams(form)
        self._defineAWPCSDParams(form)
        self._defineAWASDPOCSParams(form)

        self._defineCGLSParams(form)
        self._defineLSQRParams(form)
        self._defineLSMRParams(form)
        self._defineHYBRIDLSQRParams(form)
        self._defineABGMRESParams(form)
        self._defineBAGMRESParams(form)
        self._defineIRNTVCGLSParams(form)
        self._defineMLEMParams(form)
        self._defineFISTAParams(form)
        self._defineSARTTVParams(form)


    def _defineWBPParams(self, form):
        pass

    def _defineFBPParams(self, form):
        pass

    def _defineSARTParams(self, form):
        form.addParam('SARTiter',
                      IntParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.SART),
                      default=30,
                      label='Iterations',
                      help='N')

        form.addParam('SARTlambda',
                      FloatParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.SART),
                      default=30,
                      label='Iterations',
                      help='N')

        form.addParam('SARTlambdared',
                      FloatParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.SART),
                      default=0.9999,
                      label='lambda red',
                      help='N')

    def _defineSIRTParams(self, form):
        form.addParam('SIRTiter',
                      IntParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.SIRT),
                      default=30,
                      label='Iterations',
                      help='N')

        form.addParam('SIRTlambda',
                      FloatParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.SIRT),
                      default=30,
                      label='Iterations',
                      help='N')

        form.addParam('SIRTlambdared',
                      FloatParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.SIRT),
                      default=0.9999,
                      label='lambda red',
                      help='N')

    def _defineOSSARTParams(self, form):
        form.addParam('OSSARTiter',
                      IntParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.OSSART),
                      default=30,
                      label='Iterations',
                      help='N')

        form.addParam('OSSARTlambda',
                      FloatParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.OSSART),
                      default=30,
                      label='Iterations',
                      help='N')

        form.addParam('OSSARTlambdared',
                      FloatParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.OSSART),
                      default=0.9999,
                      label='lambda red',
                      help='N')

    def _defineASDPOCSParams(self, form):
        form.addParam('ASDPOCSiter',
                      IntParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.ASDPOCS),
                      default=30,
                      label='Iterations',
                      help='N')

        form.addParam('ASDPOCSlambda',
                      FloatParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.ASDPOCS),
                      default=30,
                      label='Iterations',
                      help='N')

        form.addParam('ASDPOCSalpha',
                      FloatParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.ASDPOCS),
                      default=0.002,
                      label='alpha',
                      help='N')

        form.addParam('ASDPOCStviter',
                      IntParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.ASDPOCS),
                      default=25,
                      label='TV iterations',
                      help='N')

        form.addParam('ASDPOCSlambdared',
                      FloatParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.ASDPOCS),
                      default=0.9999,
                      label='lambda red',
                      help='N')

        form.addParam('ASDPOCSalphared',
                      FloatParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.ASDPOCS),
                      default=0.95,
                      label='Alpha red',
                      help='N')

        form.addParam('ASDPOCSratio',
                      FloatParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.ASDPOCS),
                      default=0.94,
                      label='Ratio',
                      help='N')

    def _defineOSASDPOCSParams(self, form):
        form.addParam('OSASDPOCSiter',
                      IntParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.OSASDPOCS),
                      default=30,
                      label='Iterations',
                      help='N')

        form.addParam('OSASDPOCSlambda',
                      FloatParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.OSASDPOCS),
                      default=30,
                      label='Iterations',
                      help='N')

        form.addParam('OSASDPOCSalpha',
                      FloatParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.OSASDPOCS),
                      default=0.002,
                      label='alpha',
                      help='N')

        form.addParam('OSASDPOCStviter',
                      IntParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.OSASDPOCS),
                      default=25,
                      label='TV iterations',
                      help='N')

        form.addParam('OSASDPOCSlambdared',
                      FloatParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.OSASDPOCS),
                      default=0.9999,
                      label='lambda red',
                      help='N')

        form.addParam('OSASDPOCSalphared',
                      FloatParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.OSASDPOCS),
                      default=0.95,
                      label='Alpha red',
                      help='N')

        form.addParam('OSASDPOCSratio',
                      FloatParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.OSASDPOCS),
                      default=0.94,
                      label='Ratio',
                      help='N')


    def _definePCSDParams(self, form):
        form.addParam('PCSDiter',
                      IntParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.PCSD),
                      default=30,
                      label='Iterations',
                      help='N')

    def _defineAWPCSDParams(self, form):
        form.addParam('AWPCSDiter',
                      IntParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.AWPCSD),
                      default=30,
                      label='Iterations',
                      help='N')

    def _defineAWASDPOCSParams(self, form):
        form.addParam('AWASDPOCSiter',
                      IntParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.AWASDPOCS),
                      default=30,
                      label='Iterations',
                      help='N')

        form.addParam('AWASDPOCSalpha',
                      FloatParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.AWASDPOCS),
                      default=0.002,
                      label='alpha',
                      help='N')

        form.addParam('AWASDPOCStviter',
                      IntParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.AWASDPOCS),
                      default=25,
                      label='TV iterations',
                      help='N')

        form.addParam('AWASDPOCSlambdared',
                      FloatParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.AWASDPOCS),
                      default=0.9999,
                      label='lambda red',
                      help='N')

        form.addParam('AWASDPOCSalphared',
                      FloatParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.AWASDPOCS),
                      default=0.95,
                      label='Alpha red',
                      help='N')

        form.addParam('AWASDPOCSratio',
                      FloatParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_GRADIENT, self.AWASDPOCS),
                      default=0.94,
                      label='Ratio',
                      help='N')

    def _defineCGLSParams(self, form):
        form.addParam('CGLSiter',
                      IntParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_KRYLOV, self.CGLS),
                      default=30,
                      label='Iterations',
                      help='N')

    def _defineLSQRParams(self, form):
        form.addParam('LSQRiter',
                      IntParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_KRYLOV, self.LSQR),
                      default=30,
                      label='Iterations',
                      help='N')

    def _defineLSMRParams(self, form):
        form.addParam('LSMRiter',
                      IntParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_KRYLOV, self.LSMR),
                      default=30,
                      label='Iterations',
                      help='N')

    def _defineHYBRIDLSQRParams(self, form):
        form.addParam('HYBRIDLSQRiter',
                      IntParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_KRYLOV, self.HYBRIDLSQR),
                      default=30,
                      label='Iterations',
                      help='N')

    def _defineABGMRESParams(self, form):
        form.addParam('ABGMRESiter',
                      IntParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_KRYLOV, self.ABGMRES),
                      default=30,
                      label='Iterations',
                      help='N')

    def _defineBAGMRESParams(self, form):
        form.addParam('BAGMRESiter',
                      IntParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_KRYLOV, self.BAGMRES),
                      default=30,
                      label='Iterations',
                      help='N')


    def _defineIRNTVCGLSParams(self, form):
        form.addParam('IRNTVCGLSiter',
                      IntParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_KRYLOV, self.IRNTVCGLS),
                      default=10,
                      label='Iterations',
                      help='N')

        form.addParam('IRNTVCGLStvlambda',
                      FloatParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_KRYLOV, self.IRNTVCGLS),
                      default=5,
                      label='TV Regularizer',
                      help='Num')

    def _defineMLEMParams(self, form):
        form.addParam('MLEMiter',
                      IntParam,
                      condition='family==%d and varMethod==%d' % (self.FAMILY_STATISTICAL, self.MLEM),
                      default=500,
                      label='Iterations',
                      help='N')

    def _defineFISTAParams(self, form):
        form.addParam('FISTAiter',
                      IntParam,
                      condition='family==%d and varMethod==%d' % (self.FAMILY_VARIATIONAL, self.FISTA),
                      default=50,
                      label='Iterations',
                      help='N')

        form.addParam('FISTAtvlambda',
                      FloatParam,
                      condition='family==%d and varMethod==%d' % (self.FAMILY_VARIATIONAL, self.FISTA),
                      default=10,
                      label='TV Regularizer',
                      help='Num')
        form.addParam('FISTAtviter',
                      IntParam,
                      condition='family==%d and varMethod==%d' % (self.FAMILY_VARIATIONAL, self.FISTA),
                      default=50,
                      label='TV iterations',
                      help='Num')

    def _defineSARTTVParams(self, form):
        form.addParam('SARTTViter',
                      IntParam,
                      condition='family==%d and varMethod==%d' % (self.FAMILY_VARIATIONAL, self.SARTTV),
                      default=30,
                      label='Iterations',
                      help='N')
        form.addParam('SARTTVtvlambda',
                      FloatParam,
                      condition='family==%d and varMethod==%d' % (self.FAMILY_VARIATIONAL, self.SARTTV),
                      default=50,
                      label='TV Regularizer',
                      help='Num')
        form.addParam('SARTTVtviter',
                      IntParam,
                      condition='family==%d and varMethod==%d' % (self.FAMILY_VARIATIONAL, self.SARTTV),
                      default=50,
                      label='TV iterations',
                      help='Num')
        form.addParam('SARTTValphared',
                      FloatParam,
                      expertLevel=LEVEL_ADVANCED,
                      condition='family==%d and varMethod==%d' % (self.FAMILY_VARIATIONAL, self.SARTTV),
                      default=0.95,
                      label='alpha red',
                      help='Num')
