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
from .protocol_tigre_reconstruction_helps import *



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
    def _defineTigreParams(self, form):

        form.addParam('family',
                      EnumParam,
                      choices=self.FAMILIES,
                      default=self.FAMILY_EXACT,
                      label="Family of methods",
                      display=EnumParam.DISPLAY_COMBO,
                      help=familyHelps['family'])

        form.addParam('exactsMethod',
                      EnumParam,
                      choices=self.ALGORITHMS_EXACT,
                      default=self.WBP,
                      condition='family==%d' % self.FAMILY_EXACT,
                      display=EnumParam.DISPLAY_COMBO,
                      label="Reconstruction Algorithm",
                      help=exactsMethodHelp['exactMethod'])

        form.addParam('gradientMethod',
                      EnumParam,
                      choices=self.ALGORITHMS_GRADIENT,
                      default=self.SART,
                      condition='family==%d' % self.FAMILY_GRADIENT,
                      display=EnumParam.DISPLAY_COMBO,
                      label="Reconstruction Algorithm",
                      help=gradientMethodHelp['gradientMethod'])

        form.addParam('KrylovMethod',
                      EnumParam,
                      choices=self.ALGORITHMS_KRYLOV,
                      default=self.CGLS,
                      condition='family==%d' % self.FAMILY_KRYLOV,
                      display=EnumParam.DISPLAY_COMBO,
                      label="Reconstruction Algorithm",
                      help=KrylovMethodHelp['KrylovMethod'])

        form.addParam('statsMethod',
                      LabelParam,
                      condition='family==%d' % self.FAMILY_STATISTICAL,
                      label="MLEM algorithm",
                      help=statsMethodHelp['statsMethod'])

        form.addParam('varMethod',
                      EnumParam,
                      choices=self.ALGORITHMS_VARIATIONAL,
                      default=self.FISTA,
                      condition='family==%d' % self.FAMILY_VARIATIONAL,
                      display=EnumParam.DISPLAY_COMBO,
                      label="Reconstruction Algorithm",
                      help=varMethodHelp['varMethod'])

        filterCondition = 'family==%d' % self.FAMILY_EXACT

        form.addParam('filter',
                      EnumParam,
                      choices=self.FILTER_LIST,
                      default=self.RAMLAK,
                      condition=filterCondition,
                      display=EnumParam.DISPLAY_COMBO,
                      label="Filter",
                      help=filterHelps['filter'])




        # gradien methods form
        self._defineSARTParams(form)
        self._defineSIRTParams(form)
        self._defineOSSARTParams(form)
        self._defineASDPOCSParams(form)
        self._defineOSASDPOCSParams(form)
        self._definePCSDParams(form)
        self._defineAWPCSDParams(form)
        self._defineAWASDPOCSParams(form)

        # Krylov form
        self._defineCGLSParams(form)
        self._defineLSQRParams(form)
        self._defineLSMRParams(form)
        self._defineHYBRIDLSQRParams(form)
        self._defineABGMRESParams(form)
        self._defineBAGMRESParams(form)
        self._defineIRNTVCGLSParams(form)

        # statistical methods form
        self._defineMLEMParams(form)

        # Variational methods form
        self._defineFISTAParams(form)
        self._defineSARTTVParams(form)


    def _defineSARTParams(self, form):
        form.addParam('SARTiter',
                      IntParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.SART),
                      default=30,
                      label='Iterations',
                      help=SARThelps['SARTiter'])

        form.addParam('SARTlambda',
                      FloatParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.SART),
                      default=1.0,
                      label='λ regularizer',
                      help=SARThelps['SARTlambda'])

        form.addParam('SARTlambdared',
                      FloatParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.SART),
                      default=0.999,
                      label='λ red',
                      help=SARThelps['SARTlambdared'])

    def _defineSIRTParams(self, form):
        form.addParam('SIRTiter',
                      IntParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.SIRT),
                      default=30,
                      label='Iterations',
                      help=SIRThelps['SIRTiter'])

        form.addParam('SIRTlambda',
                      FloatParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.SIRT),
                      default=1.0,
                      label='lambda regularizer',
                      help=SIRThelps['SIRTlambda'])

        form.addParam('SIRTlambdared',
                      FloatParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.SIRT),
                      default=0.999,
                      label='lambda red',
                      help=SIRThelps['SIRTlambdared'])

    def _defineOSSARTParams(self, form):
        form.addParam('OSSARTiter',
                      IntParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.OSSART),
                      default=30,
                      label='Iterations',
                      help=OSSARThelps['OSSARTiter'])

        form.addParam('OSSARTlambda',
                      FloatParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.OSSART),
                      default=1.0,
                      label='lambda regularizer',
                      help=OSSARThelps['OSSARTlambda'])

        form.addParam('OSSARTlambdared',
                      FloatParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.OSSART),
                      default=0.999,
                      label='lambda red',
                      help=OSSARThelps['OSSARTlambdared'])

    def _defineASDPOCSParams(self, form):
        form.addParam('ASDPOCSiter',
                      IntParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.ASDPOCS),
                      default=30,
                      label='Iterations',
                      help=ASDPOCShelps['ASDPOCSiter'])

        form.addParam('ASDPOCSlambda',
                      FloatParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.ASDPOCS),
                      default=1.0,
                      label='lambda regularizer',
                      help=ASDPOCShelps['ASDPOCSlambda'])

        form.addParam('ASDPOCSalpha',
                      FloatParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.ASDPOCS),
                      default=0.002,
                      label='alpha',
                      help=ASDPOCShelps['ASDPOCSalpha'])

        form.addParam('ASDPOCStviter',
                      IntParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.ASDPOCS),
                      default=25,
                      label='TV iterations',
                      help=ASDPOCShelps['ASDPOCStviter'])

        form.addParam('ASDPOCSlambdared',
                      FloatParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.ASDPOCS),
                      default=0.9999,
                      label='lambda red',
                      help=ASDPOCShelps['ASDPOCSlambdared'])

        form.addParam('ASDPOCSalphared',
                      FloatParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.ASDPOCS),
                      default=0.95,
                      label='alpha red',
                      help=ASDPOCShelps['ASDPOCSalphared'])

        form.addParam('ASDPOCSratio',
                      FloatParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.ASDPOCS),
                      default=0.94,
                      label='Ratio',
                      help=ASDPOCShelps['ASDPOCSratio'])

    def _defineOSASDPOCSParams(self, form):
        form.addParam('OSASDPOCSiter',
                      IntParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.OSASDPOCS),
                      default=20,
                      label='Iterations',
                      help=OSASDPOCShelps['OSASDPOCSiter'])

        form.addParam('OSASDPOCSlambda',
                      FloatParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.OSASDPOCS),
                      default=1.0,
                      label='lambda regularizer',
                      help=OSASDPOCShelps['OSASDPOCSlambda'])

        form.addParam('OSASDPOCSalpha',
                      FloatParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.OSASDPOCS),
                      default=0.002,
                      label='alpha',
                      help=OSASDPOCShelps['OSASDPOCSalpha'])

        form.addParam('OSASDPOCStviter',
                      IntParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.OSASDPOCS),
                      default=25,
                      label='TV iterations',
                      help=OSASDPOCShelps['OSASDPOCStviter'])

        form.addParam('OSASDPOCSlambdared',
                      FloatParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.OSASDPOCS),
                      default=0.9999,
                      label='lambda red',
                      help=OSASDPOCShelps['OSASDPOCSlambdared'])

        form.addParam('OSASDPOCSalphared',
                      FloatParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.OSASDPOCS),
                      default=0.95,
                      label='alpha red',
                      help=OSASDPOCShelps['OSASDPOCSalphared'])

        form.addParam('OSASDPOCSratio',
                      FloatParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.OSASDPOCS),
                      default=0.94,
                      label='Ratio',
                      help=OSASDPOCShelps['OSASDPOCSratio'])


    def _definePCSDParams(self, form):
        form.addParam('PCSDiter',
                      IntParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.PCSD),
                      default=20,
                      label='Iterations',
                      help=PCSDhelps['PCSDiter'])

    def _defineAWPCSDParams(self, form):
        form.addParam('AWPCSDiter',
                      IntParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.AWPCSD),
                      default=30,
                      label='Iterations',
                      help=AWPCSDhelps['AWPCSDiter'])

    def _defineAWASDPOCSParams(self, form):
        form.addParam('AWASDPOCSiter',
                      IntParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.AWASDPOCS),
                      default=20,
                      label='Iterations',
                      help=AWASDPOCShelps['AWASDPOCSiter'])

        form.addParam('AWASDPOCSlambda',
                      FloatParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.AWASDPOCS),
                      default=1.0,
                      label='lambda regularizer',
                      help=AWASDPOCShelps['AWASDPOCSlambda'])             

        form.addParam('AWASDPOCSalpha',
                      FloatParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.AWASDPOCS),
                      default=0.002,
                      label='alpha',
                      help=AWASDPOCShelps['AWASDPOCSalpha'])

        form.addParam('AWASDPOCStviter',
                      IntParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.AWASDPOCS),
                      default=25,
                      label='TV iterations',
                      help=AWASDPOCShelps['AWASDPOCStviter'])

        form.addParam('AWASDPOCSlambdared',
                      FloatParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.AWASDPOCS),
                      default=0.9999,
                      label='lambda red',
                      help=AWASDPOCShelps['AWASDPOCSlambdared'])

        form.addParam('AWASDPOCSalphared',
                      FloatParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.AWASDPOCS),
                      default=0.95,
                      label='alpha red',
                      help=AWASDPOCShelps['AWASDPOCSalphared'])

        form.addParam('AWASDPOCSratio',
                      FloatParam,
                      condition='family==%d and gradientMethod==%d' % (self.FAMILY_GRADIENT, self.AWASDPOCS),
                      default=0.94,
                      label='Ratio',
                      help=AWASDPOCShelps['AWASDPOCSratio'])

    def _defineCGLSParams(self, form):
        form.addParam('CGLSiter',
                      IntParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_KRYLOV, self.CGLS),
                      default=30,
                      label='Iterations',
                      help=CGLShelps['CGLSiter'])

    def _defineLSQRParams(self, form):
        form.addParam('LSQRiter',
                      IntParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_KRYLOV, self.LSQR),
                      default=30,
                      label='Iterations',
                      help=LSQRhelps['LSQRiter'])

    def _defineLSMRParams(self, form):
        form.addParam('LSMRiter',
                      IntParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_KRYLOV, self.LSMR),
                      default=30,
                      label='Iterations',
                      help=LSMRhelps['LSMRiter'])
        form.addParam('LSMRlambda',
                      FloatParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_KRYLOV, self.LSMR),
                      default=1.0,
                      label='lambda',
                      help=LSMRhelps['LSMRlambda'])

    def _defineHYBRIDLSQRParams(self, form):
        form.addParam('hybridLSQRiter',
                      IntParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_KRYLOV, self.HYBRIDLSQR),
                      default=30,
                      label='Iterations',
                      help=hybridLSQRhelps['hybridLSQRiter'])

    def _defineABGMRESParams(self, form):
        form.addParam('ABGMRESiter',
                      IntParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_KRYLOV, self.ABGMRES),
                      default=30,
                      label='Iterations',
                      help=ABGMREShelps['ABGMRESiter'])

    def _defineBAGMRESParams(self, form):
        form.addParam('BAGMRESiter',
                      IntParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_KRYLOV, self.BAGMRES),
                      default=30,
                      label='Iterations',
                      help=BAGMREShelps['BAGMRESiter'])


    def _defineIRNTVCGLSParams(self, form):
        form.addParam('IRNTVCGLSiter',
                      IntParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_KRYLOV, self.IRNTVCGLS),
                      default=10,
                      label='Iterations',
                      help=IRNTVCGLShelps['IRNTVCGLSiter'])

        form.addParam('IRNTVCGLStvlambda',
                      FloatParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_KRYLOV, self.IRNTVCGLS),
                      default=5.0,
                      label='TV regularizer',
                      help=IRNTVCGLShelps['IRNTVCGLStvlambda'])
        form.addParam('IRNTVCGLSinnerouter',
                      IntParam,
                      condition='family==%d and KrylovMethod==%d' % (self.FAMILY_KRYLOV, self.IRNTVCGLS),
                      default=2,
                      label='number of estimations for weights',
                      help=IRNTVCGLShelps['IRNTVCGLSinnerouter'])

    def _defineMLEMParams(self, form):
        form.addParam('MLEMiter',
                      IntParam,
                      condition='family==%d' % self.FAMILY_STATISTICAL,
                      default=500,
                      label='Iterations',
                      help=MLEMhelps['MLEMiter'])

    def _defineFISTAParams(self, form):
        form.addParam('FISTAiter',
                      IntParam,
                      condition='family==%d and varMethod==%d' % (self.FAMILY_VARIATIONAL, self.FISTA),
                      default=50,
                      label='Iterations',
                      help=FISTAhelps['FISTAiter'])

        form.addParam('FISTAtvlambda',
                      FloatParam,
                      condition='family==%d and varMethod==%d' % (self.FAMILY_VARIATIONAL, self.FISTA),
                      default=10.0,
                      label='TV Regularizer',
                      help=FISTAhelps['FISTAtvlambda'])
        form.addParam('FISTAtviter',
                      IntParam,
                      condition='family==%d and varMethod==%d' % (self.FAMILY_VARIATIONAL, self.FISTA),
                      default=50,
                      label='TV iterations',
                      help=FISTAhelps['FISTAtviter'])

    def _defineSARTTVParams(self, form):
        form.addParam('SARTTViter',
                      IntParam,
                      condition='family==%d and varMethod==%d' % (self.FAMILY_VARIATIONAL, self.SARTTV),
                      default=30,
                      label='Iterations',
                      help=SARTTVhelps['SARTTViter'])
        form.addParam('SARTTVtvlambda',
                      FloatParam,
                      condition='family==%d and varMethod==%d' % (self.FAMILY_VARIATIONAL, self.SARTTV),
                      default=50,
                      label='TV Regularizer',
                      help=SARTTVhelps['SARTTVtvlambda'])
        form.addParam('SARTTVtviter',
                      IntParam,
                      condition='family==%d and varMethod==%d' % (self.FAMILY_VARIATIONAL, self.SARTTV),
                      default=50,
                      label='TV iterations',
                      help=SARTTVhelps['SARTTVtviter'])
