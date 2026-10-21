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
from .protocol_tigre_reconstruction_base import ProtTigreForm

from pwem.protocols import EMProtocol
import pwem.emlib.metadata as md

from tomo.protocols.protocol_base import ProtTomoBase
from tomo.objects import Tomogram, SetOfTomograms, SetOfTiltSeries
from pyworkflow import BETA
from tigre import Plugin

EXT_MRC = '.mrc'
EXT_MRCS_TS_EVEN_NAME = "_even.mrcs"
EXT_MRCS_TS_ODD_NAME = "_odd.mrcs"
EXT_MRC_EVEN_NAME = "_even.mrc"
EXT_MRC_ODD_NAME = "_odd.mrc"
OUTPUT_TOMOGRAMS_NAME = "Tomograms"

class TigreOutputs(Enum):
    tomograms = SetOfTomograms

class ProtTigreReconstruction(EMProtocol, ProtTigreForm):
    """
    'This program provides a variety of algorithms to reconstruct tomogram from a set of tilt series.\n'
    'The program will make use of a tilt series file (--tiltseries) which contains the tilt images and a set of angles given in a metadata file (--angles).\n'
    'The reconstruction algorithms can be selectec with (--method) and are grouped in families:\n\n'
    '\t F1 - Exact algorithms:\n'
    '\t\t 1) WBP or FBP    - (default) Weighted back projection or filtered back projection.\n'
    '\t\t 2) FDK           - Feldkamp-Davis-Kress algorithm.\n'
    '\t\t\t\t  See L. A. Feldkamp, L. C. Davis, and J. W. Kress,  Practical cone-beam algorithm. JOSA (1984)\n\n'


    '\t F2 - Gradient-based algorithms:\n'
    '\t\t 3) SIRT          - Simultaneus iterative reconstruction technique.\n'
    '\t\t\t\t  See P. Gilbert, Iterative methods for the three-dimensional reconstruction of an object from projections, Journal of Theoretical Biology, 35,1, 105-107, (1972)\n'
    '\t\t\t\t  See A.C. Kak and M. Slaney, Principles of Computerized Tomographic Imaging, Society of Industrial and Applied Mathematics, (2001)'
    '\t\t 4) SART          - Simultaneous algebraic reconstruction technique\n'
    '\t\t\t\t  See A. Andersen, A. Kak, Simultaneous Algebraic Reconstruction Technique (SART): A Superior Implementation of ART. Ultrasonic Imaging. 6,1, 81-94 (1984)\n'
    '\t\t 5) OSSART        - Ordered Subset Simultaneous Algebraic Reconstruction Technique\n'
    '\t\t\t\t  See Y. Censor and T. Elfving. Block-iterative algorithms with diagonally scaled oblique projections for the linear feasibility problem SIAM J. Matrix Anal. Appl. 24 40–58 (2002)\n.'
    '\t\t\t\t  See Ge, and J. Ming,  Ordered-subset simultaneous algebraic reconstruction techniques (OS-SART), Journal of X-Ray Science and Technology, 12, 3, 169-177, (2004)\n'
    '\t\t 6) ASD-POCS      - Adaptative Steepest Descent with Projections Onto the Convex Set.\n'
    '\t\t\t\t  See E.Y. Sidky, X. Pan, Image reconstruction in circular cone-beam computed tomography by constrained, total-variation minimization. Phys Med Biol, 53, 4777 (2008).\n'
    '\t\t 7) OS-ASD-POCS   - Oriented Subsets version of ASD_POCS\n'
    '\t\t\t\t  See E.Y. Sidky, X. Pan, Image reconstruction in circular cone-beam computed tomography by constrained, total-variation minimization. Phys Med Biol, 53, 4777 (2008).\n'
    '\t\t 8) PCSD          - Projection-controlled steepest descent method\n'
    '\t\t\t\t  See L. Liu, W. Lin, M. Jin. Reconstruction of sparse-view X-ray computed tomography using adaptive iterative algorithms.Computers in Biology and Medicine, 56, 97-106 (2015).\n'
    '\t\t 9) AwPCSD        - Adaptive Weighted version of PCSD\n'
    '\t\t\t\t  See M. Lohvithee, A. Biguri, M. Soleimani, Parameter selection in limited data cone-beam CT reconstruction using edge-preserving total variation algorithms, Physics in Medicine and Biology, (2017).\n'
    '\t\t 10) AwASD-POCS   - Adaptive Weighted TV (edge preserving) version of ASD-POCS\n\n'

    '\t F3 - Krylov subspace algorithms:\n'
    '\t\t 11) CGLS         - Conjugate Gradient Least Squares (CGLS) algorithm\n'
    '\t\t\t\t  See A. Bjorck, Numerical Methods for Least Squares Problems, SIAM, Philadelphia, (1996)\n'
    '\t\t 12) LSQR         - Least Squares QR factorization\n'
    '\t\t\t\t  See C.C. Paige and M.A. Saunders, A robust incomplete factorization preconditioner for positive deﬁnite matrices, ACM Trans. Math. Software, 8 (1982)\n'
    '\t\t 13) LSMR         - Least Squares MINRES\n'
    '\t\t\t\t  See D.C.L. Fong, and  M. Saunders, LSMR: an iterative algorithm for sparse least-squares problems, SIAM J. Sci. Comput. 33, 2950–2991, (2011)\n'
    '\t\t 14) hybrid-LSQR  - Hybrid Least Squares QR factorization\n'
    '\t\t\t\t  See reference is needed\n'
    '\t\t 15) AB-GMRES     - AB-generalized minimal residual method\n'
    '\t\t\t\t  See K. Hayami, J.F. Yin, and T. Ito, GMRES methods for least squares problems, SIAM J. Matrix Anal. Appl., 31, 2400–2430 (2010)\n'
    '\t\t\t\t  See Per c. Hansen, K. Hayami, K. Morikuni, GMRES methods for tomographic reconstruction with an unmatched back projector, Journal of Computational and Applied Mathematics, 413, 114352, (2022)\n'
    '\t\t 16) BA-GMRES     - BA-generalized minimal residual method\n'
    '\t\t\t\t  See K. Hayami, J.F. Yin, and T. Ito, GMRES methods for least squares problems, SIAM J. Matrix Anal. Appl., 31, 2400–2430 (2010)\n'
    '\t\t\t\t  See Per c. Hansen, K. Hayami, K. Morikuni, GMRES methods for tomographic reconstruction with an unmatched back projector, Journal of Computational and Applied Mathematics, 413, 114352, (2022)\n'
    '\t\t 17) IRN-TV-CGLS  - Using a binary mask\n'
    '\t\t\t\t  See S. Gazzola, M.E. Kilmer, J.G Nagy, O.Semerci and E.L Miller, An inner–outer iterative method for edge preservation in image restoration and reconstruction. Inverse Problems, 36(2020)\n\n'

    '\t F4 - Statistical reconstruction:\n'
    '\t\t 18) MLEM         - Maximum likelihood Expectation maximization\n'
    '\t\t\t\t See J. M. Ollinger. Maximum-likelihood reconstruction of transmission images in emission computed tomography via the EM algorithm, IEEE Transactions on Medical Imaging, 12, 89-101, (1994)\n\n'

    '\t F5 - Variational methods:\n'
    '\t\t 19) FISTA        - Fast iterative shrinkage thresholding algorithm.\n'
    '\t\t\t\t  See Q. Xu, D. Yang, J. Tan, A. Sawatzky, M.A. Anastasio, Accelerated fast iterative shrinkage thresholding algorithms for sparsity-regularized cone-beam CT image reconstruction. Med Phys. 43,4, 1849-1872 (2016)\n'
    '\t\t 20) SART-TV      - Simultaneous algebraic reconstruction technique with total variation\n'
    '\t\t\t\t  See  H. Yu and G. Wang, A soft-threshold filtering approach for reconstruction from a limited number of projections, Phys. Med. Biol. 55, 13, 3905–3916 (2010)\n').
    """

    _label = 'tigre reconstruction'
    _lastUpdateVersion = VERSION_3_0
    _devStatus = BETA
    _possibleOutputs = TigreOutputs


    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.Tomograms = None

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputSetOfTiltSeries', PointerParam, pointerClass='SetOfTiltSeries',
                      label="Tilt Series", 
                      important=True,
                      help='Select the Set of Tilt Series that will be used to reconstruct the tomograms.')

        form.addParam('tomoThickness',
                      IntParam,
                      default=300,
                      label='Tomogram thickness (voxels)',
                      important=True,
                      help='Size in voxels of the tomogram along the z axis (beam direction).')

        self._defineTigreParams(form)

        form.addHidden('useTigreInterpolation',
                      BooleanParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=False,
                      label='Use tigre interpolation?',
                      help= 'True if the interpolation is carried out by native Tigre. ' \
                            'False if the interpolation is carried out with the wrapper ' \
                            'that Tigre has (but not native) for applying alignment')

        form.addParam('processOddEven',
                      BooleanParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=True,
                      label='Reconstruct odd/even?',
                      help='If True, the full tilt series and the associated odd/even tilt series will be reconstructed. '
                           'The alignment applied to the odd/even tilt series will be exactly the same.')


        form.addHidden(GPU_LIST,
                       StringParam,
                       default='0',
                       label="Choose GPU IDs",
                       help="GPU ID. To pick the best available one set 0. "
                            "For a specific GPU set its number ID "
                            "(starting from 1).")


    # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            self._insertFunctionStep(self.reconstructTomogramStep, tsId)
            self._insertFunctionStep(self.createOutputStep, tsId)
        self._insertFunctionStep(self.closeOutputSetsStep)


    def getReconstructionMethod(self):
        family = self.family.get()
        if family == self.FAMILY_EXACT:
            recMethod = self.exactsMethod.get()
            args = ' --method %s ' % self.ALGORITHMS_EXACT[recMethod]
            args += ' --filter %s ' % self.FILTER_LIST[self.filter.get()]
        elif family == self.FAMILY_GRADIENT:
            recMethod = self.gradientMethod.get()
            args = ' --method %s ' % self.ALGORITHMS_GRADIENT[recMethod]

            if recMethod == self.SART:
                args += ' --iter %i ' % self.SARTiter.get()
                args += ' --lambda %f ' % self.SARTlambda.get()
                args += ' --lambdared %f ' % self.SARTlambdared.get()
            elif recMethod == self.SIRT:
                args += ' --iter %i ' % self.SIRTiter.get()
                args += ' --lambda %f ' % self.SIRTlambda.get()
                args += ' --lambdared %f ' % self.SIRTlambdared.get()
            elif recMethod == self.OSSART:
                args += ' --iter %i ' % self.OSSARTiter.get()
                args += ' --lambda %f ' % self.OSSARTlambda.get()
                args += ' --lambdared %f ' % self.OSSARTlambdared.get()
            elif recMethod == self.ASDPOCS:
                args += ' --iter %i ' % self.ASDPOCSiter.get()
                args += ' --tviter %i ' % self.ASDPOCStviter.get()
                args += ' --alpha %f ' % self.ASDPOCSalpha.get()
                args += ' --lambda %f ' % self.ASDPOCSlambda.get()
                args += ' --lambdared %f ' % self.ASDPOCSlambdared.get()
                args += ' --alpha_red %f ' % self.ASDPOCSalphared.get()
                args += ' --ratio %f ' % self.ASDPOCSratio.get()
            elif recMethod == self.OSASDPOCS:
                args += ' --iter %i ' % self.OSASDPOCSiter.get()
                args += ' --tviter %i ' % self.OSASDPOCStviter.get()
                args += ' --alpha %f ' % self.OSASDPOCSalpha.get()
                args += ' --lambda %f ' % self.OSASDPOCSlambda.get()
                args += ' --lambdared %f ' % self.OSASDPOCSlambdared.get()
                args += ' --alpha_red %f ' % self.OSASDPOCSalphared.get()
                args += ' --ratio %f ' % self.OSASDPOCSratio.get()
            elif recMethod == self.PCSD:
                args += ' --iter %i ' % self.PCSDiter.get()
            elif recMethod == self.AWPCSD:
                args += ' --iter %i ' % self.AWPCSDiter.get()
            elif recMethod == self.AWASDPOCS:
                args += ' --iter %i ' % self.AWASDPOCSiter.get()
                args += ' --tviter %i ' % self.AWASDPOCStviter.get()
                args += ' --alpha %f ' % self.AWASDPOCSalpha.get()
                args += ' --lambda %f ' % self.AWASDPOCSlambda.get()
                args += ' --lambdared %f ' % self.AWASDPOCSlambdared.get()
                args += ' --alpha_red %f ' % self.AWASDPOCSalphared.get()
                args += ' --ratio %f ' % self.AWASDPOCSratio.get()
            else:
                raise Exception('Combination of form parameter failed. Contact Scipion developers')
        elif family == self.FAMILY_KRYLOV:
            recMethod = self.KrylovMethod.get()
            args = ' --method %s ' % self.ALGORITHMS_KRYLOV[recMethod]
            if recMethod == self.CGLS:
                args += ' --iter %i ' % self.CGLSiter.get()
            elif recMethod == self.LSQR:
                args += ' --iter %i ' % self.LSQRiter.get()
            elif recMethod == self.LSMR:
                args += ' --iter %i ' % self.LSMRiter.get()
                args += ' --lambda %f ' % self.LSMRlambda.get()
            elif recMethod == self.HYBRIDLSQR:
                args += ' --iter %i ' % self.hybridLSQRiter.get()
            elif recMethod == self.ABGMRES:
                args += ' --iter %i ' % self.ABGMRESiter.get()
            elif recMethod == self.BAGMRES:
                args += ' --iter %i ' % self.BAGMRESiter.get()
            elif recMethod == self.IRNTVCGLS:
                args += ' --iter %i ' % self.IRNTVCGLSiter.get()
                args += ' --lambda %f ' % self.IRNTVCGLStvlambda.get()
                args += ' --innerouter %i ' % self.IRNTVCGLSinnerouter.get()
            else:
                raise Exception('Combination of form parameter failed. Contact Scipion developers')
        elif family == self.FAMILY_STATISTICAL:
            recMethod = self.MLEM
            args = ' --method %s ' % self.ALGORITHMS_STATISTICAL[recMethod]
            args += ' --iter %i ' % self.MLEMiter.get()
        elif family == self.FAMILY_VARIATIONAL:
            recMethod = self.varMethod.get()
            args = ' --method %s ' % self.ALGORITHMS_VARIATIONAL[recMethod]
            if recMethod == self.FISTA:
                args += ' --iter %i ' % self.FISTAiter.get()
                args += ' --tviter %i ' % self.FISTAtviter.get()
                args += ' --tvlambda %f ' % self.FISTAtvlambda.get()
            elif recMethod == self.SARTTV:
                args += ' --iter %i ' % self.SARTTViter.get()
                args += ' --tviter %i ' % self.SARTTVtviter.get()
                args += ' --tvlambda %f ' % self.SARTTVtvlambda.get()
            else:
                raise Exception('Combination of form parameter failed. Contact Scipion developers')
        else:
            raise Exception('Reconstruction method not properly selected.')
        return recMethod, args


    def reconstructTomogramStep(self, tsId):
        """
        This function computes the reconstructed tomogram
        """
        ts = self.inputSetOfTiltSeries.get()[{'_tsId': tsId}]

        # Defining the output folder
        tomoPath = self._getExtraPath(tsId)
        pwutils.path.makePath(tomoPath)

        fnAngles = os.path.join(tomoPath, tsId+'.tlt')
        ts.generateTltFile(fnAngles)

        firstItem = ts.getFirstItem()
        tmpPrefix = self._getTmpPath(tsId)
        fnTs = os.path.join(tmpPrefix, firstItem.parseFileName())
        fnTs = firstItem.getFileName()

        # Defining outfiles
        fullTomogramName = os.path.join(tomoPath, tsId+EXT_MRC)

        recMethod, args = self.getReconstructionMethod()
        fnXf = self._getExtraPath(tsId,tsId+'.xf')
        ts.writeXfFile(fnXf)

        paramsTS = ' --tiltseries %s' % fnTs
        otherParams = ' --angles %s ' % fnAngles
        otherParams += ' --xf %s ' % fnXf
        otherParams += ' --thickness %i ' % self.tomoThickness.get()
        otherParams += args
        otherParams += ' --gpu %s' % self.gpuList.get()
        paramsOut = ' -o %s' % fullTomogramName


        programTigre = Plugin.getTigreProgram('tigre_reconstruction.py')
        tigreArgs = programTigre + paramsTS + otherParams + paramsOut

        Plugin.runTigre(self, tigreArgs)

        if ts.hasOddEven() and self.processOddEven.get():
            fnOdd = ts.getOddFileName()
            paramsTS = ' --tiltseries %s' % fnOdd
            paramsOut = ' -o %s' % os.path.join(tomoPath, tsId+EXT_MRC_ODD_NAME)
            tigreArgs =  programTigre + paramsTS + otherParams + paramsOut
            Plugin.runTigre(self, tigreArgs)

            fnEven = ts.getEvenFileName()
            paramsTS = ' --tiltseries %s' % fnEven
            paramsOut = ' -o %s' % os.path.join(tomoPath, tsId+EXT_MRC_EVEN_NAME)
            tigreArgs =  programTigre + paramsTS + otherParams + paramsOut
            Plugin.runTigre(self, tigreArgs)

    def getOutputSetOfTomograms(self, inputSet, binning=1) -> SetOfTomograms:

        if self.Tomograms:
            getattr(self, OUTPUT_TOMOGRAMS_NAME).enableAppend()
        else:
            outputSetOfTomograms =  self._createSetOfTomograms()

            outputSetOfTomograms.setAcquisition(inputSet.getAcquisition())
            outputSetOfTomograms.setSamplingRate(inputSet.getSamplingRate())


            outputSetOfTomograms.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{OUTPUT_TOMOGRAMS_NAME: outputSetOfTomograms})
            self._defineSourceRelation(inputSet, outputSetOfTomograms)

        return self.Tomograms

    def setTomoOddEven(self, tsId: str, outTomo: Tomogram) -> None:
        if self.processOddEven.get():
            tomoPath = self._getExtraPath(tsId)
            halfMapsList = [os.path.join(tomoPath, tsId+EXT_MRC_ODD_NAME),
                            os.path.join(tomoPath, tsId+EXT_MRC_EVEN_NAME)]
            outTomo.setHalfMaps(halfMapsList)
        else:
            outTomo.setHalfMaps([])

    def createOutputStep(self, tsId):
        ts = self.inputSetOfTiltSeries.get()[{'_tsId': tsId}]
        tsId = ts.getTsId()

        fullTomogramName = self._getExtraPath(tsId, tsId+EXT_MRC)

        if os.path.exists(fullTomogramName):
            output = self.getOutputSetOfTomograms(self.inputSetOfTiltSeries.get())

            newTomogram = Tomogram()
            newTomogram.setLocation(fullTomogramName)

            if ts.hasOddEven() and self.processOddEven.get():
                self.setTomoOddEven(tsId, newTomogram)

            newTomogram.setTsId(tsId)
            newTomogram.setSamplingRate(ts.getSamplingRate())

            # Set default tomogram origin
            newTomogram.setOrigin(newOrigin=None)
            newTomogram.setAcquisition(ts.getAcquisition())

            output.append(newTomogram)
            output.update(newTomogram)
            output.write()
            self._store()

    def closeOutputSetsStep(self):
        for _, output in self.iterOutputAttributes():
            output.setStreamState(Set.STREAM_CLOSED)
            output.write()
        self._store()

    # --------------------------- INFO functions ------------------------------

    def _methods(self):
        messages = []
        return messages

    def _summary(self):
        summary = []

        return summary

    def _citations(self):
        return []
