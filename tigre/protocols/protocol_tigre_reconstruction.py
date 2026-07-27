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
from pyworkflow.protocol.params import (PointerParam, IntParam, BooleanParam, LabelParam, EnumParam, StringParam, LEVEL_ADVANCED, GPU_LIST)
import pyworkflow.utils.path as path

from pwem.protocols import EMProtocol
from pwem.emlib import lib
import pwem.emlib.metadata as md

from tomo.protocols.protocol_base import ProtTomoBase
from tomo.objects import Tomogram, SetOfTomograms, SetOfTiltSeries
#from .protocol_base import xTomoBase, OUTPUT_TOMOGRAMS_NAME
from pyworkflow import BETA
from tigre import Plugin
#from xtomo.utils import calculateRotationAngleAndShiftsFromTM

EXT_MRC = '.mrc'
EXT_MRCS_TS_EVEN_NAME = "_even.mrcs"
EXT_MRCS_TS_ODD_NAME = "_odd.mrcs"
EXT_MRC_EVEN_NAME = "_even.mrc"
EXT_MRC_ODD_NAME = "_odd.mrc"
OUTPUT_TOMOGRAMS_NAME = "Tomograms"

class TigreOutputs(Enum):
    tomograms = SetOfTomograms

class ProtTigreReconstruction(EMProtocol, ProtTomoBase):
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
        form.addSection(label='Input')

        form.addParam('inputSetOfTiltSeries', PointerParam, pointerClass='SetOfTiltSeries',
                      label="Tilt Series", 
                      important=True,
                      help='Select the Set of Tilt Series that will be used to reconstruct the tomograms.')

        form.addParam('family',
                      EnumParam,
                      choices=self.FAMILIES,
                      default=self.FAMILY_EXACT,
                      label="Family of methods",
                      isplay=EnumParam.DISPLAY_COMBO,
                      help='Select an option to reconstruct tomograms: \n '
                           '_ART_: Arithmetic reconstruction technique. \n'
                           '_SIRT_: (only with MPI) Simultaneous Iterative Reconstruction Technique. \n ')

        #TODO: Fix documentation and references to the algorithms
        form.addParam('exactsMethod',
                      EnumParam,
                      choices=self.ALGORITHMS_EXACT,
                      default=self.WBP,
                      condition='family==%d' % self.FAMILY_EXACT,
                      isplay=EnumParam.DISPLAY_COMBO,
                      label="Reconstruction Algorithm",
                      help='Select an option to reconstruct tomograms: \n '
                           '_WBP_: (default) Weighted back projection or filtered back projection. \n'
                           '_FDK_: Feldkamp-Davis-Kress algorithm. \n')

        #TODO: Fix documentation and references to the algorithms
        form.addParam('gradientMethod',
                      EnumParam,
                      choices=self.ALGORITHMS_GRADIENT,
                      default=self.SART,
                      condition='family==%d' % self.FAMILY_GRADIENT,
                      isplay=EnumParam.DISPLAY_COMBO,
                      label="Reconstruction Algorithm",
                      help='Select an option to reconstruct tomograms: \n '
                           '_SART_: (default) Simultaneous algebraic reconstruction technique. \n'
                           '_SIRT_: Simultaneus iterative reconstruction technique. \n'
                           '_OS-SART_: Ordered Subset Simultaneous Algebraic Reconstruction Technique. \n'
                           '_ASD-POCS_: Adaptative Steepest Descent with Projections Onto the Convex Set. \n'
                           '_OS-ASD-POCS_: Oriented Subsets version of ASD-POC. \n'
                           '_PCSD_: Projection-controlled steepest descent method. \n'
                           '_AW-PCSD_: Adaptive Weighted version of PCSD. \n'
                           '_AW-ASD-POCS_: Adaptive Weighted TV (edge preserving) version of ASD-POCS. \n')

        #TODO: Fix documentation and references to the algorithms
        form.addParam('KrylovMethod',
                      EnumParam,
                      choices=self.ALGORITHMS_KRYLOV,
                      default=self.CGLS,
                      condition='family==%d' % self.FAMILY_KRYLOV,
                      isplay=EnumParam.DISPLAY_COMBO,
                      label="Reconstruction Algorithm",
                      help='Select an option to reconstruct tomograms: \n '
                           '_CGLS_: (default) Conjugate Gradient Least Squares. \n'
                           '_LSQR_: Least Squares QR factorization. \n'
                           '_LSMR_: Least Squares MINRES. \n'
                           '_hybrid LSQR_: Hybrid Least Squares QR factorization. \n'
                           '_AB-GMRES_: AB-generalized minimal residual method. \n'
                           '_BA-GMRES_: BA-generalized minimal residual method. \n'
                           '_IRN-TV-CGLS_: Total variation CGLS. \n')

        #TODO: Fix documentation and references to the algorithms
        form.addParam('statsMethod',
                      LabelParam,
                      condition='family==%d' % self.FAMILY_STATISTICAL,
                      label="MLEM algorithm",
                      help='Select an option to reconstruct tomograms: \n '
                           '_MLEM_: (default) Maximum likelihood Expectation maximization. \n')

        #TODO: Fix documentation and references to the algorithms
        form.addParam('varMethod',
                      EnumParam,
                      choices=self.ALGORITHMS_VARIATIONAL,
                      default=self.FISTA,
                      condition='family==%d' % self.FAMILY_VARIATIONAL,
                      isplay=EnumParam.DISPLAY_COMBO,
                      label="Reconstruction Algorithm",
                      help='Select an option to reconstruct tomograms: \n '
                           '_FISTA_: (default) Fast iterative shrinkage thresholding algorithm.. \n'
                           '_SART-TV_: Simultaneous algebraic reconstruction technique with total variation. \n')

        filterCondition = 'family==%d' % self.FAMILY_EXACT

        form.addParam('filter',
                      EnumParam,
                      choices=self.FILTER_LIST,
                      default=self.RAMLAK,
                      condition=filterCondition,
                      isplay=EnumParam.DISPLAY_COMBO,
                      label="Filter",
                      help='Select the filter to be applied: \n '
                           '_Ram-lak_: (default) Fast iterative shrinkage thresholding algorithm.. \n'
                           '_Shepp-Logan_: Simultaneous algebraic reconstruction technique with total variation. \n'
                           '_Cosine_: Simultaneous algebraic reconstruction technique with total variation. \n'
                           '_Hamming_: Simultaneous algebraic reconstruction technique with total variation. \n'
                           '_Hann_: Simultaneous algebraic reconstruction technique with total variation. \n')

        form.addParam('tomoThickness',
                      IntParam,
                      default=300,
                      label='Tomogram thickness (voxels)',
                      important=True,
                      help='Size in voxels of the tomogram in the z axis (beam direction).')

        form.addParam('iter',
                      IntParam,
                      condition='not family==%d' % self.FAMILY_EXACT,
                      allowsNull=True,
                      label='Iterations',
                      help='Number of iterations of the reconstruction algorithm. The wizard will suggest recommended'
                           'values according to the selected algorithm. The recommended values are the next ones: \n'
                           '_SIRT_: 20 iterations \n'
                           '_SART_: 20 iterations \n '
                           '_OS-SART_: 20 iterations \n '
                           '_PCSD_: 20 iterations \n '
                           '_AW-PCSD_: 20 iterations \n '
                           '_FISTA_: 100 iterations \n '
                           '_SART-TV_: 100 iterations \n '
                           '_MLEM_: 500 iterations \n ')

        form.addParam('useTigreInterpolation',
                      BooleanParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=True,
                      label='Use tigre interpolation?',
                      help='TODO')

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
            #self._insertFunctionStep(self.convertInputStep, tsId)
            self._insertFunctionStep(self.reconstructTomogramStep, tsId)
            self._insertFunctionStep(self.createOutputStep, tsId)
        self._insertFunctionStep(self.closeOutputSetsStep)


    def convertInputStep(self, tsId):
        # Considering swapXY is required to make tilt axis vertical
        #TODO: Odd-even not checked
        oddEvenFlag = False
        if self.inputSetOfTiltSeries.get().hasOddEven() and self.processOddEven.get():
            oddEvenFlag = True

        super().convertInputStep(tsObjId, doSwap=True, oddEven=oddEvenFlag)

    def getReconstructionMethod(self):

        family = self.family.get()
        if family == self.FAMILY_EXACT:
            recMethod = self.exactsMethod.get()
            args = ' --method %s ' % self.ALGORITHMS_EXACT[recMethod]
            args += ' --filter %s ' % self.FILTER_LIST[self.filter.get()]
            #TODO: think the use of parkerweights
            #args += ' --parker '
        elif family == self.FAMILY_GRADIENT:
            recMethod = self.gradientMethod.get()
            args = ' --method %s ' % self.ALGORITHMS_GRADIENT[recMethod]
            args += ' --iter %i ' % self.getIter()
        elif family == self.FAMILY_KRYLOV:
            recMethod = self.KrylovMethod.get()
            args = ' --method %s ' % self.ALGORITHMS_KRYLOV[recMethod]
            args += ' --iter %i ' % self.getIter()
        elif family == self.FAMILY_STATISTICAL:
            recMethod = self.MLEM
            args = ' --method %s ' % self.ALGORITHMS_STATISTICAL[recMethod]
            args += ' --iter %i ' % self.getIter()
        elif family == self.FAMILY_VARIATIONAL:
            recMethod = self.varMethod.get()
            args = ' --method %s ' % self.ALGORITHMS_VARIATIONAL[recMethod]
            args += ' --iter %i ' % self.getIter()
        else:
            raise Exception('Reconstruction method not properly selected.')
        return recMethod, args

    def getIter(self):
        family = self.family.get()
        if not self.iter.hasValue():
            if family == self.FAMILY_GRADIENT:
                if self.gradientMethod.get() == self.PCSD or self.gradientMethod.get() == self.PCSD:
                    niters = 25
                else:
                    niters = 20
                return niters
            if family == self.FAMILY_KRYLOV:
                niters = 20
                return niters
            if family == self.FAMILY_STATISTICAL:
                niters = 500
                return niters
            if family == self.FAMILY_VARIATIONAL:
                niters = 200
                return niters
        else:
            return self.iter.get()
        
    def generateAlignmentFile(self, ts, fnAngles):
        mdAli = lib.MetaData()
        for ti in ts:
            newRow = md.Row()
            tilt = ti.getTiltAngle()
            rot, sx, sy = calculateRotationAngleAndShiftsFromTM(ti)
            newRow.setValue(lib.MDL_ANGLE_TILT, tilt)
            newRow.setValue(lib.MDL_ANGLE_ROT, rot)
            newRow.setValue(lib.MDL_SHIFT_X, sx)
            newRow.setValue(lib.MDL_SHIFT_Y, sy)
            newRow.addToMd(mdAli)

        mdAli.write(fnAngles)


    def reconstructTomogramStep(self, tsId):
        '''
        This function computes the reconstructed tomogram
        '''
        ts = self.inputSetOfTiltSeries.get()[{'_tsId': tsId}]

        #Defining the output folder
        tomoPath = self._getExtraPath(tsId)
        pwutils.path.makePath(tomoPath)

        fnAngles = os.path.join(tomoPath, tsId+'.tlt')
        ts.generateTltFile(fnAngles)

        # Next lines are commented on purpose. Tigre allows an internal
        # interpolation. In our hands it did not work, but the commented
        # lines prepare the data to use it
        # fnAngles = os.path.join(tomoPath, tsId+'.xmd')
        # self.generateAlignmentFile(ts, fnAngles)

        firstItem = ts.getFirstItem()
        tmpPrefix = self._getTmpPath(tsId)
        fnTs = os.path.join(tmpPrefix, firstItem.parseFileName())
        fnTs = firstItem.getFileName()
        #Defining outfiles
        fullTomogramName = os.path.join(tomoPath, tsId+EXT_MRC)

        recMethod, args = self.getReconstructionMethod()

        params = ' --tiltseries %s' % fnTs
        params += ' --angles %s ' % fnAngles
        params += ' --thickness %i ' % self.tomoThickness.get()
        params += ' -o %s' % fullTomogramName
        params += ' --normalize standard'
        params += args
        params += ' --gpu %s' % self.gpuList.get()

        programTigre = '/home/tomo/tigreBin/tigre/tigre_reconstruction.py'

        Plugin.runTigre(self, f' python3 {programTigre} ', params)


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

    def createOutputStep(self, tsId):
        ts = self.inputSetOfTiltSeries.get()[{'_tsId': tsId}]
        tsId = ts.getTsId()

        fullTomogramName = self._getExtraPath(tsId, tsId+EXT_MRC)

        if os.path.exists(fullTomogramName):
            output = self.getOutputSetOfTomograms(self.inputSetOfTiltSeries.get())

            newTomogram = Tomogram()
            newTomogram.setLocation(fullTomogramName)

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
        if hasattr(self, 'outputSetOfTomograms'):
            messages.append('')
        return messages

    def _summary(self):
        summary = []

        return summary

    def _citations(self):
        return []
