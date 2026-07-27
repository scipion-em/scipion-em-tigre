# **************************************************************************
# *
# * Authors:    Jose Luis Vilas Prieto (jlvilas@cnb.csic.es)
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
from pyworkflow.tests import *
from pyworkflow.utils import weakImport, magentaStr
from pyworkflow.utils import path
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer
from tomo.tests import RE4_STA_TUTO, DataSetRe4STATuto
from tomo.protocols.protocol_ts_import import ProtImportTs

with weakImport("imod"):
    from imod.constants import OUTPUT_TILTSERIES_NAME
    from imod.protocols import ProtImodImportTransformationMatrix, ProtImodTsNormalization

from ..protocols.protocol_tigre_reconstruction import ProtTigreReconstruction

OUTPUT_TOMOGRAMS_NAME = "Tomograms"
TS_03 = 'TS_03'
TS_54 = 'TS_54'

class TestTigreProtReconstruction(TestBaseCentralizedLayer):
    unbinnedSRate = DataSetRe4STATuto.unbinnedPixSize.value
    nTomos = 2
    tomoWidth = 300
    bin4 = 4
    bin4SRate = DataSetRe4STATuto.unbinnedPixSize.value * 4
    tomoDimsThk300 = [960, 928, 300]
    excludedViewsDict = {
        TS_03: [0, 1, 2, 38, 39],
        TS_54: [0, 1, 38, 39, 40]
    }

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet(RE4_STA_TUTO)
        cls._runPreviousProtocols()

    @classmethod
    def _runPreviousProtocols(cls):
        cls.importedTs = cls._runImportTs()
        cls.tsWithAlignment = cls._runImportTrMatrix()
        cls.tsWithAliBin4 = cls._runTsPreprocess()

    @classmethod
    def _runImportTs(cls):
        print(magentaStr("\n==> Importing the tilt series:"))
        protImportTs = cls.newProtocol(ProtImportTs,
                                       filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                       filesPattern=DataSetRe4STATuto.tsPattern.value,
                                       exclusionWords=DataSetRe4STATuto.exclusionWordsTs03ts54.value,
                                       anglesFrom=2,  # From tlt file
                                       voltage=DataSetRe4STATuto.voltage.value,
                                       magnification=DataSetRe4STATuto.magnification.value,
                                       sphericalAberration=DataSetRe4STATuto.sphericalAb.value,
                                       amplitudeContrast=DataSetRe4STATuto.amplitudeContrast.value,
                                       samplingRate=cls.unbinnedSRate,
                                       doseInitial=DataSetRe4STATuto.initialDose.value,
                                       dosePerFrame=DataSetRe4STATuto.dosePerTiltImg.value,
                                       tiltAxisAngle=DataSetRe4STATuto.tiltAxisAngle.value)

        cls.launchProtocol(protImportTs)
        tsImported = getattr(protImportTs, 'outputTiltSeries', None)
        return tsImported

    @classmethod
    def _runImportTrMatrix(cls):
        print(magentaStr("\n==> Importing the TS' transformation matrices with IMOD:"))
        protImportTrMatrix = cls.newProtocol(ProtImodImportTransformationMatrix,
                                             filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                             filesPattern=DataSetRe4STATuto.transformPattern.value,
                                             inputSetOfTiltSeries=cls.importedTs)
        cls.launchProtocol(protImportTrMatrix)
        outTsSet = getattr(protImportTrMatrix, OUTPUT_TILTSERIES_NAME, None)
        return outTsSet

    @classmethod
    def _runTsPreprocess(cls):
        print(magentaStr(f"\n==> Binning the TS with using a binning factor of {cls.bin4}:"))
        protTsPreprocess = cls.newProtocol(ProtImodTsNormalization,
                                           inputSetOfTiltSeries=cls.tsWithAlignment,
                                           binning=cls.bin4)
        cls.launchProtocol(protTsPreprocess)
        outTsSet = getattr(protTsPreprocess, OUTPUT_TILTSERIES_NAME, None)
        return outTsSet

    @classmethod
    def runTigreReconstruction(cls, test: BaseTest, inputSetOfTiltSeries, **kwargs):
        protTigre = test.newProtocol(ProtTigreReconstruction,
                                     inputSetOfTiltSeries=inputSetOfTiltSeries,
                                     **kwargs)
        test.launchProtocol(protTigre)
        return getattr(protTigre, OUTPUT_TOMOGRAMS_NAME, None)

    def testExactReconstruction(self):
        for fidx, appliedFilter in enumerate(ProtTigreReconstruction.FILTER_LIST):
            for aidx, algorithm in enumerate(ProtTigreReconstruction.ALGORITHMS_EXACT):
                label = algorithm + ' ' + appliedFilter
                print(magentaStr(f"\n==> Reconstruction with {label}"))
                tigreTomogram = self.runTigreReconstruction(self, objLabel=label,
                                                       inputSetOfTiltSeries=self.tsWithAliBin4,
                                                       family=ProtTigreReconstruction.FAMILY_EXACT,
                                                       exactsMethod=aidx,
                                                       tomoThickness=self.tomoWidth,
                                                       filter=fidx)
                self._checkTomos(tigreTomogram)

    '''
    def testGradientReconstruction(self):
        for aidx, algorithm in enumerate(ProtTigreReconstruction.ALGORITHMS_GRADIENT):
            label = 'Gradient ' + algorithm
            print(magentaStr(f"\n==> Reconstruction with {label}"))
            tigreTomogram = self.runTigreReconstruction(self, objLabel=label,
                                                   inputSetOfTiltSeries=self.tsWithAliBin4,
                                                   family=ProtTigreReconstruction.FAMILY_GRADIENT,
                                                   gradientMethod=aidx,
                                                   tomoThickness=self.tomoWidth,
                                                   iter=20)
            self._checkTomos(tigreTomogram)
    '''
    @classmethod
    def testGradientReconstruction(cls, label, ts, recmethod, thickness):
        tigreTomogram = cls.runTigreReconstruction(cls, objLabel=label,
                                                   inputSetOfTiltSeries=ts,
                                                   family=ProtTigreReconstruction.FAMILY_GRADIENT,
                                                   gradientMethod=recmethod,
                                                   tomoThickness=thickness,
                                                   iter=20)
        cls._checkTomos(tigreTomogram)

    for aidx, algorithm in enumerate(ProtTigreReconstruction.ALGORITHMS_GRADIENT):
        label = 'Gradient ' + algorithm
        testGradientReconstruction(label, cls.tsWithAliBin4, aidx, tomoWidth)
        

    def testKrylovReconstruction(self):
        for aidx, algorithm in enumerate(ProtTigreReconstruction.ALGORITHMS_KRYLOV):
            label = 'Krylov ' + algorithm
            print(magentaStr(f"\n==> Reconstruction with {label}"))
            tigreTomogram = self.runTigreReconstruction(self, objLabel=label,
                                                   inputSetOfTiltSeries=self.tsWithAliBin4,
                                                   family=ProtTigreReconstruction.FAMILY_KRYLOV,
                                                   KrylovMethod=aidx,
                                                   tomoThickness=self.tomoWidth,
                                                   iter=20)
            self._checkTomos(tigreTomogram)


    def testStatisticalReconstruction(self):
        for aidx, algorithm in enumerate(ProtTigreReconstruction.ALGORITHMS_STATISTICAL):
            label = 'Statistical ' + algorithm
            print(magentaStr(f"\n==> Reconstruction with {label}"))
            tigreTomogram = self.runTigreReconstruction(self, objLabel=label,
                                                   inputSetOfTiltSeries=self.tsWithAliBin4,
                                                   family=ProtTigreReconstruction.FAMILY_STATISTICAL,
                                                   tomoThickness=self.tomoWidth,
                                                   iter=500)
            self._checkTomos(tigreTomogram)

    def testVariationalReconstruction(self):
        for aidx, algorithm in enumerate(ProtTigreReconstruction.ALGORITHMS_VARIATIONAL):
            label = 'Variational ' + algorithm
            print(magentaStr(f"\n==> Reconstruction with {label}"))
            tigreTomogram = self.runTigreReconstruction(self, objLabel=label,
                                                   inputSetOfTiltSeries=self.tsWithAliBin4,
                                                   family=ProtTigreReconstruction.FAMILY_VARIATIONAL,
                                                   varMethod=aidx,
                                                   tomoThickness=self.tomoWidth,
                                                   iter=100)
            self._checkTomos(tigreTomogram)

    def _checkTomos(self, inTomoSet):
        self.checkTomograms(inTomoSet,
                            expectedSetSize=self.nTomos,
                            expectedSRate=self.bin4SRate,
                            expectedDimensions=self.tomoDimsThk300)