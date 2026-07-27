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
import unittest

from pyworkflow.tests import BaseTest, setupTestProject
from pyworkflow.tests import DataSet

import tomo
from pyworkflow.utils import weakImport
from tomo.tests import RE4_STA_TUTO
from tigre.protocols.protocol_tigre_denoising_tv import ProtTigreDenoisingTV, OUTPUT_TOMOGRAMS_NAME
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer
from tomo.utils import existsPlugin



@unittest.skipIf(not existsPlugin('imod'), 'IMOD plugin not available')
class TestTigreProtDenoisingTV(TestBaseCentralizedLayer):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.inputDataSet = DataSet.getDataSet(RE4_STA_TUTO)

    def importTomograms(test: BaseTest, pathToData, pattern='*.mrc', samplingRate=1.0, label='Import Tomograms'):
        protTomograms = test.newProtocol(tomo.protocols.protocol_import_tomograms.ProtImportTomograms,
                                         objLabel=label,
                                         filesPath=pathToData,
                                         filesPattern=pattern,
                                         samplingRate=samplingRate)
        test.launchProtocol(protTomograms)
        outputTomograms = getattr(protTomograms, tomo.protocols.protocol_import_tomograms.OUTPUT_NAME, None)
        test.assertIsNotNone(outputTomograms, "Test failed importing tomograms. This is a scipion-em-tomo problem")
        return outputTomograms



    def runTigreDenoisingTV(test: BaseTest, tomos, **kwargs):
        protTigreDenoising = test.newProtocol(ProtTigreDenoisingTV,
                                              tomograms=tomos,
                                              **kwargs)
        test.launchProtocol(protTigreDenoising)
        return getattr(protTigreDenoising, OUTPUT_TOMOGRAMS_NAME, None)

    def testDenoisingTV(self):
        pathToData = self.inputDataSet.getPath()

        importedTomos = self.importTomograms(self, pathToData, pattern='*.mrc', samplingRate=1.0, label='Import Tomograms')

        protLabel = 'tigre tv denoising'
        protDenoising = self.runTigreDenoisingTV(self, importedTomos, objLabel=protLabel)
        self.assertSetSize(protDenoising, size=6,
                           msg='Denoising failed')


