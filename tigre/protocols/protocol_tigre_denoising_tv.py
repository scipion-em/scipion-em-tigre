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

from pyworkflow import VERSION_3_0
from pyworkflow.object import Set, Float
from pyworkflow.protocol.params import (PointerParam, IntParam, FloatParam, StringParam, LEVEL_ADVANCED, GPU_LIST)
import pyworkflow.utils.path as path

from pwem.protocols import EMProtocol

from tomo.protocols.protocol_base import ProtTomoBase
from tomo.objects import Tomogram, SetOfTomograms
from pyworkflow import BETA
from tigre import Plugin

MRCEXT = '.mrc'
OUTPUT_TOMOGRAMS_NAME = "Tomograms"


class TigreDenoisingTvOutputs(Enum):
    tomograms = SetOfTomograms


class ProtTigreDenoisingTV(EMProtocol, ProtTomoBase):
    """
    This protocols carries out a denoising based on total variation
    """

    _label = 'denosing tv'
    _lastUpdateVersion = VERSION_3_0
    _devStatus = BETA
    _possibleOutputs = {TigreDenoisingTvOutputs}


    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.Tomograms = None
    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('tomograms', PointerParam, pointerClass='SetOfTomograms',
                      label="Tomograms", important=True,
                      help='Select the Set of Tomograms to be denoised.')

        form.addParam('iters',
                      IntParam,
                      default=50,
                      label='Iterations',
                      help='Number of total variation iterations')

        form.addParam('lamda',
                      FloatParam,
                      expeexpertLevel=LEVEL_ADVANCED,
                      default=15.0,
                      label='Lambda yperparameter',
                      help='Hyperparameter. The update will be multiplied by this number every iteration,'
                           'to make the steps bigger or smaller. Default: lmbda=15.0')

        form.addHidden(GPU_LIST,
                       StringParam,
                       default='0',
                       label="Choose GPU IDs",
                       help="GPU ID. To pick the best available one set 0. "
                            "For a specific GPU set its number ID "
                            "(starting from 1).")

    # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        for ts in self.tomograms.get():
            tsId = ts.getTsId()
            self._insertFunctionStep(self.denoisingStep, tsId)
            self._insertFunctionStep(self.createOutputStep, tsId)
        self._insertFunctionStep(self.closeOutputSetsStep)


    def denoisingStep(self, tsId):
        '''
        This function computes the reconstructed tomogram
        '''
        tomo = self.tomograms.get()[{'_tsId': tsId}]

        #Defining the output folder
        tomoPath = self._getExtraPath(tsId)
        path.makePath(tomoPath)

        fnTs = tomo.getFileName()

        #Defining outfiles
        fullTomogramName = os.path.join(tomoPath, tsId+MRCEXT)

        params = ' -i %s' % fnTs
        params += ' --iters %i ' % self.iters.get()
        params += ' --lambda %f ' % self.lamda.get()
        params += ' -o %s' % fullTomogramName
        params += ' --gpu %s' % self.gpuList.get()

        programTigre = '/home/tomo/tigreBin/tigre/tigre_denoising_tv.py'

        Plugin.runTigre(self, f' python3 {programTigre} ', params)

    def getOutputSetOfTomograms(self, inputSet, binning=1) -> SetOfTomograms:

        if self.Tomograms:
            getattr(self, OUTPUT_TOMOGRAMS_NAME).enableAppend()
        else:
            outputSetOfTomograms =  self._createSetOfTomograms()
            outputSetOfTomograms.copyInfo(inputSet)
            outputSetOfTomograms.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{OUTPUT_TOMOGRAMS_NAME: outputSetOfTomograms})
            self._defineSourceRelation(inputSet, outputSetOfTomograms)

        return self.Tomograms

    def createOutputStep(self, tsId):
        ts = self.tomograms.get()[{'_tsId': tsId}]

        fullTomogramName = os.path.join(self._getExtraPath(tsId), tsId+MRCEXT)

        if os.path.exists(fullTomogramName):
            output = self.getOutputSetOfTomograms(self.tomograms.get())

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