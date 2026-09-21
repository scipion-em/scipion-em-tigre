# **************************************************************************
# *
# * Authors:     jlvilas (jlvilas@cnb.csic.es)
# *
# * Spanish Research Council
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
from os.path import join

import pyworkflow.utils as pwutils
from pyworkflow import TOMO
import pwem

from tigre.constants import *

__version__ = '0.0.0'
_logo = "icon.png"
#_references = ['you2019']


class Plugin(pwem.Plugin):
    _url = "https://github.com/scipion-em/scipion-em-tigre"
    _homeVar = TIGRE_HOME
    _pathVars = [TIGRE_HOME, TIGRE_CUDA_LIB]
    _processingField = [TOMO]
    _supportedVersions = [V3_1_3]

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(TIGRE_HOME, f'{TIGRE}-{TIGRE_DEFAULT_VERSION}',
                         description="Root folder where tigre was extracted.")
        cls._defineVar(TIGRE_ENV_ACTIVATION, TIGRE_DEFAULT_ACTIVATION_CMD)
        cls._defineVar(TIGRE_CUDA_LIB, pwem.Config.CUDA_LIB)
        
    @classmethod
    def getTigreEnvActivation(cls):
        return cls.getVar(TIGRE_ENV_ACTIVATION)

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch Tigre. """
        environ = pwutils.Environ(os.environ)
        if 'PYTHONPATH' in environ:
            # this is required for python virtual env to work
            del environ['PYTHONPATH']
        cudaLib = cls.getVar(TIGRE_CUDA_LIB, pwem.Config.CUDA_LIB)
        environ.addLibrary(cudaLib)



    @classmethod
    def defineBinaries(cls, env):
        TIGRE_INSTALLED = '%s_%s_installed' % (TIGRE, TIGRE_DEFAULT_VERSION)
        installationCmd = cls.getCondaActivationCmd()
        
        # Cloning the repos
        installationCmd += ' git clone https://github.com/Vilax/TIGRE_forked %s && ' % TIGRE_WRAPPER

        installationCmd += ' cd %s &&' % TIGRE_WRAPPER
        installationCmd += ' git checkout fpTigre && '  
        installationCmd += ' cd Python/cli/ &&'
        # Installing tigre in the environment
        installationCmd += ' conda env create -y -n %s -f tigreEnv.yml && ' % TIGRE_ENV_NAME
       
        # Activate new the environment
        installationCmd += 'conda activate %s && ' % TIGRE_ENV_NAME
        
        installationCmd += ' cd .. && ' 
        installationCmd += ' cd .. && '
        installationCmd += ' pip install . && '

        # Flag installation finished
        installationCmd += ' cd .. && touch %s' % TIGRE_INSTALLED

        tigre_commands = [(installationCmd, TIGRE_INSTALLED)]
        envPath = os.environ.get('PATH', "")  # keep path since conda likely in there
        installEnvVars = {'PATH': envPath} if envPath else None

        env.addPackage(TIGRE,
                       version=TIGRE_DEFAULT_VERSION,
                       tar='void.tgz',
                       commands=tigre_commands,
                       neededProgs=cls.getDependencies(),
                       vars=installEnvVars,
                       default=True)

    @classmethod
    def getDependencies(cls):
        # try to get CONDA activation command
        condaActivationCmd = cls.getCondaActivationCmd()
        neededProgs = []
        if not condaActivationCmd:
            neededProgs.append('conda')
        return neededProgs

    @classmethod
    def runTigre(cls, protocol, args, cwd=None):
        """ Run Tigre command from a given protocol. """
        cmd = cls.getCondaActivationCmd() + " "
        cmd += cls.getTigreEnvActivation()
        cmd += f" && CUDA_VISIBLE_DEVICES=%(GPU)s "
        tigreCmd = f'&& python3 {args} '
        protocol.runJob(cmd, tigreCmd, env=cls.getEnviron(), cwd=cwd)

    @classmethod
    def getTigreProgram(cls, tigreProgram):
        return join(cls.getVar(TIGRE_HOME), 'tigreWrapper', 'Python/cli', tigreProgram) + ' '
