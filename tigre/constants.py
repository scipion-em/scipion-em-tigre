# **************************************************************************
# *
# * Authors:     you (you@yourinstitution.email)
# *
# * your institution
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
TIGRE = 'tigre'
TIGRE_WRAPPER = 'tigreWrapper'
TIGRE_HOME = "TIGRE_HOME"

# Suported versions
V3_1_3 = '3.1.3'
TIGRE_DEFAULT_VERSION =  V3_1_3

TIGRE_ENV_NAME = '%s-%s' % (TIGRE, TIGRE_DEFAULT_VERSION)
TIGRE_ENV_ACTIVATION = 'TIGRE_ENV_ACTIVATION'
TIGRE_DEFAULT_ACTIVATION_CMD = 'conda activate %s' % TIGRE_ENV_NAME
TIGRE_CUDA_LIB =  'TIGRE_CUDA_LIB'

TIGRE_PLUGIN_VERSION = '1.0'