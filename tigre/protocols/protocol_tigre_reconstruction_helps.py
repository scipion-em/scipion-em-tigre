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


ALGORITHMS_EXACT = ['WBP', 'FDK']

ALGORITHMS_GRADIENT = ['SART', 'SIRT', 'OS-SART', 'ASD-POCS', 'OS-ASD-POCS', 'PCSD', 'AW-PCSD', 'AW-ASD-POCS']

ALGORITHMS_KRYLOV = ['CGLS', 'LSQR', 'LSMR', 'hybrid-LSQR', 'AB-GMRES', 'BA-GMRES', 'IRN-TV-CGLS']

ALGORITHMS_VARIATIONAL = ['FISTA', 'SART-TV']



## EXACT ALGORITHMS
# WBP
WBPhelps = dict()


# FDK
FDKhelps = dict()


## GRADIENT ALGORITHMS
# SART
SARThelps = dict()

# SIRT
SIRThelps = dict()

# OS-SART
OSSARThelps = dict()

# ASD-POCS
ASDPOCShelps = dict()

# OS-ASD-POCS
OSASDPOCShelps = dict()

# PCSD
PCSDhelps = dict()

# AW-PCSD
AWPCSDhelps = dict()

# AS-ASD-POCS
AWASDPOCShelps = dict()


## KRYLOV ALGORITHMS
# CGLS
CGLShelps = dict()

# LSQR
LSQRhelps = dict()

# LSMR
LSMRhelps = dict()

# hybrid-LSQR
hybridLSQRhelps = dict()

# AB-GMRES
ABGMREShelps = dict()

# BA-GMRES
BAGMREShelps = dict()

# IRN-TV-CGLS
IRNTVCGLShelps = dict()


## VARIATIONAL ALGORITHMS
# FISTA
FISTAhelps = dict()
FISTAhelps['generalHelp'] = ''
FISTAhelps['FISTAiter'] = ''
FISTAhelps['FISTAtvlambda'] = ''
FISTAhelps['FISTAtviter'] = ''

# SART-TV
SARTTVhelps = dict()
SARTTVhelps['SARTTViter'] = ''
SARTTVhelps['SARTTVtvlambda'] = ''
SARTTVhelps['SARTTVtviter'] = ''
SARTTVhelps['SARTTValphared'] = ''



