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

# ALGORITHMS FAMILIES
familyHelps = dict()
familyHelps['family'] =   ' Select an option to reconstruct tomograms: \n ' \
                ' * _Exact algorithms_: WBP, FDK. \n ' \
                ' * _Gradient-based algorithms_: SART, SIRT, OS-SART, ASD-POCS, OS-ASD-POCS, PCSD, AW-PCSD, AW-ASD-POCS \n ' \
                ' * _Krylov algorithms_: CGLS, LSQR, LSMR, hybrid-LSQR, AB-GMRES, BA-GMRES, IRN-TV-CGLS  \n ' \
                ' * _Statistical reconstructions_: MLEM \n ' \
                ' * _Variational methods_: FISTA, SART-TV \n ' \
                ' See the help of each category once you select it. '

filterHelps = dict()
filterHelps['filter'] = 'Select the filter to be applied: \n ' \
                        '* **Ram-lak*: (default) Provides the absolute maximum resolution and sharpest possible edges. \n' \
                        ' It severely amplifies high-frequency noise, making the image look very grainy. \n' \
                        '* **Shepp-Logan*: Multiplies the pure ramp filter by a sinc function curve. \n' \
                        ' Mildly suppresses the highest frequencies. Reduces high-frequency noise while \n ' \
                        ' preserving most of the edge sharpness. It serves as an excellent standard "middle ground" \n' \
                        ' filter for general reconstruction \n ' \
                        '* **Cosine*: Multiplies the ramp filter by a cosine curve. \n' \
                        ' Starts rolling off the high frequencies earlier and more steeply than Shepp-Logan \n '\
                        ' Produces a smoother image with less noise, but structural edges will be noticeably softer \n' \
                        '* **Hamming*: Multiplies the ramp filter by a Hamming window (0.54 + 0.46 cos(pi f)). \n' \
                        ' A heavy smoothing filter that aggressively suppresses high frequencies (though it \n ' \
                        ' doesnt quite drop to zero at the absolute maximum frequency) \n ' \
                        ' Excellent for very noisy data or soft-tissue imaging where high contrast resolution \n' \
                        ' isnt critical \n ' \
                        '* **Hann*: Multiplies the ramp filter by a Hann window (0.5 + 0.5 cos(pi f)). \n ' \
                        ' Extremely similar to the Hamming filter, but its curve drops completely to zero at \n ' \
                        ' the maximum cutoff frequency.Trade-off: The most aggressive smoother on this list. \n ' \
                        ' It yields the lowest noise but the blurriest edges.'


varMethodHelp = dict()
varMethodHelp['varMethod'] = 'Select an option to reconstruct tomograms: \n ' \
                            '* **FISTA**: (default) Fast iterative shrinkage thresholding algorithm. \n' \
                            '* **SART-TV**: Simultaneous algebraic reconstruction technique with total variation. \n'

exactsMethodHelp = dict()
exactsMethodHelp['exactMethod'] = 'Select an option to reconstruct tomograms: \n ' \
                                '* **WBP**: (default) Weighted back projection or filtered back projection. \n' \
                                '* **FDK**: Feldkamp-Davis-Kress algorithm. \n'
gradientMethodHelp = dict()
gradientMethodHelp['gradientMethod'] = 'Select an option to reconstruct tomograms: \n ' \
                                        '* **SART**: (default) Simultaneous algebraic reconstruction technique. \n' \
                                        '* **SIRT**: Simultaneus iterative reconstruction technique. \n' \
                                        '* **OS-SART**: Ordered Subset Simultaneous Algebraic Reconstruction Technique. \n' \
                                        '* **ASD-POCS**: Adaptative Steepest Descent with Projections Onto the Convex Set. \n' \
                                        '* **OS-ASD-POCS**: Oriented Subsets version of ASD-POC. \n' \
                                        '* **PCSD**: Projection-controlled steepest descent method. \n' \
                                        '* **AW-PCSD**: Adaptive Weighted version of PCSD. \n' \
                                        '* **AW-ASD-POCS**: Adaptive Weighted TV (edge preserving) version of ASD-POCS. \n'
KrylovMethodHelp = dict()
KrylovMethodHelp['KrylovMethod'] =  'Select an option to reconstruct tomograms: \n ' \
                                    '* **CGLS**: (default) Conjugate Gradient Least Squares. \n' \
                                    '* **LSQR**: Least Squares QR factorization. \n' \
                                    '* **LSMR**: Least Squares MINRES. \n' \
                                    '* **hybrid LSQR_: Hybrid Least Squares QR factorization. \n' \
                                    '* **AB-GMRES**: AB-generalized minimal residual method. \n' \
                                    '* **BA-GMRES**: BA-generalized minimal residual method. \n' \
                                    '* **IRN-TV-CGLS**: Total variation CGLS. \n' 
statsMethodHelp = dict()
statsMethodHelp['statsMethod'] = 'Select an option to reconstruct tomograms: \n ' \
                                '* **MLEM**: (default) Maximum likelihood Expectation maximization. \n'

## GRADIENT ALGORITHMS
# SART
SARThelps = dict()
SARThelps['generalHelp'] = ''
SARThelps['SARTiter'] = ' The number of global iterations. One iteration is equivalent to processing ' \
                        ' and updating the 3D volume matrix using absolutely all available projection ' \
                        ' angles. It defines the global progress of the algorithm towards the optimal ' \
                        ' solution and the computational cost. It directs the trade-off between geometric ' \
                        ' sharpness and noise. Use an "early stopping" strategy. Increase it progressively ' \
                        ' until the object structure is well defined, and stop it before it begins to ' \
                        ' reconstruct and amplify the statistical background noise (semi-convergence phenomenon). ' \
                        ' * **Typical values**: from 15 to 50.'
SARThelps['SARTlambda'] =   ' The relaxation parameter or step size. It is a scalar that multiplies the error ' \
                            ' update matrix before adding it to the current volume. It controls the speed and ' \
                            ' stability of convergence. A high value makes abrupt corrections to reach the ' \
                            ' solution quickly; a low value takes conservative steps to avoid instability in ' \
                            ' problems with a lot of noise. If you notice that the reconstruction has blotches ' \
                            ' (artifacts), diverges, or oscillates between iterations, you should reduce it. ' \
                            ' If your data is clean, you can keep it high to save computation time. ' \
                            ' * **Typical values**: from 0.8 to 1.0.' 
SARThelps['SARTlambdared'] =    'A dynamic decay factor. Upon completing each iteration, the algorithm recalculates ' \
                                ' the relaxation by making lambda_new = lambda_current} * lambda_red. ' \
                                ' It allows the algorithm to take "large steps" in the early iterations (when the ' \
                                ' error is massive) and "tiny steps" in the final iterations. This refines subtle ' \
                                ' details and prevents the algorithm from bouncing or oscillating around the optimal ' \
                                ' value. Adjust it slightly below 1. If your reconstruction stops improving or ' \
                                ' "stalls" before recovering fine details, bring it closer to 1 so that the ' \
                                ' algorithms deceleration is slower. ' \
                                ' * **Typical values**: from 0.95 to 0.99.'



# SIRT
SIRThelps = dict()
SIRThelps['generalHelp'] = ''
SIRThelps['SIRTiter'] = ' The total number of global cycles in which the entire set of projections is processed ' \
                        ' simultaneously and the complete 3D volume is updated. It defines the limit of the ' \
                        ' iterative process, controlling the level of image convergence and execution time on the GPU. ' \
                        ' Since SIRT converges very smoothly and stably, you can adjust it flexibly. Increase it ' \
                        ' if you notice that the general image still lacks sufficient contrast, but keep in mind ' \
                        ' that an excessive number will not improve sharpness due to the inherently diffusive ' \
                        ' nature of the method. ' \
                        ' * **Typical values**: from 20 to 100.'
SIRThelps['SIRTlambda'] =   ' The initial relaxation parameter that weights the size of the update step applied to ' \
                            ' the volume in each global iteration. It regulates the aggressiveness with which the ' \
                            ' algorithm corrects projection errors. High values force faster convergence; low values ' \
                            ' ensure a stable and gradual evolution. Keep it within a safe, standard range. If you ' \
                            ' observe strange oscillations or artifacts in the reconstruction, lower it slightly; ' \
                            ' if the process is too slow to converge, bring it closer to 1.0.' \
                            ' * **Typical values**: from 0.5 to 1.0 (values 0.5 - 0.9 are sometimes used to maximize stability).'
SIRThelps['SIRTlambdared'] =    ' A multiplicative decay coefficient that dynamically reduces the relaxation parameter ' \
                                ' at the end of each iteration (λnew ​= λcurrent*λred​). It allows fine-tuning at the end ' \
                                ' of the reconstruction. It makes the algorithm apply large corrections at the beginning ' \
                                ' to quickly find the global structure and reduces changes in the final iterations to ' \
                                ' prevent oscillations or numerical instabilities. Set it close to 1 if you want the ' \
                                ' relaxation to decay very slowly, or adjust it a bit lower if the algorithm needs to ' \
                                ' slow down its rate of change more quickly.' \
                                ' *Typical values*: 0.95 to 0.99.'

# OS-SART
OSSARThelps = dict()
OSSARThelps['generalHelp'] = ''
OSSARThelps['OSSARTiter'] = ' The total number of global cycles in which all subsets into which the projections have ' \
                            ' been divided are traversed. It controls the total calculation time and the convergence ' \
                            ' level of the 3D volume towards the final solution. Due to the acceleration provided by ' \
                            ' the subsets, OSSART usually requires fewer global iterations than SART to achieve an ' \
                            ' equivalent level of sharpness. ' \
                            ' *Typical values*: from 10 to 30.'

OSSARThelps['OSSARTlambda'] = ' The initial relaxation parameter that determines the step size applied to the error '\
                            ' updates calculated by each subset. It modulates the intensity of the corrections. ' \
                            ' Since subsets contain less data than the complete set, the changes are more local ' \
                            ' and require careful control of stability. Keep it at a moderate or high value. ' \
                            ' If you notice excessive noise, circular artifacts, or strange oscillations due to ' \
                            ' jumps between subsets, reduce it slightly. ' \
                            ' *Typical values*: from 0.2 to 0.5.'
OSSARThelps['OSSARTlambdared'] = ' A multiplicative decay coefficient that gradually reduces the relaxation parameter ' \
                                ' at the end of each global iteration (lambda_new = lambda_current *lambda_red).'\
                                ' It ensures long-term stability. It allows the algorithm to progress quickly at the ' \
                                ' beginning and smoothly slow down changes as it approaches the solution, preventing ' \
                                ' it from oscillating indefinitely among the different subsets. Bring it closer to 1.0 ' \
                                ' if you want the relaxation factor to decay very slowly, or lower it slightly if you ' \
                                ' notice that the reconstruction oscillates or fails to stabilize in the final ' \
                                ' iterations.' \
                                ' *Typical values*: from 0.95 to 0.99.'

# ASD-POCS
ASDPOCShelps = dict()
ASDPOCShelps['generalHelp'] = ''
ASDPOCShelps['ASDPOCSiter'] =   ' The total number of cycles that jointly execute the algebraic projection phase ' \
                                ' and the Total Variation minimization phase. It controlls how many times the ' \
                                ' volume is completely refined. Since ASD-POCS converges efficiently thanks to ' \
                                ' TV regularization, a moderate value is usually enough to achieve an excellent ' \
                                ' result without over-executing. ' \
                                ' *Typical values*: from 20 to 50.'
ASDPOCShelps['ASDPOCStviter'] = ' The number of inner steepest descent gradient iterations applied exclusively to ' \
                                ' minimize the Total Variation norm in each global cycle. It controls how exhaustive ' \
                                ' the noise suppression and texture smoothing (piecewise-constant effect) are at ' \
                                ' each step before re-verifying the model with the actual data. If the resulting ' \
                                ' image shows too much residual noise, increase this value moderately. If the image ' \
                                ' starts to look like a painting or loses fine details, reduce it. '\
                                ' *Typical values*: from 10 to 20.'
ASDPOCShelps['ASDPOCSalpha'] =  ' The scalar step size used in gradient descent for Total Variation minimization. ' \
                                ' It regulates the magnitude of the correction with which Total Variation is reduced ' \
                                ' in each internal tviter step. Adjust it with caution. If it is too large, the ' \
                                ' TV minimization will diverge or generate geometric artifacts; if it is too small, ' \
                                ' the noise filtering will be practically imperceptible.' \
                                ' *Typical values*: from 0.001 to 0.02.'
ASDPOCShelps['ASDPOCSlambda'] = ' The initial relaxation parameter or step size that weights the algebraic update ' \
                                ' derived from the projection phase. It modulates the aggressiveness of the corrections ' \
                                ' based on geometric backprojection to ensure that the solution respects the consistency ' \
                                ' of the measured data. Keep it within a safe, standard range. If you notice noticeable ' \
                                ' oscillations or instability in the first global iterations, reduce it slightly.' \
                                ' *Typical values*: from 0.5 to 1.0.'
ASDPOCShelps['ASDPOCSlambdared'] =  ' A multiplicative factor that dynamically decreases the relaxation parameter upon ' \
                                    ' completing each global iteration (lambda_new = lambda_current*lambda_red$). ' \
                                    ' It guarantees long-term stability, allowing the algorithm to make large corrections ' \
                                    ' at the beginning and smoothly slow down toward the end to prevent oscillations. ' \
                                    ' Bring it close to 1.0 so that the decay of the relaxation is slow and controlled. ' \
                                    ' *Typical values*: from 0.95 to 0.99.'
ASDPOCShelps['ASDPOCSalphared'] =   ' The alpha_red parameter acts exclusively on the Total Variation phase, controlling ' \
                                    ' how the filtering intensity evolves over time.' \
                                    ' A multiplicative decay factor that gradually reduces the gradient step size of ' \
                                    ' the Total Variation (alpha) at the end of each global cycle ' \
                                    ' (alpha_new = alpha_current * alpha_red). t guarantees convergence and fine stability.' \
                                    ' It allows the algorithm to apply very aggressive noise reduction in the early ' \
                                    ' iterations (when the initial image is chaotic) and reduces the smoothing intensity ' \
                                    ' towards the end to preserve subtle textures and details, preventing oscillations ' \
                                    ' in the solution. Bring it closer to 1.0 if you want the TV filter to maintain its ' \
                                    ' strength throughout most of the process. If you notice that in the final iterations ' \
                                    ' the image gets too "washed out" (watercolor effect) or unstable artifacts appear, ' \
                                    ' lower it slightly (e.g., 0.95) so the algorithm reduces the smoothing step much faster' \
                                    ' *Typical values*: from 0.95 to 0.99.'
ASDPOCShelps['ASDPOCSratio'] =  ' The maximum relative change ratio allowed between successive iterations of the Total ' \
                                ' Variation norm. It works as an adaptive control mechanism or internal stopping criterion ' \
                                ' to evaluate the progress of TV minimization within the loop. It is rarely necessary ' \
                                ' to modify it manually. ' \
                                ' *Typical values*: 0.95 to 0.99.'


# OS-ASD-POCS
OSASDPOCShelps = dict()
OSASDPOCShelps['generalHelp'] = ''
OSASDPOCShelps['OSASDPOCSiter'] =   ' The total number of global cycles in which all projection subsets are traversed and ' \
                                    ' Total Variation minimization is applied. It controls the progress of the global ' \
                                    ' convergence of the 3D volume. Thanks to subset acceleration, you can use lower ' \
                                    ' values than in standard ASD-POCS to achieve equivalent results. ' \
                                    ' *Typical values*: from 10 to 25.'
OSASDPOCShelps['OSASDPOCStviter'] = ' The number of internal gradient descent iterations executed exclusively to minimize ' \
                                    ' the Total Variation norm in each global cycle. It defines how exhaustive the noise ' \
                                    ' suppression and region smoothing (piecewise-constant effect) are at each step. ' \
                                    ' Increase this value if you notice too much residual noise in the background; ' \
                                    ' reduce it if the image starts losing fine details of the object. ' \
                                    ' *Typical values*: from 10 to 20.'
OSASDPOCShelps['OSASDPOCSalpha'] =  ' The scalar step size used in gradient descent for Total Variation minimization. ' \
                                    ' It regulates the intensity of the correction when applying TV smoothing to the ' \
                                    ' volume. Excessively high values generate geometric distortions or artifacts, ' \
                                    ' while very low values make the TV filter ineffective. ' \
                                    ' *Typical values*: 0.001 to 0.02.'
OSASDPOCShelps['OSASDPOCSlambda'] = ' The initial relaxation parameter that weights the subset-based algebraic updates. ' \
                                    ' It modulates the stability of the process against the rapid local changes introduced ' \
                                    ' by ordered subsets. You should keep it moderate or low to avoid oscillations or ' \
                                    ' numerical instabilities caused by jumps between subsets. ' \
                                    ' *Typical values*: from 0.2 to 0.5.'

OSASDPOCShelps['OSASDPOCSlambdared'] =  ' A multiplicative decay factor that gradually reduces the relaxation parameter at ' \
                                        ' the end of each global iteration (lambda_new = lambda_current*lambda_red). ' \
                                        ' It ensures long-term stability, allowing large corrections at the beginning and ' \
                                        ' a smooth slowdown at the end. Bring it closer to 1.0 so that the relaxation decay ' \
                                        ' is slow and controlled. ' \
                                        ' *Typical values*: from 0.95 to 0.99'
OSASDPOCShelps['OSASDPOCSalphared'] = ASDPOCShelps['ASDPOCSalphared']
OSASDPOCShelps['OSASDPOCSratio'] =  ' The maximum relative change ratio allowed between successive iterations of the Total ' \
                                    ' Variation operator. It acts as an internal adaptive stopping criterion to regulate ' \
                                    ' the convergence of the TV phase. It rarely requires manual modification. ' \
                                    ' *Typical values*: from 0.95 to 0.99.'

# PCSD
PCSDhelps = dict()
PCSDhelps['generalHelp'] = ''
PCSDhelps['PCSDiter'] = ' The total number of global cycles executed by the algorithm, where each cycle alternates ' \
                        ' between updating algebraic consistency and minimizing the Total Variation (TV) term. ' \
                        ' It defines the limit of the global iterative process, controlling the degree of convergence ' \
                        ' of the reconstructed image. As it is a highly efficient algorithm for solving problems ' \
                        ' with TV constraints, you can increase it moderately until the image achieves the desired ' \
                        ' visual sharpness. ' \
                        ' *Typical values*: from 20 to 50.'

# AW-PCSD
AWPCSDhelps = dict()
AWPCSDhelps['AWPCSDiter'] = ' The total number of global primal-dual optimization cycles executed by the algorithm ' \
                            ' to converge towards the optimal solution. ' \
                            ' *Typical values*: from 20 to 50.'


# AS-ASD-POCS
AWASDPOCShelps = dict()
AWASDPOCShelps['AWASDPOCSiter'] =  ' The total number of global cycles in which the algebraic projection phase and ' \
                                    ' the adaptive Total Variation minimization phase alternate. By incorporating ' \
                                    ' adaptive weights, a moderate number is usually enough to stabilize the image. ' \
                                    ' *Typical values*: from 20 to 50.'
AWASDPOCShelps['AWASDPOCStviter'] =    ' The number of internal gradient descent iterations executed to solve the ' \
                                        ' minimization of the weighted Total Variation norm in each global cycle. ' \
                                        ' It controls how exhaustive the noise suppression and adaptive weight ' \
                                        ' adjustment process is at each iterative step. Increase it slightly if the ' \
                                        ' background noise persists; reduce it if you notice the image losing fine textures.' \
                                        ' *Typical values*: from 10 to 20.'

AWASDPOCShelps['AWASDPOCSalpha'] = ' The scalar step size used in gradient descent for adaptive Total Variation ' \
                                    ' minimization. It regulates the intensity of the correction when applying ' \
                                    ' weighted smoothing over the volume. Very high values generate distortions, ' \
                                    ' while very low values nullify the effect of the adaptive filter. ' \
                                    ' *Typical values*: 0.001 to 0.02. '
AWASDPOCShelps['AWASDPOCSlambda'] =    ' The initial relaxation parameter that weights the algebraic update derived ' \
                                        ' from the projection phase. It modulates the speed and stability of corrections ' \
                                        ' based on geometric backprojection. Keep it in a moderate or safe range. ' \
                                        ' If you observe strange oscillations in the early iterations, reduce it slightly. ' \
                                        ' *Typical values*: from 0.5 to 1.0.'
AWASDPOCShelps['AWASDPOCSlambdared'] = ' A multiplicative decay factor that dynamically reduces the relaxation parameter ' \
                                        ' upon completing each global iteration (lambda_new = lambda_current*lambda_red). ' \
                                        ' It guarantees long-term stability, allowing large corrections at the beginning ' \
                                        ' and a smooth slowdown at the end of the process. Bring it closer to 1.0 so that ' \
                                        ' the relaxation decay is slow and controlled. ' \
                                        ' *Typical values*: from 0.95 to 0.99.'
AWASDPOCShelps['AWASDPOCSalphared'] = ASDPOCShelps['ASDPOCSalphared']
AWASDPOCShelps['AWASDPOCSratio'] = ' The maximum relative change ratio allowed between successive iterations of the Total ' \
                                    ' Variation operator. It acts as an internal adaptive stopping criterion to regulate ' \
                                    ' the convergence of the TV minimization phase. In the vast majority of practical cases, ' \
                                    ' the internal default value of the CERN TIGRE repository offers optimal performance, so ' \
                                    ' it rarely requires manual modification. ' \
                                    ' *Typical values*: from 0.95 to 0.99.'


## KRYLOV ALGORITHMS
# CGLS
CGLShelps = dict()
CGLShelps['generalHelp'] = ''
CGLShelps['CGLSiter'] = ' It defines the number of conjugate directions that the algorithm will calculate. ' \
                        ' At iteration k, the algorithm expands the mathematical subspace (Krylov subspace) ' \
                        ' to k dimensions and finds the strict mathematical optimal image within that space.' \
                        ' In CGLS, the number of iterations acts as your only regularizer. This is known in ' \
                        ' inverse mathematics as "regularization by early stopping".' \
                        ' The semi-convergence phenomenon: This is the critical behavior of CGLS. In the first ' \
                        ' few iterations, the algorithm reconstructs the main and global features of the object ' \
                        ' (low frequencies). As you progress, it begins to recover fine details and edges. ' \
                        ' However, if you perform too many iterations, CGLS will attempt to fit the image to the ' \
                        ' noise of the original data, causing the image to "blow up" and become filled with ' \
                        ' high-frequency noise.How to adjust it: You must find the exact point right before the ' \
                        ' noise begins to dominate. If the image looks smeared or blurry: increase the iterations.' \
                        ' If the image has excessive grainy noise: reduce the iterations immediately. ' \
                        ' *Typical values*: Values between 10 and 30 iterations are the most common.'

# LSQR
LSQRhelps = dict()
LSQRhelps['generalHelp'] = ''
LSQRhelps['LSQRiter'] = ' It defines the number of steps in the bidiagonalization process. In each iteration, ' \
                        ' the algorithm expands the mathematical subspace (the search space) by one additional ' \
                        ' dimension and extracts the best possible image from that new space.' \
                        ' The number of iterations is your only regularizer through the principle of "early stopping" ' \
                        ' Lacking a built-in noise filter, the algorithm relies on you stopping it before it starts  ' \
                        ' reconstructing the noise. ' \
                        ' The first few iterations reconstruct the large shapes (low frequencies). Middle iterations ' \
                        ' outline the edges (high frequencies). If you perform too many iterations, LSQR will attempt ' \
                        ' to make the image fit perfectly with the detector noise, ruining the reconstruction ' \
                        ' with an intense grainy pattern. ' \
                        ' The golden rule is if the tomogram is blurry, increase the iterations; and if it looks grainy, ' \
                        ' decrease the iterations. ' \
                        ' *Typical values*: Values between 10 and 50 iterations.'

# LSMR
LSMRhelps = dict()
LSMRhelps['generalHelp'] = ''
LSMRhelps['LSMRiter'] = ' It defines the number of Golub-Kahan bidiagonalization steps to expand the search subspace. ' \
                        ' It determines how much detail you extract from the data. LSMR is very tolerant to an excess ' \
                        ' of iterations. If you set a good lambda value, you can afford to let the algorithm run for ' \
                        ' more iterations (e.g., 30 to 60) to recover fine details, without the constant fear that ' \
                        ' the image will blow up in noise. ' \
                        ' *Typical Values*: Between 10 to 50'
LSMRhelps['LSMRlambda'] =   ' It is the weight lambda in the previous equation. It penalizes the L2 norm of the ' \
                            ' image (||x||^2), which mathematically means it punishes pixels that acquire absurdly ' \
                            ' high numerical values. In tomography problems, noise usually manifests as pixels with ' \
                            ' extreme positive and negative values that cancel each other out.' \
                            ' It acts as a damper (damping). Tikhonov regularization simply assumes that "the total ' \
                            ' energy of the image must not be infinite". It serves to stabilize the reconstruction ' \
                            ' when you have very noisy data or very few projection angles.' \
                            ' *Typical Values*: In the range 0.001 to 1.0' \
                            ' If you set it to 0: LSMR behaves almost identically to LSQR, ' \
                            ' depending purely on stopping the iterations in time. If you have noise: Start ' \
                            ' increasing it slowly (e.g., 0.1, 0.5, 1.0). You will notice that the noise is ' \
                            ' attenuated without creating the "oil painting" effect that Total Variation sometimes ' \
                            ' causes. If you go too far  (very high values): The algorithm will penalize the pixel ' \
                            ' values so much that the entire image will lose contrast, look dark or "washed out", ' \
                            ' and the edges will become excessively blurred.'

# hybrid-LSQR
hybridLSQRhelps = dict()
hybridLSQRhelps['generalHelp'] = ''
hybridLSQRhelps['hybridLSQRiter'] = ' The higher the number of iterations, the closer the image gets to the exact ' \
                                    ' mathematical solution, improving the geometric resolution of the edges and better ' \
                                    ' recovering subtle contrasts. ' \
                                    ' It represents the number of steps in the algorithm, which mathematically defines ' \
                                    ' the dimension k of the Krylov subspace. In each iteration, the  algorithm performs ' \
                                    ' a multiplication by matrix A (simulated projection) and by its transpose At ' \
                                    ' (backprojection). With each step, it mathematically expands the search space by ' \
                                    ' adding a new orthogonal basis vector and recalculates the volume x that best fits ' \
                                    ' the data under the imposed regularization. You should use it as your stopping ' \
                                    ' criterion (when to halt the reconstruction so as not to waste computational ' \
                                    ' time in vain). In a standard LSQR, you would have to stop it strategically before the ' \
                                    ' image becomes filled with noise. However, since you are using the Hybrid version, the ' \
                                    ' algorithm already controls the noise for you. This means you can play with this ' \
                                    ' parameter without fear: increase it freely and check the resulting image. When you notice ' \
                                    ' that raising the iterations value no longer produces a visually perceptible improvement ' \
                                    ' in the sharpness of the reconstruction, that is the value you should leave fixed in your ' \
                                    ' script.' \
                                    ' *Typical values*: from 30 to 60.'

# AB-GMRES
ABGMREShelps = dict()
ABGMREShelps['ABGMRESiter'] =   ' It defines the Krylov subspace dimension k and the number of Arnoldi orthogonalization steps. ' \
                                ' Each iteration projects data (A · A^t), generates a new orthogonal basis vector, and updates ' \
                                ' a Hessenberg matrix to find the 3D volume with the minimum error in that k-dimensional space. ' \
                                ' It acts as your core regularization mechanism, balancing detail and noise. Iterations act as a ' \
                                ' frequency filter, early steps recover global shapes (low frequencies), while later steps ' \
                                ' resolve fine details and sensor noise (high frequencies). ' \
                                ' *Typical values*: from 5 to 50 iterations. Values greater than  50 are rarely used.'

# BA-GMRES
BAGMREShelps = dict()
BAGMREShelps['BAGMRESiter'] =  ABGMREShelps['ABGMRESiter']

# IRN-TV-CGLS
IRNTVCGLShelps = dict()
IRNTVCGLShelps['IRNTVCGLSiter'] =   ' The number of steps of the Conjugate Gradient Least Squares (CGLS) method that are ' \
                                    ' executed to solve the weighted least squares subproblem within each cycle. ' \
                                    ' It defines how much computational precision is invested in solving the matrix before ' \
                                    ' the outer loop updates the Total Variation weights. ' \
                                    ' Keep it relatively low. You do not need the CGLS to converge perfectly at each step, ' \
                                    ' as the system of equations will change slightly in the next outer iteration.' \
                                    ' *Typical Values*: from 5 to 20.'
IRNTVCGLShelps['IRNTVCGLStvlambda'] =   ' The regularization hyperparameter that weights the penalty term (Total Variation) ' \
                                        ' against the data fidelity term (the geometric projection). It controls the ' \
                                        ' aggressiveness of the filter. A high value forces a very flat image with very ' \
                                        ' hard edges (drastically reducing noise). A low value relies more on the measured ' \
                                        ' data (preserving fine textures but leaving the noise intact). You have to adjust it ' \
                                        ' empirically through trial and error for your specific scanner or sample. You must ' \
                                        ' look for the exact balance point: if you raise it too much, your image will look like ' \
                                        ' a watercolor painting with no inner texture ("staircasing effect"); if you lower it ' \
                                        ' too much, the algorithm will be equivalent to a standard CGLS full of noise.' \
                                        ' *Typical values*:  in the range from  0.1, 1, 10, 50.'
IRNTVCGLShelps['IRNTVCGLSinnerouter'] = ' The number of times the weights of the regularization matrix are recalculated and ' \
                                        ' applied (the reweighting step of the IRN norm). It is the main driver of convergence ' \
                                        ' towards the Total Variation solution. Each outer cycle refines which parts of the ' \
                                        ' image are noise (and should be flattened) and which parts are structural edges (and ' \
                                        ' should be preserved). Increase it progressively if you see that the algorithm has not ' \
                                        ' yet finished cleaning up the background noise, but without overdoing it to avoid ' \
                                        ' multiplying reconstruction times, since each outer iteration multiplies the cost ' \
                                        ' of the inner iterations.'


# STATISTICAL ALGORITHMS
# MLEM
MLEMhelps = dict()
MLEMhelps['generalHelp'] = ''
MLEMhelps['MLEMiter'] =  ''


## VARIATIONAL ALGORITHMS
# FISTA
FISTAhelps = dict()
FISTAhelps['generalHelp'] = ''
FISTAhelps['FISTAiter'] =   ' This parameter defines the number of times the main loop of the algorithm will be executed. ' \
                            ' In each iteration, the algorithm takes a step in the direction of the gradient to reduce the ' \
                            ' projection error and then applies regularization. The more iterations the higher probability ' \
                            ' to to find the minimum of the objective function. ' \
                            ' *Typical values*: It is usually between 50 and 200 ' \
                            ' If the final image keeps changing significantly or looks blurry, more iterations will be needed ' \
                            ' If the image no longer improves between iteration 50 and 60, setting it to 100 is a waste of ' \
                            ' computation time. '
FISTAhelps['FISTAtvlambda'] =   ' It represents the scalar tvlambda in objective function. It is the regularization hyperparameter. ' \
                                ' It determines how much importance is given to minimizing noise (Total Variation - TV) versus ' \
                                ' respecting the pure measured data. It controls the level of image smoothing. ' \
                                ' The TV penalty is excellent at removing noise and artifacts while preserving sharp edges. ' \
                                ' If your image has a lot of noise or streaking artifacts: Increase the value of tvlambda. ' \
                                ' If your image looks too artificial, like an "oil painting," with a blocky effect, or loses ' \
                                ' very fine textures, then decrease tvlambda. ' \
                                ' *Typical values*: It often requires a sweep testing with 1, 10, 50.'
FISTAhelps['FISTAtviter']  =    ' In each main iteration the algorithm has to evaluate the TV term. However, for 3D TV, there is ' \
                                ' no direct mathematical solution. Therefore, TIGRE has to run an internal iterative mini-algorithm ' \
                                ' to solve that part of the equation. tviter is the number of iterations of that internal sub-loop. ' \
                                ' It ensures that the mathematical regularization dictated by tvlambda is calculated accurately ' \
                                ' at each step of the main loop. If you notice that increasing tvlambda is not smoothing the ' \
                                ' image as it should, then the internal loop might not be converging then try to  increase tviter. ' \
                                ' *Caution*: This parameter has a massive impact on computation time. It multiplies the algorithms ' \
                                ' workload (total iterations = iterations * tviter). ' \
                                ' *Typical values*: Values between 20 and 100 are the standard in the community. A good rule of thumb '\
                                ' is to keep a moderate value, 50, and try not to touch it unless issues with the tvlambda effect.'

# SART-TV
SARTTVhelps = dict()
SARTTVhelps['generalHelp'] = ''
SARTTVhelps['SARTTViter'] = ' It is the number of times the entire algorithm alternates between the algebraic reconstruction step (SART) ' \
                            ' and the filtering step (TV). It defines the overall convergence of the algorithm toward the final tomogram. ' \
                            ' SART usually converges quite fast in its early steps compared to gradient-based methods. ' \
                            ' *Typical values*: Values between 15 and 50 are common. Increasing it much further does not always improve ' \
                            ' the tomogram significantly and consumes a lot of time, since SART is computationally heavy.'
SARTTVhelps['SARTTVtvlambda'] = ' tvlambda works as the step size in the gradient descent of the TV. It controls how aggressive the noise ' \
                                ' correction applied to the tomogram is after each SART step. It dictates how much you change the pixel ' \
                                ' values to flatten out homogeneous areas. ' \
                                ' If the image is noisy, you need to increase tvlambda. However, since it acts as a step size here, ' \
                                ' if you set it too high, the TV step can become mathematically unstable (diverge), creating strange ' \
                                ' artifacts or failing completely.' \
                                ' *Typical values*: It often requires values around 10 to 50.'
SARTTVhelps['SARTTVtviter'] =   ' Minimizing Total Variation is not done in a single mathematical calculation.' \
                                ' The operation x - tvlambda · nabla TV(x) is repeated in an inner loop. tviter determines how many times ' \
                                ' this TV gradient descent is repeated before returning the image to SART for the next main iteration. ' \
                                ' It ensures that noise is properly removed in each global iteration before re-injecting new geometric ' \
                                ' information. This parameter is the main performance bottleneck in SART-TV. If you notice that the image ' \
                                ' is still noisy but increasing tvlambda destabilizes the reconstruction, the solution is to keep ' \
                                ' tvlambda moderate and increase tviter.' \
                                ' *Typical values*: Between 20 and 50. If you have few projections (a very incomplete scan), you may need ' \
                                ' a higher tviter (e.g., 100) to force the image to be "piecewise flat" and compensate for the lack of views.'