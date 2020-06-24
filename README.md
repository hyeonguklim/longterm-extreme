# Long-term extreme response prediction by MCS and PCE
## Schematic
![](https://github.com/hyeonguklim/longterm_extreme/blob/master/scheme/scheme.png)  
**Fig**: A simple SDOF structural system subject to wave loading and a quantity of interest (QoI)

## Authors/Collaborators
[HyeongUk Lim](https://hyeonguk.wordpress.com), [Lance Manuel](https://lancemanuel.netlify.com), and [Ying Min Low](https://www.eng.nus.edu.sg/cee/staff/low-ying-min/)

## Description
This study focuses on the development of efficient surrogate models by polynomial chaos expansion (PCE) for prediction of the long-term extreme surge motion of a moored floating offshore structure. The structure is subjected to first-order and second-order (difference-frequency) wave loading. Uncertainty in the long-term response arises from contrasting sea state conditions, characterized by significant wave height, Hs, and spectral peak period, Tp, and their relative likelihood of occurrence; these two variables are explicitly included in the PCE-based uncertainty quantification (UQ). In a given sea state, however, response simulations must be run for the associated Hs and Tp; in such simulations, typically, a set of random amplitudes and phases define an irregular wave train consistent with that sea state.  These random amplitudes and phases for all the frequency components in the wave train introduce additional uncertainty in the simulated waves and in the response.  The UQ framework treats these two sources of uncertainty---from Hs and Tp on the one hand, and the amplitude and phase vectors on the other---in a manner that efficiently yields long-term surge motion extreme predictions consistent with more expensive Monte Carlo simulations (MCS) that serve as the "truth" system. To reduce uncertainty in response extremes that result from sea states with a low likelihood of occurrence, importance sampling is employed with both MCS- and PCE-based extreme response predictions. Satisfactory performance with such efficient surrogate models can help in assessing the long-term response of various offshore structures.  

## Related Publications/Presentations
- Lim, H., Manuel, L., and Low, Y. M. "On Efficient Surrogate Model Development for Prediction of the Long-Term Extreme Response of a Moored Floating Structure." ASME. J. Offshore Mech. Arct. Eng. [[presentation]](https://github.com/hyeonguklim/longterm_extreme/blob/master/presentation/OMAE_2018-78763.pdf) [[paper]](https://doi.org/10.1115/1.4047545)

## MATLAB Codes
The developed MATLAB codes are explained in brief; you can also find line-by-line explanations inside the codes.

### prerequisite
You may need the following MATLAB toolboxes:
- Statistics and Machine Learning Toolbox
- Symbolic Math Toolbox

You may type the following command to see the toolboxes required for the analysis:
```Matlab
[fList,pList] = matlab.codetools.requiredFilesAndProducts('PCE.m')
disp({pList.Name}')
```

### data files

- You may need the following structural/hydrodynamic parameters for the analysis:
```Matlab
'M: mass (kg)'
'C: damping (N/m/s)'
'K: stiffness (N/m)'
'zeta: damping ratio to critical'
'N: number of harnomics'
'dw: width of a frequency bin'
'TF: transfer function (N/m)'
'wmin: minimum frequency (rad/sec)'
'wmax: maximum frequency (rad/sec)'
'w: frequency range'
'w_LF: frequency range for low frequency response'
'N_LF: number of harmonics for low frequency response'
'Diag_surge: diagonal terms in QTF'
```

### main functions

- `MCS_30min.m`: Monte Carlo simulation (MCS) of 30-min extreme surge response
- `IS_30min.m`: MCS by importance sampling (IS)
- `PCE.m`: surrogate model building using polynomial chaos expansion (PCE)

### subfunctions
The subfuctions needed for running the main functions are:

- `incdfHs`: inverse cumulative distribution function (cdf) of Hs
- `pdfHs`: probability density function (pdf) of Hs
- `pdfHsIS`: importance sampling density function of Hs
- `incdfTp`: inverse cdf of Tp
- `surge_max_mex`: extreme surge response calculation by [IFFT](https://www.mathworks.com/help/matlab/ref/ifft.html)

	#### usage for MCS:
	```Matlab
	surge_max_mex(q1,q2,0) % q1 and q2 are wave height and peak period in standard normal variable space
	```
	#### usage for IS:
	```Matlab
	surge_max_mex(h,t,1) % h and t are wave height and peak period in physical variable space
	```
- `Hermite_PC`: Hermite orthogonal polynomial family
- `multi_index`: multi-index for multivariate orthogonal functions

### post-processing
- `plot_exprob_mcs.m`: exceedance probability plots for MCS and IS results
- `plot_exprob_pce.m`: exceedance probability plots for MCS and PCE results

## Example
1. Clone this repository to your directory
2. Open `MCS_30min.m` in MATLAB and run as it is or you can change variables (ex: number of samples, `N_T`)
3. You can also run the others: `IS_30min.m` and `PCE.m`
4. The results will be saved in `.\mat_file\`
5. You can create exceedance probability plots using the saved .mat files by `plot_exprob_mcs.m` and `plot_exprob_pce.m`

## Figures
<img src="https://github.com/hyeonguklim/longterm_extreme/blob/master/figures/exprob_mcs.png" width="350" height="">

**Fig**: Exceedance probability estimation by MCS and importance sampling

Variance in exceedance probability estimation is reduced with the help of importance sampling.

<img src="https://github.com/hyeonguklim/longterm_extreme/blob/master/figures/exprob_pce.png" width="350" height="">

**Fig**: Exceedance probability estimation by MCS and PCE

## Contact
For any questions or comments, please email me at: hyeonguklim@gmail.com.

