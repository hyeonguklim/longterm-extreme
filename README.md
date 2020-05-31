# Long-term extreme response prediction by MCS and PCE
## Schematic
![](https://github.com/hyeonguklim/longterm_extreme/blob/master/scheme/scheme.png)  
**Fig**: A simple SDOF structural system subject to wave loading and a quantity of interest (QoI)

## Authors/Collaborators
[HyeongUk Lim](https://hyeonguk.wordpress.com), [Lance Manuel](https://lancemanuel.netlify.com), and [Ying Min Low](https://www.eng.nus.edu.sg/cee/staff/low-ying-min/)

## Description
This study focuses on the development of efficient surrogate models by polynomial chaos expansion (PCE) for prediction of the long-term extreme surge motion of a moored floating offshore structure. The structure is subjected to first-order and second-order (difference-frequency) wave loading. Uncertainty in the long-term response arises from contrasting sea state conditions, characterized by significant wave height, Hs, and spectral peak period, Tp, and their relative likelihood of occurrence; these two variables are explicitly included in the PCE-based uncertainty quantification (UQ). In a given sea state, however, response simulations must be run for the associated Hs and Tp; in such simulations, typically, a set of random amplitudes and phases define an irregular wave train consistent with that sea state.  These random amplitudes and phases for all the frequency components in the wave train introduce additional uncertainty in the simulated waves and in the response.  The UQ framework treats these two sources of uncertainty---from Hs and Tp on the one hand, and the amplitude and phase vectors on the other---in a manner that efficiently yields long-term surge motion extreme predictions consistent with more expensive Monte Carlo simulations (MCS) that serve as the "truth" system. To reduce uncertainty in response extremes that result from sea states with a low likelihood of occurrence, importance sampling is employed with both MCS- and PCE-based extreme response predictions. Satisfactory performance with such efficient surrogate models can help in assessing the long-term response of various offshore structures.  

## Related Publications/Presentations
- HyeongUk Lim, Lance Manuel and Ying Min Low (2018), On Efficient Long-term Extreme Response Estimation for a Moored Floating Structure, Madrid, Spain, June 17-22, Proceedings of the ASME 2018 37th International Conference on Ocean, Offshore and Arctic Engineering, OMAE2018-78763. [[presentation]](https://github.com/hyeonguklim/longterm_extreme/blob/master/presentation/OMAE_2018-78763.pdf) [[paper]](https://doi.org/10.1115/OMAE2018-78763)

## MATLAB Codes
Short explanation of the developed MATLAB codes. You can also find line-by-line explanations inside the codes.

### data files

- `data.mat`: structural/hydrodynamic parameters and initial setting for the analysis; please see `input_desc` in the .mat file for description. You may type the following lines to see the description:

```Matlab
load data.mat input_desc
disp(input_desc)
```

### main functions

- `MCS_30min.m`: Monte Carlo simulation of 30-min extreme surge response
- `IS_30min.m`: MCS by importance sampling
- `pce.m`: surrogate model building using polynomial chaos expansion (PCE)

### subfunctions
The subfuctions needed for running the main functions are:

- `incdfHs`: inverse cumulative distribution function (cdf) of Hs
- `pdfHs`: probability density function (pdf) of Hs
- `pdfHsIS`: importance sampling density function of Hs
- `incdfTp`: inverse cdf of Tp
- `Jonswap`: [JONSWAP](https://wikiwaves.org/Ocean-Wave_Spectra) wave spectrum
- `surge_max`: extreme surge response calculation by [IFFT](https://www.mathworks.com/help/matlab/ref/ifft.html)

### post-processing
- `plot_exprob_mcs.m`: exceedance probability plots for MCS and IS results
- `plot_exprob_pce.m`: exceedance probability plots for MCS and PCE results

## Example
1. Clone this repository to your directory
2. Open `MCS_30min.m` in MATLAB and run as it is or you can change variables (ex: number of samples, `N_T`)
3. You can also run the others: `IS_30min.m` and `pce.m`
4. The results after the runs are saved in `.\mat_file\`
5. You can create exceedance probability plots using the saved .mat files by `plot_exprob_mcs.m` and `plot_exprob_pce.m`

## Figures
<img src="https://github.com/hyeonguklim/longterm_extreme/blob/master/figures/exprob_mcs.png" width="350" height="">

**Fig**: Exceedance probability estimation by MCS and importance sampling

Variance in exceedance probability estimation is reduced with the help of importance sampling.

<img src="https://github.com/hyeonguklim/longterm_extreme/blob/master/figures/exprob_pce.png" width="350" height="">

**Fig**: Exceedance probability estimation by MCS and PCE

## Contact
For any questions or comments, please email me at: hyeonguklim@gmail.com.

