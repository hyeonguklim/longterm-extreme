# MCS (importance sampling as well) and PCE for long-term extreme response prediction
## Summary
![](https://github.com/hyeonguklim/longterm_extreme/blob/master/scheme/scheme.png)  
**Fig**: A simple SDOF structural system subject to wave loading and a quantity of interest (QoI)

## Authors/Collaborators
[HyeongUk Lim](https://hyeonguk.wordpress.com), [Lance Manuel](https://lancemanuel.netlify.com), and [Ying Min Low](https://www.eng.nus.edu.sg/cee/staff/low-ying-min/)

## Description
This study focuses on the development of efficient surrogate models by polynomial chaos expansion (PCE) for prediction of the long-term extreme surge motion of a moored floating offshore structure. The structure is subjected to first-order and second-order (difference-frequency) wave loading. Uncertainty in the long-term response arises from contrasting sea state conditions, characterized by significant wave height, Hs, and spectral peak period, Tp, and their relative likelihood of occurrence; these two variables are explicitly included in the PCE-based uncertainty quantification (UQ). In a given sea state, however, response simulations must be run for the associated Hs and Tp; in such simulations, typically, a set of random amplitudes and phases define an irregular wave train consistent with that sea state.  These random amplitudes and phases for all the frequency components in the wave train introduce additional uncertainty in the simulated waves and in the response.  The UQ framework treats these two sources of uncertainty---from Hs and Tp on the one hand, and the amplitude and phase vectors on the other---in a manner that efficiently yields long-term surge motion extreme predictions consistent with more expensive Monte Carlo simulations (MCS) that serve as the "truth" system. To reduce uncertainty in response extremes that result from sea states with a low likelihood of occurrence, importance sampling is employed with both MCS- and PCE-based extreme response predictions. Satisfactory performance with such efficient surrogate models can help in assessing the long-term response of various offshore structures.  

## Related Publications/Presentations
- HyeongUk Lim, Lance Manuel and Ying Min Low (2018), On Efficient Long-term Extreme Response Estimation for a Moored Floating Structure, Madrid, Spain, June 17-22, Proceedings of the ASME 2018 37th International Conference on Ocean, Offshore and Arctic Engineering, OMAE2018-78763. [[presentation]](https://github.com/hyeonguklim/longterm_extreme/blob/master/presentation/OMAE_2018-78763.pdf) [[paper]](https://doi.org/10.1115/OMAE2018-78763)

## MATLAB Codes
Short explanation of MATLAB codes. You can also find line-by-line explanations inside the codes.

### data files

- `data.mat`: structural/hydrodynamic parameters and initial setting for the analysis; please see `input_desc` in the file for more information. For example, type the following lines in Matlab:

```Matlab
load data.mat input_desc
disp(input_desc)
```

### main functions

- `MCS_30min.m`: Monte Carlo simulation of 30-min extreme surge response
- `IS_30min.m`: MCS with the help of importance sampling

### subfunctions
The subfuctions needed for MCS and IS are attached at the bottom of the main functions:

- `incdfHs`: inverse cumulative distribution function (cdf) of Hs
- `pdfHs`: probability density function (pdf) of Hs
- `pdfHsIS`: importance sampling density function of Hs
- `incdfTp`: inverse cdf of Tp
- `Jonswap`: [JONSWAP](https://wikiwaves.org/Ocean-Wave_Spectra) wave spectrum
- `surge_max`: extreme surge response calculation by [IFFT](https://www.mathworks.com/help/matlab/ref/ifft.html)

### post-processing
- `plot_exprob.m`: exceedance probability plots for MCS and IS results

## How to Run an Example
1. Clone this repository to your directory
2. Run `aPCE_Ishigami.m` in MATLAB
3. You can change parameters, e.g, a polynomial order (`p`)
4. You will get an exceedance probability plot

## Sample Figures
![](https://github.com/hyeonguklim/aPCE/blob/master/figures/exceedance_plot.png)  
**Fig**: Exceedance probability estimation by aPCE for the Ishigami function

Ten sets of order-8 aPCE surrogate models estimate exceedance probabilities well when compared with ten sets of Monte Carlo simulations in the Ishigami function.

## Contact
For any questions or comments, please email me at: hyeonguklim@gmail.com.

