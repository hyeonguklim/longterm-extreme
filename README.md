# MCS (importance sampling as well) and PCE for long-term extreme response prediction
## Summary
![](https://github.com/hyeonguklim/aPCE/blob/master/figures/scheme.png)  
**Fig**: A schematic of using arbitrary polynomial chaos expansion (aPCE)

## Authors/Collaborators
[HyeongUk Lim](https://hyeonguk.wordpress.com) and [Lance Manuel](https://lancemanuel.netlify.com)

## Description
This study focuses on the development of efficient surrogate models by polynomial chaos expansion (PCE) for prediction of the long-term extreme surge motion of a moored floating offshore structure. The structure is subjected to first-order and second-order (difference-frequency) wave loading. Uncertainty in the long-term response arises from contrasting sea state conditions, characterized by significant wave height, Hs, and spectral peak period, Tp, and their relative likelihood of occurrence; these two variables are explicitly included in the PCE-based uncertainty quantification (UQ). In a given sea state, however, response simulations must be run for the associated Hs and Tp; in such simulations, typically, a set of random amplitudes and phases define an irregular wave train consistent with that sea state.  These random amplitudes and phases for all the frequency components in the wave train introduce additional uncertainty in the simulated waves and in the response.  The UQ framework treats these two sources of uncertainty---from Hs and Tp on the one hand, and the amplitude and phase vectors on the other---in a manner that efficiently yields long-term surge motion extreme predictions consistent with more expensive Monte Carlo simulations (MCS) that serve as the "truth" system. To reduce uncertainty in response extremes that result from sea states with a low likelihood of occurrence, importance sampling is employed with both MCS- and PCE-based extreme response predictions. Satisfactory performance with such efficient surrogate models can help in assessing the long-term response of various offshore structures.  
[see slides](https://github.com/hyeonguklim/longterm_extreme/blob/master/presentation/OMAE_2018-78763.pdf)

## Related Publications/Presentations
- Lim, H and Manuel, L, Distribution-Free Polynomial Chaos Expansion Surrogate Models for Efficient Structural Reliability Analysis, *Engineering Mechanics Institute Conference*, Pasadena, CA, June 18-21, 2019. [[presentation](https://hyeonguk.files.wordpress.com/2019/07/emi19_presentation.pdf)]

## Codes
### examples
This folder contains examples of using aPCE:
- [Ishigami function](https://www.sfu.ca/~ssurjano/ishigami.html): `aPCE_Ishigami.m`

### subfunctions
This folder contains the subfunctions needed for running aPCE:

- `aPCE.m`: builds an aPCE model 
- `aPCE_coef.m`: calculates the coefficients of a polynomial function by the Gram-Schmidt orthogonalization
- `ishigami.m`: Ishigami function evaluation
- `multi_index.m`: gives multi-indices needed for multi-variate polynomial functions

### other functions
- `load_path.m`: sets the path where subfunctions are located

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

$-b \pm \sqrt{b^2 - 4ac} \over 2a$
$x = a_0 + \frac{1}{a_1 + \frac{1}{a_2 + \frac{1}{a_3 + a_4}}}$
$\forall x \in X, \quad \exists y \leq \epsilon$