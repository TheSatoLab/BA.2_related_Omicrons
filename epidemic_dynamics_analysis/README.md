# Modeling the epidemic dynamics of SARS-CoV-2 lineages

To quantify the spread rate of each SARS-CoV-2 lineage in the human population, we estimated the relative effective reproduction number (R<sub>e</sub>) of each viral lineage according to the epidemic dynamics, calculated on the basis of viral genomic surveillance data. The data were downloaded from the GISAID database (https://www.gisaid.org/) on May 15, 2022. We excluded the data of viral strains with the following features from the analysis: i) a lack of collection date information; ii) sampling in animals other than humans; or iii) sampling by quarantine. We analyzed the datasets of the five countries (South Africa, the USA, France, Denmark and Belgium) where BA.4/5, BA.2.12.1, BA.2.11, BA.2.9.1, and BA.2.13 were most detected, respectively (Table Sx5). The BA.2 sublineages without amino acid mutations at position 452 in S were summarized as BA.2. In addition, the Delta sublineages were also summarized as Delta. The dynamics of up to five most predominant viral lineages in each country from February 5, 2022, to May 15, 2022, were analyzed. The number of viral sequences of each viral lineage collected on each day in each country was counted, and the count matrix was constructed as an input for the statistical model below. 
We constructed a Bayesian statistical model to represent relative lineage growth dynamics with multinomial logistic regression, as described in our previous study (Suzuki et al., 2022). In the present study, the epidemic dynamics in respective countries were independently estimated. Arrays in the model index over one or more indices: viral lineages l and days t. The model is:

![\begin{align*}
\mu_{lt}=\alpha_{l}+\beta_{l}t
\end{align*}
](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0A%5Cmu_%7Blct%7D%3D%5Calpha_%7Blc%7D%2B%5Cbeta_%7Blc%7Dt%0A%5Cend%7Balign%2A%7D%0A)

![\begin{align*}
\theta_{.t}=softmax(\mu_{.t})
\end{align*}
](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0A%5Ctheta_%7B.ct%7D%3Dsoftmax%28%5Cmu_%7B.ct%7D%29%0A%5Cend%7Balign%2A%7D%0A)

![\begin{align*}
y_{lt}\sim Multinomial(\sum_{l}y_{lt},\theta_{.t})
\end{align*}
](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0Ay_%7Blct%7D%5Csim+Multinomial%28%5Csum_%7Bl%7Dy_%7Blct%7D%2C%5Ctheta_%7B.ct%7D%29%0A%5Cend%7Balign%2A%7D%0A)


The explanatory variable was time t, and the outcome variable was y<sub>lt</sub>, which represented the count of viral lineage l at time t. In the model, the linear estimator &mu;<sub>t</sub>, consisting of the intercept &alpha;<sub>l</sub> and the slope &beta;<sub>l</sub>, was converted to the simplex &theta;<sub>t</sub>, which represented the probability of occurrence of each viral lineage at time t, based on the softmax link function defined as:

$$softmax(x)=\frac{exp(x)}{\sum_{i}exp(x_i)}$$


y<sub>lt</sub> is generated from theta<sub>t</sub> and the total count of all lineages at time t according to a multinomial distribution.
	The relative R<sub>e</sub> of each viral lineage (r<sub>l</sub>) was calculated according to the slope parameter &beta;<sub>l</sub> as:

![\begin{align*}
r_{l}=exp(\gamma\beta_{l})
\end{align*}
](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0Ar_%7Blc%7D%3Dexp%28%5Cgamma%5Cbeta_%7Blc%7D%29%0A%5Cend%7Balign%2A%7D%0A)


where gamma is the average viral generation time (2.1 days) (http://sonorouschocolate.com/covid19/index.php?title=Estimating_Generation_Time_Of_Omicron).
For parameter estimation, the intercept and slope parameters of the BA.2 variant were fixed at 0. Consequently, the relative R<sub>e</sub> of BA.2 was fixed at 1, and those of the other lineages were estimated relative to that of BA.2.
Parameter estimation was performed via the MCMC approach implemented in CmdStan v2.28.1 (https://mc-stan.org) with CmdStanr v0.4.0 (https://mc-stan.org/cmdstanr/). Noninformative priors were set for all parameters. Four independent MCMC chains were run with 500 and 2,000 steps in the warmup and sampling iterations, respectively. We confirmed that all estimated parameters showed <1.01 R-hat convergence diagnostic values and >200 effective sampling size values, indicating that the MCMC runs were successfully convergent. The above analyses were performed in R v4.1.3 (https://www.r-project.org/).

## Contents:
*  **transmissibility.for_Fig.pango.R:** Main script
*  **multinomial_independent.stan:** Stan model file
*  **summarize_mut_info.py:** Python script to generate **metadata.mut_long.tsv** used in **transmissibility.for_Fig.pango.R**

## Not included:
*  **metadata.tsv:** Please download from **Download Packages** in GISAID (https://www.gisaid.org/.
*  **metadata.mut_long.tsv:** Please generate using **summarize_mut_info.py** and **metadata.tsv** as below:

```bash
python3 summarize_mut_info.py \
        metadata.tsv \
        > metadata.mut_long.tsv
```

