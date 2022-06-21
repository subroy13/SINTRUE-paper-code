# A Generalized Epidemiological Model for COVID-19 with Dynamic and Asymptomatic Population

## Abstract
In this paper we develop an extension of compartmental epidemiological models which is suitable for COVID-19. 
The model presented in this paper comprises seven compartments in the progression of the disease. This model, named as the SINTRUE (**S**usceptible, **I**nfected and pre-symptomatic, Infected and Symptomatic but **N**ot Tested, **T**ested Positive, Recorded **R**ecovered, **U**nrecorded Recovered, and **E**xpired) model. The proposed model incorporates transmission due to asymptomatic carriers and captures the spread of the disease due to the movement of people to/from different administrative boundaries within a country. In addition, the model allows estimating the number of undocumented infections in the population and the number of unrecorded recoveries. The associated parameters in the model can help architect the public health policy and operational management of the pandemic. The results show that the testing rate of the asymptomatic patients is a crucial parameter to fight against the pandemic. The model also is shown to have significantly better predictive capability than the other epidemiological models. 

## Code Description

All the codes are written in `R` programming language.

1. `main.R` - Contains the main code for the parameter estimation for the SINTRUE model for Chattishgarh's data.
2. `get_data.R` - The data cleaning and preprocessing code.
3. `estimation.R` - The partial parameter estimation code, to calculate the fatality and survival rates for different states of India for Covid-19.
4. `Italy_prediction.R` - Contains the code for the comparison of SINTRUE model with SIDARTHE model for covid-19 data of Italy.
5. `residual_diagnostics.R` - Contains the code for residual diagnostics, verifying assumptions of the SINTRUE model etc.
6. `datasets/` - folder contains the datasets used for the parameter estimation.
7. `figures/` - folder contains the generated plots used in the paper.


## Authors
- Anirban Ghatak
- Shivshanker Singh Patel
- Soham Bonnerjee
- Subhrajyoty Roy
