# Generalized-additive-models-to-capture-the-death-rates-in-Canada-COVID-19

https://arxiv.org/pdf/2008.01030.pdf

Abstract. To capture the death rates and strong weekly, biweekly and probably monthly patterns in the Canada COVID-19, we utilize the generalized
additive models in the absence of direct statistically based measurement of
infection rates. By examining the death rates of Canada in general and Quebec, Ontario and Alberta in particular, one can easily figured out that there
are substantial overdispersion relative to the Poisson so that the negative binomial distribution is an appropriate choice for the analysis. Generalized additive
models (GAMs) are one of the main modeling tools for data analysis. GAMs
can efficiently combine different types of fixed, random and smooth terms in
the linear predictor of a regression model to account for different types of effects. GAMs are a semi-parametric extension of the generalized linear models
(GLMs), used often for the case when there is no a priori reason for choosing a
particular response function such as linear, quadratic, etc. and need the data
to ’speak for themselves’. GAMs do this via the smoothing functions and take
each predictor variable in the model and separate it into sections delimited
by ’knots’, and then fit polynomial functions to each section separately, with
the constraint that there are no links at the knots - second derivatives of the
separate functions are equal at the knots.


