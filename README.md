# nimbleEco

Repository for the NimbleEcology R package as part of Google Summer of Code 2017

# Motivation

The goal of this Google Summer of Code project was to build a new R package (nimbleEcology) providing high-level user interfaces to common ecological models and implement the analyses internally using [NIMBLE](https://r-nimble.org/), a system for building and estimating Bayesian hierarchical models.  There were two primary motivations for this project. The first: to provide better tools for widely used ecological models.  The second, perhaps equally important, was to showcase how NIMBLE can serve as a computational engine “under the hood” of an R package, and to develop extensions to NIMBLE that will support its use in this way.
 
# Work Completed 

I set out to implement 3 broad classes of ecological models: abundance/occupancy models, capture/recapture models, and population dynamic models.  Of the proposed models, I was able to design, and test, functions for abundance models, single-season occupancy models, and a multi-season single species occupancy model.  After loading the nimbleEcology package into R, each of the models are run using one of three functions: nimble.abund, nimble.occ, and nimble.dynamic.occ.  The source code for each function is contained within the /R folder of the nimbleEco Git repository.  Each function allows the analyst to specify the model formula, provide data, and select a variety of options depending on the model class.  The functions then munge the data into a [tidy](http://vita.had.co.nz/papers/tidy-data.html) format, dynamically generate BUGS code to fit the model, and finally call NIMBLE to perform the estimation.  The user is returned a summary of parameter estimates, the generated BUGS code, and the full set of posterior MCMC draws.  Details about the function arguments are available in Help once the package is loaded into R.

In addition to writing the model code, I also coded unit tests for all 3 models using the [testthat](http://r-pkgs.had.co.nz/tests.html) package.  Writing unit tests using testthat is a common part of R package development.  For each of the functions, I developed a gold file with a subset of posterior draws from a variety of model specifications estimated using simulated data.  When the package is tested the models are rerun using a set seed and the new posterior draws compared against the gold file.  Each of the tests are contained in the /tests/testthat folder; the test data and gold files are in /tests/testdata.

In addition to the ecological models, along with my mentors, I helped to implement a new code expansion feature that allows for internal NIMBLE dynamic code expansion.  Currently called "model macros", this function allows standard pieces of code to get full expanded internally within NIMBLE.  Currently, this feature is only implemented within the nimbleEcology package.  The goal, however, is to get these model macros, and additional code expansion macros, into core nimble.

# Future Work

While the nimbleEcology package is fully operational and estimates 3 main classes of ecological statistical models, additional work is required to develop it to the point where nimbleEcology could be a viable alternative to unmarked, ESURGE, and other more popular ecological software.  Below are some of the next steps I see for the package: 

* Add additional ecological models, such as the Cormack-Jolly-Seber Model, Dynamic Population Models, and spatial models.
* Allow for different model parameterizations to leverage nimble's multiple MCMC algorithms.
* Further develop the model macros and make them part of nimbleEcology, as well as nimble.
* Add functionality to return the full nimbleModel object so users can simulate data, perform model comparisons, and calculate functionals of interest.

# Conclusion

I had a great experience with the Google Summer of Code program.  All of my mentors (Perry, Daniel, and Nick) were a pleasure to work with, both professionally and personally.  In addition to providing useful feedback, they really incorporated me into the NIMBLE team, which I highly valued.  My R coding skills and knowledge of ecological statistical models improved immensely over the summer.  I hope I can continue to work on the project, and would highly recommend Google Summer of Code to future statistics graduate students.






