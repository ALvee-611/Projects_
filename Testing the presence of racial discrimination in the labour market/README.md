Testing the presence of racial discrimination in the labour market
---

The dataset **reg_data** is a subset of the dataset used by M. Bertrand and S. Mullainathan for their research paper "Are Emily and Greg More Employable Than Lakisha and Jamal? A Field Experiment on Labour Market Discrimination", published in 2004 in the American Economic Review. It contains the following 27 variables:

* **name**: Factor indicating the applicant's first name.
* **gender**: Factor indicating gender.
* **ethnicity**: Factor indicating ethnicity (i.e., Caucasian-sounding vs. African-American sounding first name).
* **quality**: Factor indicating the quality of the resume.
* **call**: Was the applicant called back?
* **city**: Factor indicating city: Boston or Chicago.
* **jobs**: Number of jobs listed on the resume.
* **experience**: Number of years of work experience on the resume.
* **honors**: Did the resume mention some honors?
* **volunteer**: Did the resume mention some volunteering experience?
* **military**: Does the applicant have military experience?
* **holes**: Does the resume have some employment holes?
* **school**: Does the resume mention some work experience while at school?
* **email**: Was the e-mail address on the applicant's resume?
* **computer**: Does the resume mention some computer skills?
* **special**: Does the resume mention some special skills?
* **college**: Does the applicant have a college degree or more?
* **minimum**: Factor indicating the minimum experience requirement of the employer.
* **equal**: Is the employer EOE (equal opportunity employment)?
* **wanted**: Factor indicating the type of position wanted by the employer.
* **requirements**: Does the ad mention some requirement for the job?
* **reqexp**: Does the ad mention some experience requirement?
* **reqcomm**: Does the ad mention some communication skills requirement?
* **reqeduc**: Does the ad mention some educational requirement?
* **reqcomp**: Does the ad mention some computer skills requirement?
* **reqorg**: Does the ad mention some organizational skills requirement?
* **industry**: Factor indicating the type of employer industry.

The purpose of this project is to test the presence of racial discrimination in the labour market. In this project, I formulated 4 different questions related to discrimination and build one model for each question. The questions were:

1. Is discrimination different for people from Boston and Chicago?

2. Is discrimination different for male and female individuals with prior military experience having African-American sounding name?

3. Is the effect of experience on the probability of being called back the same for people provided email address and those that didnot?

4. Is discrimination different for people who mention some volunteering experience VS those that didnot for people applying to EOE?

A detailed analysis of this can be found here: https://alvee-611.github.io/portfolio/Regression_model_2/
