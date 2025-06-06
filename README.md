# Executable Environment for OSF Project [ecjy6](https://osf.io/ecjy6/)

This repository was automatically generated as part of a project to test the reproducibility of open science projects hosted on the Open Science Framework (OSF).

**Project Title:** How to disentangle psychobiological stress reactivity and recovery: A comparison of model-based and non-compartmental analyses of the cortisol concentrations

**Project Description:**
> This article seeks to address the prevailing issue of how to measure specific process components of psychobiological stress responses. Particularly the change of cortisol secretion due to stress exposure has been discussed as an endophenotype of many psychosomatic health outcomes. To assess its process components, a large variety of non-compartmental parameters (i.e., composite measures of substance concentrations at different points in time) like the area under the concentration-time curve (AUC) are commonly utilized. However, a systematic evaluation and validation of these parameters based on a physiologically plausible model of cortisol secretion has not been performed so far.
Thus, a population pharmacokinetic (mixed-effects stochastic differential equation) model was developed and fitted to densely sampled salivary cortisol data of 10 males from Montreal, Canada, and sparsely sampled data of 200 mixed-sex participants from Dresden, Germany, who completed the Trier Social Stress Test (TSST). Besides the two major process components representing (1) stress-related cortisol secretion (reactivity) and (2) cortisol elimination (recovery), the model incorporates two additional, often disregarded components: (3) the secretory delay after stress onset, and (4) deviations from the projected steady-state concentration due to stress-unrelated fluctuations of cortisol secretion.
The fitted model (R2 = 99%) was thereafter used to investigate the correlation structure of the four individually varying, and readily interpretable model parameters and eleven popular non-compartmental parameters. Based on these analyses, we recommend to use the minimum-maximum cortisol difference and the minimum concentration as proxy measures of reactivity and recovery, respectively. Finally, statistical power analyses of the reactivity-related sex effect illustrate the consequences of using impure non-compartmental measures of the different process components that underlie the cortisol stress response.

**Original OSF Page:** [https://osf.io/ecjy6/](https://osf.io/ecjy6/)

---

**Important Note:** The contents of the `ecjy6_src` folder were cloned from the OSF project on **12-03-2025**. Any changes made to the original OSF project after this date will not be reflected in this repository.

The `DESCRIPTION` file was automatically added to make this project Binder-ready. For more information on how R-based OSF projects are containerized, please refer to the `osf-to-binder` GitHub repository: [https://github.com/Code-Inspect/osf-to-binder](https://github.com/Code-Inspect/osf-to-binder)

## flowR Integration

This version of the repository has the **[flowR Addin](https://github.com/flowr-analysis/rstudio-addin-flowr)** preinstalled. flowR allows visual design and execution of data analysis workflows within RStudio, supporting better reproducibility and modular analysis pipelines.

To use flowR, open the project in RStudio and go to `Addins` > `flowR`.

## How to Launch:

**Launch in your Browser:**

🚀 **MyBinder:** [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/code-inspect-binder/osf_ecjy6-f/HEAD?urlpath=rstudio)

   * This will launch the project in an interactive RStudio environment in your web browser.
   * Please note that Binder may take a few minutes to build the environment.

🚀 **NFDI JupyterHub:** [![NFDI](https://nfdi-jupyter.de/images/nfdi_badge.svg)](https://hub.nfdi-jupyter.de/r2d/gh/code-inspect-binder/osf_ecjy6-f/HEAD?urlpath=rstudio)

   * This will launch the project in an interactive RStudio environment on the NFDI JupyterHub platform.

**Access Downloaded Data:**
The downloaded data from the OSF project is located in the `ecjy6_src` folder.

## Run via Docker for Long-Term Reproducibility

In addition to launching this project using Binder or NFDI JupyterHub, you can reproduce the environment locally using Docker. This is especially useful for long-term access, offline use, or high-performance computing environments.

### Pull the Docker Image

```bash
docker pull meet261/repo2docker-ecjy6-f:latest
```

### Launch RStudio Server

Run the container (with a name, e.g. `rstudio-dev`):
```bash
docker run -it --name rstudio-dev --platform linux/amd64 -p 8888:8787 --user root meet261/repo2docker-ecjy6-f bash
```

Inside the container, start RStudio Server with no authentication:
```bash
/usr/lib/rstudio-server/bin/rserver --www-port 8787 --auth-none=1
```

Then, open your browser and go to: [http://localhost:8888](http://localhost:8888)

> **Note:** If you're running the container on a remote server (e.g., via SSH), replace `localhost` with your server's IP address.
> For example: `http://<your-server-ip>:8888`

## Looking for the Base Version?

For the original Binder-ready repository **without flowR**, visit:
[osf_ecjy6](https://github.com/code-inspect-binder/osf_ecjy6)

