# 2024_pvpmc_self_shade

------
__UPDATE - 2024-05-28:__ An issue was found with the obeservation data from plant 20/21. This data and the plant specs were updated. This resulted in the `haydavies` model looking better than `isotropic`, so many of the examples for this plant and the others were updated. This results in the notebooks will now be slightly different from what was presented at PVPMC. The presentation files _have not been updated_. 

------

This repository contains code, data, and supporting files for the 2024 PV Performance Modeling Collaborative (PVPMC) presentation "An approach to modeling linear and non-linear self-shading losses with pvlib" by William B. Hobbs (Southern Company), Kevin S. Anderson (Sandia), Mark A. Mikofski (DNV), and Madison Ghiz (DNV).

You can check out the presentation [PDF](presentation/2024_PVPMC_hobbs_pvlib_self-shade.pdf) or [PPTX](presentation/2024_PVPMC_hobbs_pvlib_self-shade.pptx) files in the [presentation](presentation) folder. 

The python functions created for this project are in [self_shade.py](self_shade.py). 

To walk through examples, see these notebooks:
- [2024_pvpmc_self_shade.ipynb](2024_pvpmc_self_shade.ipynb)
- [additional_analyses.ipynb](additional_analyses.ipynb)
- [more_plots_for_the_slides.ipynb](more_plots_for_the_slides.ipynb)

---

<img src="images\pvlib_powered_logo_horiz.png" width="600"/>


This repository uses pvlib [1]. Check it out at https://pvlib-python.readthedocs.io/en/stable/ or https://github.com/pvlib/pvlib-python. 

[1] Anderson, K., Hansen, C., Holmgren, W., Jensen, A., Mikofski, M., and Driesse, A. “pvlib python: 2023 project update.” Journal of Open Source Software, 8(92), 5994, (2023). https://doi.org/10.21105/joss.05994