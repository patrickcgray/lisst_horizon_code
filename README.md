This repository contains the code for the manuscript "*The LISST-Horizon: continuous measurements of polarized and unpolarized angular scattering*" under review at Applied Optics.

Authors: Patrick Clifton Gray<sup>\*,1</sup>, Wayne Slade<sup>2</sup>,  Daniel Koestner<sup>3</sup>, Margaret Estapa<sup>1</sup>, Thomas Leeuw<sup>4</sup>, Kirby Simon<sup>4</sup>, Emmanuel Boss<sup>1</sup> 

1. School of Marine Sciences, University of Maine, Orono, ME, USA
2. Harbor Branch Oceanographic Institute, Florida Atlantic University, Fort Pierce, FL, USA
4. University of Bergen, Bergen, Norway
5. Sequoia Scientific, Inc., Bellevue, WA 98005, USA

*Corresponding authors: Patrick Gray (patrick.gray@maine.edu)

We introduce the LISST-Horizon, a new bench-top instrument for measuring angular polarized scattering which provides the volume scattering function from 0.05° to 150° and the degree of linear polarization from 30° to 150°. This instrument can measure continuously for extended periods in a flow-through setting with minimal operator intervention. We present the instrument design, calibration procedure, processing of the data, and validation. Following calibration, the instrument shows excellent agreement with Mie theory for polystyrene beads and with the LISST-VSF for both beads and natural samples, though with higher uncertainty in the back angles (120°-150°) due to internal reflections. The instrument calibration factors are stable over 5 years including one full year where it was deployed at sea and collecting data continuously while the ship was underway. We suggest that this instrument provides a valuable new addition to the optical oceanographer's toolkit by rapidly measuring angular and polarized scattering properties, supporting marine particle characterization across space and time, and validating ocean products from new satellite multi-angle polarimeters such as PACE's HARP2 and SPEXone.

## Code access

hrzn_utils.py contains the code to process the LISST-Horizon raw data into calibrated VSF and DoLP measurements. The use of this code is shown in a few Jupyter notebooks included here:
1.  `primary_paper_figs_and_analysis.ipynb` recreates nearly all analysis and figures from the paper.
2.  `lisst_fwd_optics_cal.ipynb` shows a simple validation of the forward scattering array.
3.  `nd_filter_test.ipynb` compares data collected with and without a neutral density filter to eliminate the internal reflection in the instrument.
4.  `sample_tara_data.ipynb` processes and visualizes a section of the data collected by the LISST-Horizon while on a year long oceanographic expedition.

All python code can be run on Docker and is likely to work on most scientific computing environments but an exactly a matching environment can be obtained by copying the CryoCloud Docker image: https://github.com/CryoInTheCloud/hub-image.

# Data Access
All data is included in this repository in the data directory except for the Tara ancillary data which can be accessed at https://preprints.opticaopen.org/s/6c2eec720ffe27913fff or https://github.com/patrickcgray/spatial_patchiness_tara/releases/tag/v1.0.
