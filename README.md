# CFFDRS-FWI_India - Gridded Implementation of CFFDRS-FWI System for India
Computes and plots Fire weather indices over India / tropical climate

This repository provides a MATLAB-based gridded implementation of the Canadian Forest Fire Danger Rating System (CFFDRS-FWI) customized for use over the Indian region at 10 km spatial resolution. The implementation retains the original scientific formulation of the FWI system but adapts specific components to ensure regional relevance and computational scalability.

The CFFDRS-FWI system is a widely used model for assessing forest fire danger based on meteorological inputs. However, the original implementation is point-based and optimized for Canadian latitudes. To address limitations in computational efficiency and to support high-resolution analysis over India, this project:

Rewrites the model in MATLAB for gridded data support.

Enables fast, parallel computation over large spatial domains.

Calibrates latitude-sensitive components for Indian conditions.

Key Features
FFMC (Fine Fuel Moisture Code): Retained without modification as it responds directly to meteorological inputs.

DMC (Duff Moisture Code) & DC (Drought Code):
Adjusted to use spatially explicit day-length values derived from latitude and day-of-year at each grid point, replacing fixed daylight factors.

PET (Potential Evapotranspiration) for DC:
Modified to reflect Indian tropical conditions:

For latitudes <20°N: consistent heating assumed throughout the year.

For higher latitudes: monthly PET variation considered.

Applications
Forest fire danger estimation under current and future climate scenarios.

Integration with climate model outputs (e.g., CESM, CMIP6).

Regional forest fire monitoring

Disclaimer
This implementation maintains the original scientific formulation of the CFFDRS-FWI system. No changes were made to danger threshold values or empirical coefficients. Calibration steps were only applied to enhance the spatial relevance and computational feasibility for Indian conditions.

