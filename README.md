## GeoScience Ionospheric Tool (GSIT)

### Overview
This is a Python software framework for ionospheric community, which includes data processing processing tools for GNSS, optical and magnetometer sensors. 

### GNSS
`pyGps.py` GNSS data processing tool includes:
- dual frequency pseudo-range TEC
- dual frequency phase TEC
- dual frequency phase corrected TEC
- rate of TEC (ROTI)
- vertical TEC with respect to thin shell approximation
- Ionospheric Piercing Point (IPP)
- Satellite position in AER and WGS84

### Optical data
`asi.py` Optical data procesing is done with a GeoData object by @jswoboda.
Optical data tool includes:
- All-sky Imager intensity with respect to IPP
- Keogram

### 3-axis earth-based magnetometer
`magnetometer.py` Magentometer data processing tool includes:
- Conversion tool between HDZ and XYZ rfame of reference

### Space weather data
`imf.py` Reads and plots space weather data from NASA Goddard space Flight Center OMNIWeb page
