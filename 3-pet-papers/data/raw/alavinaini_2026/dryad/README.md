# Data and code from: Diel activity patterns of bobcats and domestic cats in the Houston metropolitan area

Dataset DOI: [10.5061/dryad.s4mw6m9jg](10.5061/dryad.s4mw6m9jg)

## Description of the data and file structure

This dataset supports a study examining diel activity patterns and temporal niche partitioning between bobcats (*Lynx rufus*) and domestic cats (*Felis catus*) along 2 transects in the Houston, Texas metropolitan area from October 1, 2020 to August 31, 2024. Motion-activated trail cameras were deployed at 33 unique sites. Each detection event includes timestamp data recorded in Universal Coordinated Time (UTC) and processed to local standard time (US Central Standard Time, CST). All data processing and statistical analysis were conducted using MATLAB and R.

This dataset includes two CSV files listing raw and processed timestamps of bobcat and domestic cat detections, site names and coordinates, a MATLAB script to plot polar histograms of diel activity and calculate Watson's U^2^, and an R script to create temporal overlap plots and calculate overlap coefficients.

## File list

*  Bobcat_Data.csv
*  Domestic_Cat_Data.csv
*  Watsons U2 and Polar Plot.m
*  Dhat R Script.R
*  suntime_overlap_both_sites.png
*  suntime_overlap_solo_sites.png

## File Descriptions

### Bobcat_v_Domestic_Cat_MicroPub_7_31_26.zip

#### File: Bobcat_Data.csv

**Description:** Timestamp and location data for bobcat observations from October 1, 2020, to August 31, 2024.

#### File: Domestic_Cat_Data.csv

**Description:** Timestamp and location data for domestic cat observations from October 1, 2020, to August 31, 2024.  

##### Variables in Bobcat Timestamp Data and in Domestic Cat Timestamp Data

* PhotoName: unique identifier for each animal detection listed as the name of the JPG file.
* Timestamp (UTC): Original Date and Time in YEAR-DAY-MONTH and Time (T) HOUR:MINUTE:SECOND in Coordinated Universal Time (UTC).
* UTC Date: Calendar date of observation (text extracted from Timestamp (UTC)); YEAR-MONTH-DAY.
* UTC Time: Universal Coordinated Time of observation (text extracted from Timestamp (UTC)); HOUR:MINUTE:SECOND.
* UTC DateValue: Microsoft Excel serial number for the UTC date, calculated as days since January 1, 1900
* UTC TimeValue: text value for time converted to fraction of a 24-hr day
* UTC DateTimeValue: time of day on a particular date calculated as the sum of UTC DateValue and UTC TimeValue
* Local TimeValue: local time of day calculated as UTC TimeValue minus 0.25, the difference between US Central Standard Time and Universal Coordinated Time, expressed as a decimal day.
* Local TimeValue(ABS): Local TimeValue adjusted to  a value between 0 and 1.  Negative values are adjusted to a fraction of the previous day.  Example:  LocalTimeValue = -0.21indicates midnight minus 0.21 days = 0.79 of the previous day.
* Local DateTime Value: UTC DateTimeValue adjusted to local standard time, US Central Standard Time (Time zone for Houston, Texas, USA).
* Local DateTime: local date and time expressed as a text string
* Species: Common name of the species.
* SpeciesDetail: Extra information about the individual detected.
* Count: Number of that species in the photo.
* Status: Confirmation status of the species identification.  All photos  in the timestamp data sets are "confirmed" detections. Either  two separate volunteers tagged the animal as that species or an expert reviewer identified the animal.
* SiteAbbr: Site name where photo was taken (abbreviation).
* Latitude: Latitude of the Harris County 1910 Courthouse as a representative point in the center of the Houston metropolitan area.
* Longitude: Longitude of the Harris County 1910 Courthouse as a representative point in the center of the Houston metropolitan area.
* Visit ID: Unique identifier assigned to each instance of a camera being visited to check batteries and collect an SD card.
* Visit Datetime (UTC): Date and time when a visit occurred.
* Lure: Indicates whether bait was used. "Null" and "No lure" indicate that no bait was used.
* Camera Sensitivity: Sensitivity of camera to detect a difference between air temperature and animal temperature. (Set to normal for all cameras).

#### File: Watsons U2 and Polar Plot.m

Author: Natalie M. Linde

Contact: [nmlinde@uh.edu](mailto:nmlinde@uh.edu)

MATLAB Version: 24.1.0.2537033 (R2024a)

**Overview:**

This script compares diel activity patterns between bobcats (*Lynx rufus*) and domestic cats (*Felis catus*) using timestamped data. Timestamps were converted to minutes of a 1440-minute day, then to a fraction of the day (0 to 1), then  to radians before testing. The analysis includes: polar histograms of activity patterns and Watson’s U² test to compare distributions.

**Data Preparation:**

Time values were manually curated from the original bobcat and domestic cat CSV files based on site classification (co-occurrence vs. solo (only one species occurred at a particular site)), making the data smaller, fixed subsets of the complete dataset.

No required files.

**Statistical Comparison:**

A user-defined function implements Watson’s U² test to assess differences in circular time distributions between:

1. Bobcat vs. cat detections at shared (co-occurrence) sites
2. Bobcat vs. cat detections at solo (independent occurrence) sites

**Visualizations:**

Two polar histograms are generated: One histogram for solo sites and one for co-occurrence sites.  Time is plotted clockwise from midnight (0) with 24 hourly bins.

Color:

- Bobcat = Blue
- Domestic Cat = Magenta

**Citations:**

MATLAB:

The MathWorks Inc. (2024). MATLAB version: 24.1.0.2537033 (R2024a).

Polar Visualization Toolbox:

Berens, P. (2009). CircStat: A MATLAB Toolbox for Circular Statistics.

Journal of Statistical Software, 31(10), 1–21. [https://doi.org/10.18637/jss.v031.i10](https://doi.org/10.18637/jss.v031.i10)

#### File: Dhat R Script.R

Author: Leila Alavi Naini

**Description:** The script converts UTC timestamps to local time, transforms detection times to solar time using the sunTime() function in the overlap package, classifies sites as solo or shared, and calculates Dhat values for bobcat/domestic cat temporal overlap.

This dataset contains detection timestamp data for bobcats and domestic cats.

**Files:**

- Bobcat_Data.csv:   Bobcat detection records.
- Domestic_Cat_Data.csv:   Domestic cat detection records.

**To run the analysis:**

1. Open the R script in RStudio.
2. Place the script, Bobcat_Data.csv, and Domestic_Cat_Data.csv in the same folder.
3. Set the working directory to the script location.
4. Run the script from top to bottom.

#### File: suntime_overlap_both_sites.png

**Description:** Density plot depicting activity patterns of bobcats and domestic cats at  sites where both species occur (both), with the shaded area representing temporal overlap and a Dhat1 overlap coefficient of 0.646.

#### File: suntime_overlap_solo_sites.png

**Description:** Density plot depicting activity patterns of bobcats and domestic cats at sites where they occur independently (solo), with the shaded area representing temporal overlap and a Dhat1 overlap coefficient of 0.852.

## Map generation

We used Google Earth Pro on Desktop to create a map showing the location of all 33 sites included in Bobcat_Data.csv and Domestic_Cat_Data.csv.

**Generation of a map inset showing the location of Texas, USA using Gemini 3.5 Flash:** 

Prompt Used: “make a black and white map inlay of Texas within the US. Bold Texas' borders and give only Texas a darker hue with a darker border. Use WGS 84 datum to make the map but don't indicate this on the map. the scale bar should be a single line with end brackets, km not miles. Put a north indicator next to the scalebar. Put the scale bar and north indicator on the bottom left. Provide a 8k image PNG format and put a red square with sharp edges over Houston. Give the map a medium gray background.  All other states other than Texas should be the same shade. Texas should be a darker shade.”  

The generated image was then downloaded in 8k PNG format and uploaded into Microsoft PowerPoint and retouched as needed.