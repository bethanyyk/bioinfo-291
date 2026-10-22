# Domestic dogs (*Canis familiaris*) distinguish conspecific emotional chemosignals

## Description of the data and file structure

Data was collected from 43 dogs with owners at Barnard College's Dog Cognition Lab. Odor stimuli was from an unfamiliar dog in three induced emotional states: joy (*J*), stress (*S*), and baseline (*B*). Part 1: habituation-discrimination paradigm tested dogs’ ability to distinguish these emotional odors in conspecifics. Part 2: data was collected from a subset of 24 dogs, assessing their behavioral reactivity to conspecific emotional odors with their owner and a stranger present; measured stress/joy behaviors, proximity to owner/stranger/odor, and body posture. Videos were coded in BARC and analyzed in R. We provide raw datasets (Part_1.xlsx, Part_2.xlsx) and an R script (statistical_analyses.R) for reproducing all analyses conducted in this study using these spreadsheets.

## Files and variables

#### **File: Part_1.xlsx**

Description: Subject dog demographics and investigation times (in seconds) in habituation-discrimination trials.

Cells labeled "null" represent data that was not collected and are   intentionally missing from the spreadsheet. Red-highlighted "null" cells mean the dog did not participate in that portion of the study. Unhighlighted "null" cells indicate data that was not applicable or never intended to be collected (e.g., in a joy vs. baseline comparison from the *specificorders* sheet, only one condition can serve as the habituated odor, so the corresponding columns for the other condition remain "null").

Demographic variables:

* *Season*: Data collection period (April 2024, n=14; November 2024, n=17; and February–March 2025, n=12)
* *Subject*: Anonymized ID for each subject dog
* *Sex*: Subject dog’s sex
* *Birthday*: Owner-reported date of birth of subject dog
* *Age_months*: Age of subject dog in months as calculated from their birthday
* *Spayedneutered*: Spay/neuter status of subject dog
* *Breed*: Mixed-breed or purebred (owner-reported)
* *Breed_specific*: Specific breed (owner-reported)

Emotions are coded as *J*=joy, *S*=stress, *B*=baseline.

Sheets overview:

* *specificorders*: investigation times (sec) organized by emotional pairing and presentation order. Examples:
  * *Jhabit1_vB*: time investigating first presentation of joy as the habituated odor (*Jhabit1*) in a joy vs. baseline comparison (*_vB*; ‘versus B’)
  * *Bnovel_vJ*: time investigating baseline as the novel odor (*Bnovel*) in a joy vs. baseline comparison (*_vJ*)
  * *Bside_vJ*: side (left/right) of baseline as the novel odor (*Bside*) in joy vs. baseline comparison (*_vJ*)

- *emotionalpairings*: investigation times (sec) organized by emotional pairing only (order ignored).
  Examples:
  * *JBhabit3*: time investigating third habituation presentation (*habit3*) in a joy vs. baseline comparison (*JB*)
  * *JBnovel*: time investigating the novel odor in a joy vs. baseline comparison
  * *JBside*: side of novel odor in a joy vs. baseline comparison
  * *JBdishabitdifference*: difference between time investigating the novel odor (*JBnovel*) and the third habituation presentation (*JBhabit3*) in a joy vs. baseline comparison
  * *JBdiscrmdifference*: difference between time investigating the novel odor (*JBnovel*) and the simultaneous habituation presentation (*JBhabit4*) in a joy vs. baseline comparison
- *all*: investigation times (sec) for each presentation, pooled across all emotional odor pairings
  * *habit1*: time investigating first habituation presentation
  * *habit2*: time investigating second habituation presentation
  * *habit3*: time investigating third habituation presentation
  * *habit4*: time investigating habituated odor presented alongside the novel odor
  * *novel*: time investigating novel odor presented alongside the habituated odor (i.e., fourth habituation presentation)

#### **File: Part_2.xlsx_**

Description: Behavioral frequencies in subject dogs 1 min after exposure to conspecific emotional odors: joy (*J*), stress (*S*), or control (empty canister; *C*). Recorded every 5 sec and then averaged across dogs for each condition. 

Cells labeled "null" represent data that could not be collected (e.g., due to early trial termination) and are intentionally missing from the spreadsheet.

Raw data variables:

* *Subject*: Anonymized ID for each subject dog (consistent with those used in Part 1)
* *Sex*: Subject dog’s sex
* *Odor:* Conspecific emotional odor condition (*J*=joy, *S*=stress, *C*=control)
* *5s – 60s*: Behavioral frequencies recorded in 5-second intervals. Each column represents the frequency of a behavior observed during the 5 seconds leading up to the indicated time point (e.g., “*10s*” = frequency between 5s to 10s after the dog investigated the odor stimulus). See Data Analysis – Part 2 in the Methods for details on the definition of frequency for each behavioral category. 

Processed data variables (used for statistical analyses):

* *Condition*: Conspecific emotional odor condition (*J*=joy, *S*=stress, *C*=control)
* *Frequency*: Mean behavioral frequency across all dogs for a sample point (i.e., 5-second interval).

Sheets overview: Each sheet represents a specific behavioral category.

* *stressbehav*: stress behaviors
* *joybehav*: joy behaviors
* *owner*: proximity to owner
* *box*: proximity to odor box
* *stranger*: proximity to stranger
* *upright*: upright posture
* *sit*: sitting posture
* *prone*: lying down posture
* For sex-specific data, sheet names containing "F" indicate female dog behavioral frequencies, "M" for male dogs. Columns with names containing "_MF" are for assessing sex × odor condition interaction on behavior.

#### **File: statistical_analyses.R**

Description: Statistical analyses conducted on the datasets (Part_1.xlsx, Part_2.xlsx) performed in R.

## Code/software

BARC v.1.0.0

RStudio 2023.06.2+561
