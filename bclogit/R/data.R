#' Framingham Heart Study Dataset
#'
#' A subset of the Framingham Heart Study data.
#' 
#' This dataset was constructed by running the following code: 
#' 
#' ```r
#' pacman::p_load(riskCommunicator, data.table)
#' data("framingham")
#' D = data.table(framingham)
#' D = D[!is.na(CIGPDAY)] #we drop missing data in covariates
#' D = D[!is.na(BMI)]
#' D = D[!is.na(HEARTRTE)]
#' D = D[!is.na(TOTCHOL)]
#' Dba = D[PERIOD %in% c(1,3)] #we drop intermediate periods so we have matched pairs
#' Dba[PERIOD == 1, PERIOD := 0]
#' Dba[PERIOD == 3, PERIOD := 1]
#' Dba[, num_periods_per_id := .N, by = RANDID]
#' Dba = Dba[num_periods_per_id == 2] #we drop intermediate periods so we have matched pairs
#' Dba[, num_periods_per_id := NULL]
#' ```
#' 
#' @format A data frame with 5944 rows and 39 variables:
#' \describe{
#'   \item{RANDID}{Unique identification number for each participant. (This is the strata for the matched pairs).}
#'   \item{PERIOD}{Examination Cycle where 0 = baseline, 1 = endpoint (the treatment variable)}
#'   \item{SEX}{Participant sex (1 = Male, 2 = Female).}
#'   \item{TOTCHOL}{Total serum cholesterol (mg/dL).}
#'   \item{AGE}{Age at exam (years).}
#'   \item{SYSBP}{Systolic blood pressure (mmHg).}
#'   \item{DIABP}{Diastolic blood pressure (mmHg).}
#'   \item{CURSMOKE}{Current smoking status (0 = No, 1 = Yes).}
#'   \item{CIGPDAY}{Number of cigarettes smoked per day.}
#'   \item{BMI}{Body Mass Index (kg/m^2).}
#'   \item{DIABETES}{Diabetes status (0 = No, 1 = Yes).}
#'   \item{BPMEDS}{Use of Anti-hypertensive medication (0 = No, 1 = Yes).}
#'   \item{HEARTRTE}{Heart rate (beats/minute).}
#'   \item{GLUCOSE}{Fast blood glucose (mg/dL).}
#'   \item{educ}{Education level.}
#'   \item{PREVCHD}{Prevalence of Coronary Heart Disease.}
#'   \item{PREVAP}{Prevalence of Angina Pectoris.}
#'   \item{PREVMI}{Prevalence of Myocardial Infarction.}
#'   \item{PREVSTRK}{Prevalence of Stroke.}
#'   \item{PREVHYP}{Prevalence of Hypertension.}
#'   \item{TIME}{Number of days since baseline exam.}
#'   \item{HDLC}{High Density Lipoprotein Cholesterol.}
#'   \item{LDLC}{Low Density Lipoprotein Cholesterol.}
#'   \item{DEATH}{Death status.}
#'   \item{ANGINA}{Angina Pectoris status.}
#'   \item{HOSPMI}{Hospitalized Myocardial Infarction.}
#'   \item{MI_FCHD}{Myocardial Infarction or Fatal Coronary Heart Disease.}
#'   \item{ANYCHD}{Any Coronary Heart Disease event.}
#'   \item{STROKE}{Stroke status.}
#'   \item{CVD}{Cardiovascular Disease status.}
#'   \item{HYPERTEN}{Hypertension status.}
#'   \item{TIMEAP}{Time to Angina Pectoris.}
#'   \item{TIMEMI}{Time to Myocardial Infarction.}
#'   \item{TIMEMIFC}{Time to Myocardial Infarction or Fatal CHD.}
#'   \item{TIMECHD}{Time to Any CHD.}
#'   \item{TIMESTRK}{Time to Stroke.}
#'   \item{TIMECVD}{Time to Cardiovascular Disease.}
#'   \item{TIMEDTH}{Time to Death.}
#'   \item{TIMEHYP}{Time to Hypertension.}
#' }
#' @source \url{https://biolincc.nhlbi.nih.gov/teaching/}
"fhs"
