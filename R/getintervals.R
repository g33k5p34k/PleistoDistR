
#' Generate interval file, binning by time
#'
#' This function generates an interval file, which is the master file for controlling all downstream PleistoDist calculations. This file essentially simplifies Pleistocene-era sea level change into a number of discrete time intervals.
#'
#' @param time Maximum time cutoff in kya. PleistoDist calculates the average distance between islands over a user-defined period of time. Users will have to specify an upper time bound (in thousands of years [kya]) for their PleistoDist analysis, which by default ranges from 0.1 kya to 3000 kya (i.e. 100 to 3,000,000 years ago). The lower time bound is fixed at the present day (0 kya).
#' @param intervals Number of time interval bins to generate between the present day and the time-cutoff. More time intervals will lead to more fine-grained geomorphology calculations, at the expense of increased computational time. Under default conditions, the number of intervals should not exceed time-cutoff*10.
#' @param sealvl Optional parameter defining the sea level reconstruction file. By default, PleistoDist uses the sea level reconstruction from Bintanja & van de Wal (2008). If you use a custom sea level reconstruction file, make sure the file has a Time column (in kya) and a Sealevel_Corrected column (in metres), and is saved in CSV format.
#' @return An interval file in the automatically-generated 'output' folder. This interval file is needed for all downstream PleistoDist analyses
#' @examples
#' #if using the default sea level reconstruction
#' getintervals_time(100,10) #time cutoff of 100,000 years, for 10 time intervals
#' #if using a custom sea level reconstruction
#' sealvlfile <- read.csv("reconstruction.csv")
#' getintervals_time(100,10,sealvl=sealvlfile) #time cutoff of 100,000 years, for 10 time intervals, using a custom sea level reconstruction file.
#' @export
getintervals_time <- function(time,intervals,sealvl=bintanja_vandewal_2008) {

  #time should be in kya (i.e. 10,000 years = 10 kya)

  message("Preparing interval file... ")
  time = as.numeric(time)

  #check time value to see if it makes sense
  if (time > max(sealvl$Time)) {
    stop("Time cutoff exceeds temporal range of sea level reconstruction")
  }

  intervals = as.numeric(intervals)

  #error checking for interval value
  if (intervals > time*10) {
    stop("Number of intervals is too damn high!")
  }

  #create outputfile
  if (base::dir.exists("output")==FALSE) {
    base::dir.create("output")
  }

  #subset sealvl table to only include rows within time cutoff
  sealvl_subset <- sealvl %>%
    dplyr::filter(Time <= time) %>%
    dplyr::filter(Time > 0) %>%
    dplyr::select(Time,Sealevel_Corrected)

  timeinterval = time/intervals

  #create intervalfile dataframe
  intervalfile <- tibble::tibble(
    Interval = numeric(),
    MinDepth = numeric(),
    MaxDepth = numeric(),
    MeanDepth = numeric(),
    LowerTimeBound = numeric(),
    UpperTimeBound = numeric(),
    TimeInterval = numeric()
  )

  #fill in row 0 (for time = 0)
  intervalfile <- intervalfile %>% tibble::add_row(
    Interval = 0,
    MinDepth = 0,
    MaxDepth = 0,
    MeanDepth = 0,
    LowerTimeBound = 0,
    UpperTimeBound = 0,
    TimeInterval = 0.1
  )

  #this is the main for-loop that calculates the mean sea level for each time bin
  for (x in 1:intervals) {
    intvl <- sealvl_subset %>%
      dplyr::filter(Time > timeinterval*(x-1)) %>%
      dplyr::filter(Time <= timeinterval*x)
    intervalfile <- intervalfile %>% tibble::add_row(
      Interval = x,
      MinDepth =  min(intvl$Sealevel_Corrected),
      MaxDepth = max(intvl$Sealevel_Corrected),
      MeanDepth = base::mean(intvl$Sealevel_Corrected),
      LowerTimeBound = min(intvl$Time),
      UpperTimeBound = max(intvl$Time),
      TimeInterval = length(intvl$Time)*0.1
    )
  }

  #write output file
  write.csv(intervalfile,"output/intervals.csv")
  message("Done!")
}

#' Generate interval file, binning by sea level
#'
#' #' This function generates an interval file, which is the master file for controlling all downstream PleistoDist calculations. This file essentially simplifies Pleistocene-era sea level change into a number of discrete sea level intervals.
#'
#' @param time Maximum time cutoff in kya. PleistoDist calculates the average distance between islands over a user-defined period of time. Users will have to specify an upper time bound (in thousands of years [kya]) for their PleistoDist analysis, which by default can range from 0.1 kya to 3000 kya (i.e. 100 to 3,000,000 years ago). The lower time bound is fixed at the present day (0 kya).
#' @param intervals Number of sea level interval bins to generate between the present day and the time-cutoff. More sea level intervals will lead to more fine-grained geomorphology calculations, at the expense of increased computational time.
#' @param sealvl Optional parameter defining the sea level reconstruction file. By default, PleistoDist uses the sea level reconstruction from Bintanja & van de Wal (2008). If you use a custom sea level reconstruction file, make sure the file has a Time column (in kya) and a Sealevel_Corrected column (in metres), and is saved in CSV format.
#' @return An interval file in the automatically-generated 'output' folder. This interval file is needed for all downstream PleistoDist analyses
#' @examples
#' getintervals_sealvl(100,10) #time cutoff of 100,000 years, for 10 sea level depth intervals
#' getintervals_sealvl(100,10,sealvl="customfile.csv") #time cutoff of 100,000 years, for 10 sea level depth intervals, using a custom sea level reconstruction file
#' @export
getintervals_sealvl <- function(time,intervals,sealvl=bintanja_vandewal_2008) {

  #time should be in kya (i.e. 10,000 years = 10 kya)

  message("Preparing interval file... ")
  time = as.numeric(time)

  #check time value to see if it makes sense
  if (time > max(sealvl$Time)) {
    stop("Time cutoff exceeds temporal range of sea level reconstruction")
  }
  
  #create outputfile
  if (base::dir.exists("output")==FALSE) {
    base::dir.create("output")
  }

  intervals = as.numeric(intervals)

  #subset sealvl table to only include rows within time cutoff
  sealvl_subset <- sealvl %>%
    dplyr::filter(Time <= time) %>%
    dplyr::filter(Time > 0) %>%
    dplyr::select(Time,Sealevel_Corrected)

  rangemax = max(sealvl_subset$Sealevel_Corrected)
  rangemin = min(sealvl_subset$Sealevel_Corrected)

  depthinterval = rangemin/intervals

  #create intervalfile dataframe
  intervalfile <- tibble::tibble(
    Interval = numeric(),
    MinDepth = numeric(),
    MaxDepth = numeric(),
    MeanDepth = numeric(),
    TimeInterval = numeric()
  )

  #fill in row 0 (for time = 0)
  intervalfile <- intervalfile %>% tibble::add_row(
    Interval = 0,
    MinDepth = 0,
    MaxDepth = 0,
    MeanDepth = 0,
    TimeInterval = 0.1
  )

  for (x in 1:intervals) {
    if (x == 1) {
      intvl <- sealvl_subset %>%
        dplyr::filter(Sealevel_Corrected <= rangemax) %>%
        dplyr::filter(Sealevel_Corrected >= depthinterval)
      intervalfile <- intervalfile %>% tibble::add_row(
        Interval = x,
        MinDepth = max(intvl$Sealevel_Corrected),
        MaxDepth = min(intvl$Sealevel_Corrected),
        MeanDepth = (max(intvl$Sealevel_Corrected) + min(intvl$Sealevel_Corrected))/2,
        TimeInterval = length(intvl$Sealevel_Corrected)*0.1
      )
    }
    else {
      intvl <- sealvl_subset %>%
        dplyr::filter(Sealevel_Corrected < depthinterval*(x-1)) %>%
        dplyr::filter(Sealevel_Corrected >= depthinterval*x)
      intervalfile <- intervalfile %>% tibble::add_row(
        Interval = x,
        MinDepth = max(intvl$Sealevel_Corrected),
        MaxDepth = min(intvl$Sealevel_Corrected),
        MeanDepth = (max(intvl$Sealevel_Corrected) + min(intvl$Sealevel_Corrected))/2,
        TimeInterval = length(intvl$Sealevel_Corrected)*0.1
      )
    }
  }
  write.csv(intervalfile,"output/intervals.csv")
  message("Done!")
}
