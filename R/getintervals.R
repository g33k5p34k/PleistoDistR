
#' Generate interval file, binning by time
#'
#' This function generates an interval file, which is the master file for controlling all downstream PleistoDist calculations. This file essentially simplifies Pleistocene-era sea level change into a number of discrete time intervals.
#'
#' @param time Maximum time cutoff in kya. PleistoDist calculates the average distance between islands over a user-defined period of time. Users will have to specify an upper time bound (in thousands of years (kya)) for their PleistoDist analysis, which by default ranges from 0.1 kya to 3000 kya (i.e. 100 to 3,000,000 years ago). The lower time bound is fixed at the present day (0 kya).
#' @param intervals Number of time interval bins to generate between the present day and the time-cutoff. More time intervals will lead to more fine-grained geomorphology calculations, at the expense of increased computational time. Under default conditions, the number of intervals should not exceed time-cutoff*10.
#' @param outdir Output directory for the interval file. If the specified output directory doesn't already exist, PleistoDist will create the output directory.
#' @param sealvl Optional parameter defining the sea level reconstruction file. By default, PleistoDist uses the sea level reconstruction from Bintanja & van de Wal (2008). If you use a custom sea level reconstruction file, make sure the file has a Time column (in kya) and a Sealevel_Corrected column (in metres), and is saved in CSV format.
#' @return An interval file in the automatically-generated 'output' folder. This interval file is needed for all downstream PleistoDist analyses
#' @examples
#' #create temp directory
#' path <- file.path(tempdir())
#' #generate 5 time-binned intervals, with a cutoff of 5kya,
#' #using default sea-level reconstruction by Bintanja & Van de Wal (2008)
#' getintervals_time(time=5,intervals=5,outdir=path)
#' #if using sea-level reconstruction by Spratt & Lisecki (2016)
#' getintervals_time(time=5,intervals=5,outdir=path,sealvl=PleistoDist::spratt_lisecki_2016)
#' @export
getintervals_time <- function(time,intervals,outdir,sealvl=PleistoDist::bintanja_vandewal_2008) {

  #create output folder
  if (base::dir.exists(outdir)==FALSE) {
    base::dir.create(outdir)
  }

  #time should be in kya (i.e. 10,000 years = 10 kya)
  message("Preparing interval file... ")
  time = as.numeric(time)

  #check time value to see if it makes sense
  if (time > max(sealvl$Time)) {
    stop("Time cutoff exceeds temporal range of sea level reconstruction")
  }

  intervals = as.numeric(intervals)

  #calculate the time break between each row in the sea level reconstruction file
  unittime = base::round(base::max(sealvl$Time)/base::length(sealvl$Time),digits=2)

  #error checking for interval value
  if (intervals > base::nrow(dplyr::filter(sealvl, Time <= time))) {
    stop("The number of intervals is too damn high!")
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
    TimeInterval = unittime
  )

  #this is the main for-loop that calculates the mean sea level for each time bin
  for (x in 1:intervals) {
    intvl <- sealvl_subset %>%
      dplyr::filter(Time > timeinterval*(x-1)) %>%
      dplyr::filter(Time <= timeinterval*x)
    intervalfile <- intervalfile %>% tibble::add_row(
      Interval = x,
      MinDepth =  base::max(intvl$Sealevel_Corrected),
      MaxDepth = base::min(intvl$Sealevel_Corrected),
      MeanDepth = base::mean(intvl$Sealevel_Corrected),
      LowerTimeBound = base::min(intvl$Time),
      UpperTimeBound = base::max(intvl$Time),
      TimeInterval = length(intvl$Time)*unittime
    )
  }

  #write output file
  utils::write.csv(intervalfile,base::paste0(outdir,"/intervals.csv"))
  message("Done!")
}

#' Generate interval file, binning by sea level
#'
#' #' This function generates an interval file, which is the master file for controlling all downstream PleistoDist calculations. This file essentially simplifies Pleistocene-era sea level change into a number of discrete sea level intervals.
#'
#' @param time Maximum time cutoff in kya. PleistoDist calculates the average distance between islands over a user-defined period of time. Users will have to specify an upper time bound (in thousands of years (kya)) for their PleistoDist analysis, which by default can range from 0.1 kya to 3000 kya (i.e. 100 to 3,000,000 years ago). The lower time bound is fixed at the present day (0 kya).
#' @param intervals Number of sea level interval bins to generate between the present day and the time-cutoff. More sea level intervals will lead to more fine-grained geomorphology calculations, at the expense of increased computational time.
#' @param outdir Output directory for the interval file. If the specified output directory doesn't already exist, PleistoDist will create the output directory.
#' @param sealvl Optional parameter defining the sea level reconstruction file. By default, PleistoDist uses the sea level reconstruction from Bintanja & van de Wal (2008). If you use a custom sea level reconstruction file, make sure the file has a Time column (in kya) and a Sealevel_Corrected column (in metres), and is saved in CSV format.
#' @return An interval file in the automatically-generated 'output' folder. This interval file is needed for all downstream PleistoDist analyses
#' @examples
#' #create temp directory
#' path <- file.path(tempdir())
#' #generate 5 sea level-binned intervals, with a cutoff of 20kya,
#' #using default sea-level reconstruction by Bintanja & Van de Wal (2008)
#' getintervals_sealvl(time=20,intervals=5,outdir=path)
#' #if using sea-level reconstruction by Spratt & Lisecki (2016)
#' getintervals_sealvl(time=20,intervals=5,outdir=path,sealvl=PleistoDist::spratt_lisecki_2016)
#' @export
getintervals_sealvl <- function(time,intervals,outdir,sealvl=PleistoDist::bintanja_vandewal_2008) {

  #time should be in kya (i.e. 10,000 years = 10 kya)

  message("Preparing interval file... ")
  time = as.numeric(time)

  #check time value to see if it makes sense
  if (time > max(sealvl$Time)) {
    stop("Time cutoff exceeds temporal range of sea level reconstruction")
  }

  intervals = as.numeric(intervals)

  #calculate the time break between each row in the sea level reconstruction file
  unittime = base::round(base::max(sealvl$Time)/base::length(sealvl$Time),digits=2)

  #create outputfile
  if (base::dir.exists(outdir)==FALSE) {
    base::dir.create(outdir)
  }

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
    TimeInterval = unittime
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
        TimeInterval = length(intvl$Sealevel_Corrected)*unittime
      )
    }
    else {
      intvl <- sealvl_subset %>%
        dplyr::filter(Sealevel_Corrected < depthinterval*(x-1)) %>%
        dplyr::filter(Sealevel_Corrected >= depthinterval*x)
      intervalfile <- intervalfile %>% tibble::add_row(
        Interval = x,
        MinDepth = depthinterval*x,
        MaxDepth = depthinterval*(x-1),
        MeanDepth = (depthinterval*x + depthinterval*(x-1))/2,
        TimeInterval = length(intvl$Sealevel_Corrected)*unittime
      )
    }
  }
  utils::write.csv(intervalfile,base::paste0(outdir,"/intervals.csv"))
  message("Done!")
}
