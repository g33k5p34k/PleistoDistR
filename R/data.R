#' Default Sea Level Reconstruction Dataset (Bintanja & van de Wal, 2008)
#'
#' Sea level reconstruction file from Bintanja & van de Wal (2008)
#'
#' @format ## `bintanja_vandewal_2008`
#' A data frame with 30,002 rows and 3 columns:
#' \describe{
#'   \item{Time}{Time in the past (kya)}
#'   \item{Sealevel_Uncorrected}{Raw sea level estimate from original dataset
#'    (metres)}
#'   \item{Sealevel_Corrected}{Sea level estimate, with present day set to 0
#'   (metres)}
#' }
#' @source <https://www.nature.com/articles/nature07158>
"bintanja_vandewal_2008"

#' Alternative Sea Level Reconstruction Dataset (Spratt & Lisecki, 2016)
#'
#' Sea level reconstruction file from Spratt & Lisecki (2016)
#'
#' @format ## `spratt_lisecki_2016`
#' A data frame with 800 rows and 2 columns:
#' \describe{
#'   \item{Time}{Time in the past (kya)}
#'   \item{Sealevel_Corrected}{Sea level estimate, with present day set to 0
#'   (metres)}
#' }
#' @source <https://cp.copernicus.org/articles/12/1079/2016/>
"spratt_lisecki_2016"
