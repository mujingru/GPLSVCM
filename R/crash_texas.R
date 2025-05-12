#' Crash data in Texas
#'
#' A dataset containing the number of car crashes within each census tract in
#' Texas of year 2107 and other variables of 4771 locations.
#'
#' @format A data frame with 4771 rows and 7 variables:
#' \describe{
#'   \item{count}{off-roadway crash frequencies}
#'   \item{vmt}{log of vehicle miles traveled}
#'   \item{pop}{log of total population}
#'   \item{old}{proportion of people age 65 and older}
#'   \item{hispanics}{proportion of Hispanics}
#'   \item{lon}{longitude of a location}
#'   \item{lat}{latitude of a location}
#' }
#'
#' @docType data
#'
#' @usage data(Crash_Texas)
#'
#' @keywords datasets
#'
#' @examples
#' data(Crash_Texas)
#' count <- Crash_Texas$count
#' hist(count)
#' summary(count)

"Crash_Texas"
