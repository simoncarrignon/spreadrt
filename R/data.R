
#' All cascades from Vosoughi et al
#'
#' A dataset containing a selection of properties of the cascades analysed in Vosoughi et al 2018.
#'
#' \itemize{
#'   \item date. data of apparition of the cascade
#'   \item rumor_id. Id of the rumor associated to the cascade
#'   \item breadth. breadth of the cascade (\emph{ie} max number of parallel retweets) 
#'   \item size. size of the cascade (\emph{ie} number of unique retweets) 
#'   \item depth. depth of the cascade (\emph{ie} max number successive of unique retweets)
#'   \item veracity. veracity of the associated rumor (in  \{"TRUE","FALSE","MIXED"\})
#'   \item time_last. number of second between the apparition of the cascade and the time of the last retweet
#' }
#'
#' @docType data
#' @keywords datasets
#' @name allca
#' @usage data(allca)
#' @format A data frame with 126301 rows and 7 variables
#' @export allca 
NULL


#' True Cascades from Vosoughi et al
#'
#' Only the true cascades from in Vosoughi et al 2018.
#'
#' @docType data
#' @keywords datasets
#' @name trueca
#' @usage data(trueca)
#' @seealso \code{\link{allca}}
#' @export trueca 
NULL


#' Mixed Cascades from Vosoughi et al
#'
#' Only the mixed cascades from in Vosoughi et al 2018.
#'
#' @docType data
#' @keywords datasets
#' @name mixedca
#' @usage data(mixedca)
#' @seealso \code{\link{allca}}
#' @export mixedca 
NULL

#' False Cascades from Vosoughi et al
#'
#' Only the false cascades from in Vosoughi et al 2018.
#'
#' @docType data
#' @keywords datasets
#' @name falseca
#' @usage data(falseca)
#' @seealso \code{\link{allca}}
#' @export falseca
NULL

