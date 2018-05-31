#' @export
#' @importFrom utils packageVersion
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom S4Vectors DataFrame SimpleList
setClass("FacileSummarizedExperiment",
         slots=c(stuff = "list",
                 more_stuff = "list"),
         contains = "SummarizedExperiment",
         prototype = prototype(stuff = list(),
                               more_stuf = list()))

#' @export
#' @rdname facilitate
#' @method facilitate SummarizedExperiment
facilitate.SummarizedExperiment <- function(x, ...) {

}
