#' @export
#' @importFrom utils packageVersion
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom S4Vectors DataFrame SimpleList
setClass("FacileSingleCellExperiment",
         slots=c(stuff = "list",
                 more_stuff = "list"),
         contains = "SingleCellExperiment",
         prototype = prototype(stuff = list(),
                               more_stuf = list()))

#' @export
#' @rdname facilitate
#' @method facilitate SingleCellExperiment
facilitate.SingleCellExperiment <- function(x, ...) {

}

