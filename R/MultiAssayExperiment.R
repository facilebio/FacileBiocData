#' @export
#' @importFrom utils packageVersion
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @importClassesFrom S4Vectors DataFrame SimpleList
setClass("FacileMultiAssayExperiment",
         slots=c(stuff = "list",
                 more_stuff = "list"),
         contains = "MultiAssayExperiment",
         prototype = prototype(stuff = list(),
                               more_stuf = list()))


#' @export
#' @rdname facilitate
#' @method facilitate MultiAssayExperiment
facilitate.MultiAssayExperiment <- function(x, ...) {

}
