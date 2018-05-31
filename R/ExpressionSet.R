#' @export
#' @importFrom utils packageVersion
#' @importClassesFrom Biobase ExpressionSet
setClass("FacileExpressionSet",
         slots=c(stuff = "list",
                 more_stuff = "list"),
         contains = "ExpressionSet",
         prototype = prototype(stuff = list(),
                               more_stuf = list()))


#' @export
#' @rdname facilitate
#' @method facilitate ExpressionSet
facilitate.ExpressionSet <- function(x, ...) {

}
