#' Sum-posterior-log-density
#'
#' Nimble Sampler function to record the sum of the (unnormalized) model posterior log-density, on every iteration of an MCMC. This is based on the vignette https://danielturek.github.io/public/sumPostLogDens/sumPostLogDens.html
#'
#' @param model .
#' @param mvSaved .
#' @param target .
#' @param control .
#'
#' @importFrom nimble nimbleFunction
#'
#' @export
sumLogPostDens <- nimbleFunction(
  name = 'sumLogPostDens',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    stochNodes <- setdiff(model$getNodeNames(stochOnly = TRUE), target)
  },
  run = function() {
    model[[target]] <<- model$getLogProb(stochNodes)
  },
  methods = list( reset = function() {} )
)