# R/link_groups.R
# linking across groups, approach chooses chain, joint, or
# restricted joint, variant the simultaneous or pairwise form

link_groups <- function(x, approach = c("chain", "joint",
                                        "joint_restricted"),
                        method = NULL,
                        variant = c("simultaneous", "pairwise"),
                        ref = NULL, ...) {
  approach <- match.arg(approach)
  variant <- match.arg(variant)
  if (approach == "chain") {
    if (is.null(method)) method <- "mgm"
    if (!method %in% c("mgm", "mm", "haberman", "haebara", "sl")) {
      stop("`method` must be one of \"mgm\", \"mm\", \"haberman\", ",
           "\"haebara\", \"sl\".", call. = FALSE)
    }
    if (variant == "pairwise") {
      stop("variant = \"pairwise\" needs approach = \"joint\" or ",
           "\"joint_restricted\".", call. = FALSE)
    }
    return(link_chain(x, method = method, ref = ref, ...))
  }
  if (is.null(method)) method <- "haberman"
  if (!method %in% c("haberman", "haebara", "sl")) {
    stop("`method` must be one of \"haberman\", \"haebara\", \"sl\" ",
         "for approach = \"", approach, "\".", call. = FALSE)
  }
  if (variant == "pairwise") {
    method <- switch(method, haberman = "phl", haebara = "haebara_pw",
                     sl = "sl_pw")
  }
  link_joint(x, method = method,
             restricted = (approach == "joint_restricted"),
             ref = ref, ...)
}
