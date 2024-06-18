### Report Generation

#' test1.
#'
#' This function takes two numeric inputs and returns their sum.
#'
#' @param x A numeric value.
#' @param y A numeric value.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add_numbers(2, 3)
#' @export
Gen_QQC_report_DDA <- function(data_path, linker_path, plex_exp,output_path = "output_report.html", isolation, CV_thresh = .4) {

    # Set up parameters for the report
    file_paths <- list(data_path = data_path, linker_path = linker_path, plex = plex_exp, isolation = isolation, CV_thresh = CV_thresh)


    # Generate the report using knitr
    rmarkdown::render(input = system.file("rmarkdown/QuantQC_DDA.Rmd", package = "QuantQC"),
                      output_format = "html_document",
                      output_file = output_path,
                      params = file_paths)

}

#' test1.
#'
#' This function takes two numeric inputs and returns their sum.
#'
#' @param x A numeric value.
#' @param y A numeric value.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add_numbers(2, 3)
#' @export
Gen_QQC_report_DIA <- function(data_path, linker_path, isolation, prot_vis_umap = NULL, output_path = "output_report.html", ChQ = .1,
                               plex_exp = 3,carrier_used = F) {

  # Set up parameters for the report
  param_vals <- list(data_path = data_path, linker_path = linker_path, isolation = isolation,
                     prot_vis_umap = prot_vis_umap, ChQ = ChQ,plex_exp = plex_exp,
                     carrier_used = carrier_used)


  # Generate the report using knitr
  rmarkdown::render(input = system.file("rmarkdown/QuantQC_DIA.Rmd", package = "QuantQC"),
                    output_format = "html_document",
                    output_file = output_path,
                    params = param_vals)

}


