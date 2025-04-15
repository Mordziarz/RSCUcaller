#' @title Sample Grouping Metadata
#' 
#' @description 
#' Demonstration dataset showing required structure for group comparisons in RSCU analysis.
#' Used with functions like \code{\link{boxplot_between_groups}} and \code{\link{PCA_RSCU}}.
#'
#' @format A data frame with 8 observations and 2 variables:
#' \describe{
#'   \item{Species}{Character vector matching sample names from \code{get_RSCU} output.
#'                 Must follow "number_Species" format (e.g. "1_Homo_sapiens")}
#'   \item{group}{Character vector specifying experimental groups or categories}
#' }
#'
#' @usage data(grouping_table)
#'
#' @source Synthetic data created for demonstration purposes. Not from actual research data.
#'
#' @examples
#' # Load dataset
#' data(grouping_table)
#' 
#' # View structure
#' str(grouping_table)
#' 
#' # Typical usage with analysis functions
#' \dontrun{
#' boxplot_between_groups(get_RSCU_out = rscu_data,
#'                        grouping_table = grouping_table)
#' }
#' 
#' @seealso
#' Related functions:
#' \code{\link{boxplot_between_groups}}, \code{\link{PCA_RSCU}}, \code{\link{neutrality_pr2}}
"grouping_table"