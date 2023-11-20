#'Align mass spectra
#'
#'GCxGC-MS data of mixture of standards including: linear hydrocarbons (nC7-nC40), Grob mixture, Lime oil
#'measured on LECO BT machine with split/splitless injection sampled on red SPME fiber automatically procesed
#'with no further data curation
#'
#'
#' @format matrix with 705 rows and 502 columns where in columns are mass to charge ratios from 1-502
#' and in rows are individual detected peaks
#' \describe{
#'  \item{m.z1:m.z502}{mass to charge ratios from 1:502}
#' }
#' @source home measured dataset
#' @examples data(Tab_align_ms)   # Lazy loading.
"Tab_align_ms"
