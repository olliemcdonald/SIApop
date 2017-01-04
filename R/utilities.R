### UTILITY FUNCTIONS
##------------------------------------------------------------------------
#' pops off last instance of pattern in a string
#'
#' @param string - character string where matches are sought.
#' @param pattern - character string to search for.
#'
#' @return .pop_off - string with everything after final pattern removed.
#'
.pop_off <- function(string, pattern = ">"){
  string <- unlist(strsplit(string, pattern))
  string <- string[1:length(string)-1]
  paste(string, collapse = pattern)
}

##------------------------------------------------------------------------
#' pops off last element of a string after final pattern in a vector
#'
#' @param string - character string where matches are sought.
#' @param pattern - character string to search for.
#'
#' @return .pop - string of everything after final pattern.
#'
.pop <- function(string, pattern = ">"){
  string <- unlist(strsplit(string, pattern))
  string[length(string)]
}

##------------------------------------------------------------------------
#' matches an allele index to the respective parent index
#'
#' @param ids - vector of allele ids.
#'
#' @return .match_parents - returns a vector of parent indices where 0 indicates an ancestor
#'
.match_parents <- function(ids){
  parents <- unlist(lapply(ids, .pop_off, ">"))
  parents <- match(parents, ids)
  ifelse(is.na(parents), 0, parents)
}

##------------------------------------------------------------------------
# replaces an element of a vector with a particular pattern with the
# replacement given
#'
#' @param vec - character vector
#' @param pattern - character element to search for
#' @param replacement - character to replace pattern with
#'
#' @return .replace - returns a vector with the replaced elements
#'
.replace <- function(vec, pattern, replacement) {
  vec[vec == pattern] <- replacement
  vec
}

##------------------------------------------------------------------------
#' for a list of clone ids each separated by ">" a specific allele can be
#' searched for and replaced by splitting the elements of a vector,
#' replacing individual elements, and collapsing the list back to a vector
#'
#' @param vec - character vector of ids
#' @param pattern - character element to search for
#' @param replacement - character to replace pattern with
#'
#' @return .split_replace_collapse - returns a vector with the replaced ids
#'
.split_replace_collapse <- function(vec, pattern, replacement) {
  vec <- strsplit(vec, ">")
  vec <- lapply(vec, .replace, pattern, replacement)
  unlist(lapply(vec, paste, collapse = ">"))
}

##------------------------------------------------------------------------
#' Converts data to Muller_df object in package ggmuller for plotting
#'
#' @param id_list list of clone ids
#' @param reduce if true we label a clone by its number rather than the list of ancestors
#' 
#'
#' @return create_adj_matrix - returns an edge list for the population
#'
.create_adj_matrix <- function(id_list, reduce = TRUE){
  # "0" is reserved for the root - rename any 0 allele to x0x
  id_list <- .split_replace_collapse(id_list, "0", "x0x")
  
  edgelist <- data.frame(Parent = unlist(lapply(id_list, .pop_off, ">")),
                         Identity = id_list, stringsAsFactors = F)
  edgelist$Parent[edgelist$Parent == ""] <- "0"
  
  if(reduce){
    edgelist$Parent <- unlist(lapply(edgelist$Parent, .pop, ">"))
    edgelist$Identity <- unlist(lapply(edgelist$Identity, .pop, ">"))
  }
  edgelist
  
}
