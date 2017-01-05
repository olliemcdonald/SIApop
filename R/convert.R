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


##------------------------------------------------------------------------
#' Converts data to fish object in package fishplot for plotting
#'
#' @param time_data dataframe of time course data
#' @param timepoints timepoints to use
#' @param threshold minimum allele frequency to count a clone at
#' 
#'
#' @return convert_fishplot - returns a fish object for use with fishPlot()
#'
convert_fishplot <- function(time_data, threshold = 0.0001, timepoints = NULL){
  
  if(is.null(timepoints)) timepoints <- unique(time_data$time)

  max_cell_count <- max((time_data %>% group_by(time) %>%
                         summarize(numcells = sum(numcells)))$numcells)

  time_data <- time_data %>% filter(time %in% timepoints) %>%
    select(time, unique_id, allelefreq)
  

  
  to_keep <- (time_data %>% mutate(allelefreq = allelefreq / max_cell_count) %>%
                group_by(unique_id) %>% summarize(maxfreq = max(allelefreq)) %>% 
                filter(maxfreq > threshold))$unique_id
  
  parents <- .match_parents(to_keep)
  
  time_data <- time_data %>% filter(unique_id %in% to_keep) %>%
    mutate(allelefreq = allelefreq / max_cell_count * 100)
  
  frac.table <- as.matrix((time_data %>% spread(time, allelefreq, fill = 0))[,-1])
  
  fish = createFishObject(frac.table, parents, timepoints = timepoints,
                          col = rainbow(nrow(frac.table)),
                          clone.labels = to_keep)
  fish = layoutClones(fish)
  return(fish)
}

##------------------------------------------------------------------------
#' Converts data to Muller_df object in package ggmuller for plotting
#'
#' @param time_data dataframe of time course data
#' @param timepoints timepoints to use
#' @param threshold minimum allele frequency to count a clone at
#' @param freqplot if TRUE then plot as a frequency relative to each time point
#' @param reduce if TRUE then reduce clone names to their final ID rather than their ancestry
#' 
#'
#' @return convert_fishplot - returns a fish object for use with fishPlot()
#'
convert_ggmuller <- function(time_data, threshold = 0.001, timepoints = NULL, freqplot = FALSE, reduce = TRUE, ...){

  if(is.null(timepoints)) timepoints <- unique(time_data$time)
  time_data <- time_data %>% filter(time %in% timepoints)
  
  ##################################################################
  ### if want to call allelefreq as freq w.r.t total number of cells
  ### throughout time instead of at a single time point
  # max_cell_count <- max((time_data %>% group_by(time) %>%
  #                          summarize(numcells = sum(numcells)))$numcells)
  # to_keep <- (time_data %>%
  #               mutate(allelefreq = allelefreq / max_cell_count) %>%
  #               group_by(unique_id) %>%
  #               summarize(maxfreq = max(allelefreq))%>%
  #               filter(maxfreq >= threshold))$unique_id
  ##################################################################
 
  to_keep <- (time_data %>% group_by(time) %>%
                mutate(allelefreq = allelefreq / sum(numcells)) %>%
                ungroup %>% group_by(unique_id) %>%
                summarize(maxfreq = max(allelefreq)) %>%
                filter(maxfreq >= threshold))$unique_id
  

  edgelist <- .create_adj_matrix(to_keep, reduce = reduce)
  
  time_data <- time_data %>% select(time, unique_id, numcells)
  pop_df <- time_data %>% filter(unique_id %in% to_keep) %>%
    rename(Generation = time, Identity = unique_id, Population = numcells)
  
  if(reduce) pop_df$Identity <- unlist(lapply(pop_df$Identity, .pop, ">"))
  
  max_cell_count <- max((pop_df %>% group_by(Generation) %>%
                           summarize(numcells = sum(Population)))$numcells)
  
  pop_df <- rbind(data.frame(Generation = -1, Identity = "0", Population = 0), pop_df)
  pop_df <- pop_df %>% spread(Identity, Population, fill = 0)
  
  if(freqplot) {
    pop_df$'0' <- 0
  } else {
    pop_df$'0' <- max_cell_count - rowSums(pop_df[,-1]) + 1e-8
  }
  
  pop_df <- pop_df %>% gather(Identity, "Population", -1) %>% arrange(Generation)
  
  return(list(edgelist = edgelist, pop_df = pop_df))
}

##------------------------------------------------------------------------
#' Wrapper to convert from list containing edgelist and pop_df to Muller_df
#'
#' @param Muller_list list as output from convert_ggmuller containing edge list and pop_df
#'
#' @return create_Muller_df - returns a Muller_df tibble for Muller_plot
#'
create_Muller_df <- function(Muller_list){
  get_Muller_df(Muller_list$edgelist, Muller_list$pop_df)
}

##------------------------------------------------------------------------
#' Wrapper for ggmuller:adj_matrix_to_tree to return the edgelist output from
#' convert_ggmuller as a tree structure from ape package
#'
#' @param Muller_list list as output from convert_ggmuller containing edge list and pop_df
#'
#' @return adj_matrix_to_tree - wrapper to return the edgelist as a tree structure
#'
adj_matrix_to_tree <- function(Muller_list){
  edgelist <- Muller_list$edgelist

  ggmuller::adj_matrix_to_tree(edgelist)
}
