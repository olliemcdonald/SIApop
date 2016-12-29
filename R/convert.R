##------------------------------------------------------------------------
#' pops off string after the last instance of pattern in a vector
#'
#' @param string - character string where matches are sought.
#' @param pattern - character string to search for.
#'
#' @return string - string with final pattern removed.
#'
.pop <- function(string, pattern = ">"){
  string <- unlist(strsplit(string, pattern))
  string <- string[1:length(string)-1]
  paste(string, collapse = pattern)
}

# matches clone index to their respective
# parent clone index
##------------------------------------------------------------------------
#' matches an allele index to the respective parent index
#'
#' @param ids - vector of allele ids.
#'
#' @return .match_parents - returns a vector of parent indices where 0 indicates an ancestor
#'
.match_parents <- function(ids){
  parents <- unlist(lapply(ids, .pop, ">"))
  parents <- match(parents, ids)
  ifelse(is.na(parents), 0, parents)
}

##------------------------------------------------------------------------
#' Converts data to fish object in package fishplot for plotting
#'
#' @param time_data dataframe of time course data
#' @param timepoints timepoints to use
#' @param min_freq minimum allele frequency to count a clone at
#' 
#'
#' @return convert_fishplot - returns a fish object for use with fishPlot()
#'
convert_fishplot <- function(time_data, min_freq = 0.0001, timepoints = NULL){
  
  if(is.null(timepoints)) timepoints <- unique(time_data$time)

  max_cell_count <- max((time_data %>% group_by(time) %>%
                         summarize(numcells = sum(numcells)))$numcells)

  time_data <- time_data %>% filter(time %in% timepoints) %>%
    select(time, unique_id, allelefreq)
  

  
  to_keep <- (time_data %>% mutate(allelefreq = allelefreq / max_cell_count) %>%
                group_by(unique_id) %>% summarize(maxfreq = max(allelefreq)) %>% 
                filter(maxfreq > min_freq))$unique_id
  
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



convert_ggmuller <- function(){
}