#' Make a trial
#'
#' @param trial_num integer, the trial number
#' @param event_names character vector, names from event_vectors matrix for the events that appear in the trial
#' @param onsets integer vector, onset times for each event listed in event_names
#' @param durations integer vector, durations for each event listed in event_names
#'
#' @return data.frame
#' @export
#'
#' @examples
#' make_trial(1,event_names=c("A","B"),c(1,100),duration=50)

make_trial <- function(trial_num,
                       event_names,
                       onsets,
                       durations) {

  trial <- data.frame(
    trial_num = trial_num,
    event = event_names,
    onset = onsets,
    duration = durations,
    offset = onsets+durations
  )

  return(trial)
}


#' Make AB timeline
#'
#' @param num_trials integer, the number of trials to add to timeline
#' @param events character vector, defaults to A and B events
#' @param A_onsets integer vector, the onset times for each A event on each trial
#' @param A_durations integer vector, the duration times for each A event on each trial
#' @param B_onsets integer vector, the onset times for each B event on each trial
#' @param B_durations integer vector, the duration times for each B event on each trial
#'
#' @return data.frame
#' @export
#'
#' @examples
#' make_AB_timeline(
#'  num_trials = 5,
#'  events = c('A', 'B'),
#'  A_onsets = 1,
#'  A_durations = 50,
#'  B_onsets = round(rnorm(5,100,10)),
#'  B_durations = 50
#' )

make_AB_timeline <- function(num_trials = 5,
                             events = c('A', 'B'),
                             A_onsets = 1,
                             A_durations = 50,
                             B_onsets = 50,
                             B_durations = 50) {
  timeline <- data.frame()

  for (i in 1:num_trials) {

    A_on <- A_onsets[1]
    if (length(A_onsets) > 1) A_on <- A_onsets[i]
    B_on <- B_onsets[1]
    if (length(B_onsets) > 1) B_on <- B_onsets[i]

    new_trial <- make_trial(
      i,
      event_names = events,
      onsets = c(A_on, B_on),
      durations = c(A_durations, B_durations)
    )

    timeline <- rbind(timeline, new_trial)
  }

  return(timeline)

}

#' Make timeline of events
#'
#' @param num_trials integer, the number of trials
#' @param events character vector, names of each event
#' @param onsets list, containing an integer, or vector, defining onsets for each event
#' @param durations list, containing an integer, or vector, defining durations for each event
#'
#' @return dataframe
#' @export
#'
#' @examples
#' make_event_timeline(
#'  num_trials = 5,
#'  events = c('A', 'B', 'C'),
#'  onsets = list(1,round(rnorm(5,50,10)),1),
#'  durations = list(25,25,199)
#'  )
make_event_timeline <- function(num_trials = 5,
                                events = c('A','B','C'),
                                onsets = list(1,round(rnorm(5,50,10)),1),
                                durations = list(25,25,199)
){
  timeline <- data.frame()

  for(i in 1:num_trials){

    # could functionalize this later
    process_onsets <- onsets
    for(o in 1:length(process_onsets)){
      if(length(process_onsets[[o]]) == 1) {
        process_onsets[[o]] <- rep(process_onsets[[o]],num_trials)
      }
    }
    process_onsets <- matrix(unlist(process_onsets),nrow=num_trials, byrow=FALSE)

    process_durations <- durations
    for(d in 1:length(process_durations)){
      if(length(process_durations[[d]]) == 1) {
        process_durations[[d]] <- rep(process_durations[[d]],num_trials)
      }
    }
    process_durations <- matrix(unlist(process_durations),nrow=num_trials, byrow=FALSE)

    new_trial <- make_trial(i,
                            event_names=events,
                            onsets=process_onsets[i,],
                            durations=process_durations[i,])

    timeline <- rbind(timeline,new_trial)
  }

  return(timeline)

}


#' Make Temporal Vectors
#'
#' @param num_timesteps integer, the number of timesteps within a trial
#' @param overlap integer, controls overlap between vectors
#'
#' @return matrix, each row is a temporal vector
#' @export
#'
#' @examples
#' make_temporal_vectors(10,2)
make_temporal_vectors <- function(num_timesteps = 200, overlap = 2){

  temporal_matrix <- matrix(1,
                            nrow = num_timesteps,
                            ncol = num_timesteps)

  lower <- -1*overlap
  upper <- overlap
  delta <- row(temporal_matrix) - col(temporal_matrix)
  temporal_matrix[delta < lower | delta > upper] <- 0

  return(temporal_matrix)

}


#' Make Event Vectors
#'
#' @param event_names character vector, names for each event
#' @param vector_length integer, number of 1s for each vector
#' @param type string, 'field' sets all values in each event vector to 1, which is used for the concatenation approach, 'riv' creates sparse random index vectors used in the superposition approach.
#' @param sparsity numeric, proportion of non-zero elements (1s, or -1s) in a random index vector
#'
#' @return matrix
#' @export
#'
#' @examples
#' make_event_vectors()
make_event_vectors <- function(event_names = LETTERS[1:10],
                               vector_length = 10,
                               type = 'field',
                               sparsity = .1) {
  if(type == 'field'){
    event_vectors <- matrix(1,
                            ncol = vector_length,
                            nrow = length(event_names))

    row.names(event_vectors) <- event_names
  }

  if(type == 'riv'){

    index_vector <- c(rep(1,round((sparsity/2)*vector_length)),
                      rep(-1,round((sparsity/2)*vector_length)),
                      rep(0,vector_length-(round(sparsity*vector_length)))
    )

    event_vectors <-t(replicate(length(event_names),
                                sample(index_vector)))

    row.names(event_vectors) <- event_names
  }

  if(type == 'binomial') {

    index_vector <- c(rep(1,round(vector_length/2)),
                      rep(-1,round(vector_length/2))
    )

    event_vectors <-t(replicate(length(event_names),
                                sample(index_vector)))

    row.names(event_vectors) <- event_names

  }

  return(event_vectors)

}


#' Convert timeline to lists of environment vectors
#'
#' @param timeline dataframe, created by make_timeline()
#' @param event_vectors matrix, created by make_event_vectors()
#' @param temporal_vectors matrix, created by make_temporal_vectors()
#' @param timesteps integer, number of timesteps per trial
#'
#' @return list, containing a environment vector matrix for each trial.
#' @export
#'
#' @examples
timeline_to_vector <- function(timeline,
                               event_vectors,
                               temporal_vectors,
                               timesteps = 200) {
  timeline_list <- list()

  # iterate through each trial
  for(t in 1:max(timeline$trial_num)){

    # get events for current trial
    current_trial <- timeline[timeline$trial_num == t,]

    # create event vector matrix
    all_events_matrix <- matrix(0,
                                nrow = timesteps,
                                ncol = dim(event_vectors)[1]*dim(event_vectors)[2])

    # loop through each event name in trial
    for (e in 1:length(current_trial$event)){

      # create blank matrix for event
      e_matrix <- matrix(0,
                         nrow = timesteps,
                         ncol = length(event_vectors[current_trial[e,'event'],]))

      # add event vector for specified durations in matrix
      e_matrix[(current_trial[e,'onset']:current_trial[e,'offset']),] <- event_vectors[current_trial[e,'event'],]

      #append matrix for each event to whole trial matrix
      #all_events_matrix <-  cbind(all_events_matrix,e_matrix)

      #insert event matrix to it's field position
      e_name <- current_trial[e,'event']
      e_row <- which(row.names(event_vectors) %in% e_name)
      e_length <- length(event_vectors[e_row,])
      first_ind <- ((e_row-1)*e_length)+1
      last_ind <- first_ind + (e_length-1)
      all_events_matrix[ ,first_ind:last_ind] <- e_matrix

    }

    # all_events_matrix <- all_events_matrix[,-1] # delete column of NAs

    # append temporal vectors
    all_events_matrix <- cbind(all_events_matrix,temporal_vectors)

    # save to list
    timeline_list[[t]] <- all_events_matrix
  }

  return(timeline_list)
}

#' Convert timeline to lists of riv environment vectors
#'
#' @param timeline dataframe, created by make_timeline()
#' @param event_vectors matrix, created by make_event_vectors()
#' @param temporal_vectors matrix, created by make_temporal_vectors()
#' @param timesteps integer, number of timesteps per trial
#'
#' @return list, containing a environment vector matrix for each trial.
#' @export
#'
#' @examples
timeline_to_riv_vector <- function(timeline,
                               event_vectors,
                               temporal_vectors,
                               timesteps = 200) {
  timeline_list <- list()

  # iterate through each trial
  for(t in 1:max(timeline$trial_num)){

    # get events for current trial
    current_trial <- timeline[timeline$trial_num == t,]

    # create blank matrix for event
    e_matrix <- temporal_vectors

    # loop through each event name in trial
    for (e in 1:length(current_trial$event)){

      # add event vector for specified durations in matrix

      the_event <- t(replicate(length((current_trial[e,'onset']:current_trial[e,'offset'])),event_vectors[current_trial[e,'event'],]))

      the_event_by_time <- t(replicate(length((current_trial[e,'onset']:current_trial[e,'offset'])),event_vectors[current_trial[e,'event'],])) * temporal_vectors[(current_trial[e,'onset']:current_trial[e,'offset']),]

      e_matrix[(current_trial[e,'onset']:current_trial[e,'offset']),] <- e_matrix[(current_trial[e,'onset']:current_trial[e,'offset']),] + the_event + the_event_by_time

    }

    # save to list
    timeline_list[[t]] <- e_matrix
  }

  return(timeline_list)
}

#' Convert timeline to environment vectors with noise
#'
#' @param timeline dataframe, created by make_timeline()
#' @param event_vectors matrix, created by make_event_vectors()
#' @param temporal_vectors matrix, created by make_temporal_vectors()
#' @param timesteps integer, number of timesteps per trial
#'
#' @return list, containing a environment vector matrix for each trial.
#' @export
#'
#' @examples
timeline_to_vector_noise <- function(timeline,
                                   event_vectors,
                                   temporal_vectors,
                                   timesteps = 200,
                                   noise = TRUE,
                                   L = .1) {
  timeline_list <- list()

  # iterate through each trial
  for(t in 1:max(timeline$trial_num)){

    # get events for current trial
    current_trial <- timeline[timeline$trial_num == t,]

    # get vectors
    t_vectors <- temporal_vectors
    e_vectors <- event_vectors[current_trial[,'event'],]

    te_vectors <- matrix(0,nrow=1, ncol=dim(temporal_vectors)[2])
    for (e in 1:dim(e_vectors)[1]){
      e_onset <- current_trial[e,'onset']
      e_offset <- current_trial[e,'offset']
      te_vectors <- rbind(te_vectors,
                          vec_x_mat(e_vectors[e,],
                                    t_vectors[e_onset:e_offset,]))
    }
    te_vectors <- te_vectors[-1,]

    e_matrix <- rbind(t_vectors,
                      e_vectors,
                      te_vectors)

    if(noise == TRUE){
      e_matrix <- add_noise_to_matrix(e_matrix,
                                      L,
                                      merge=TRUE)
    }

    # save to list
    timeline_list[[t]] <- e_matrix
  }

  return(timeline_list)
}

#' Normalize a vector by max absolute value
#'
#' @param x numeric vector, a vector of numbers
#'
#' @return numeric vector, the vector divided by the largest absolute value in the vector
#' @examples
#' a <- c(1,2,3)
#' normalize_vector(a)
#'
#' @export
normalize_vector <- function (x) {return(x/abs(max(x)))}

vec_x_mat <- function(v,m){
  return(t(t(m)*v))
}

#' Add noise or create noise matrix
#'
#' @param x matrix
#' @param learning_rate numeric from 0 to 1, determining probability that each element in the row is sampled perfectly
#' @param merge logical default is TRUE and return original matrix multiplied by noise matrix. FALSE will return only the noise matrix
#'
#' @return
#' @export
#'
#' @examples
#'

add_noise_to_matrix <- function(x, learning_rate= .9, merge= TRUE){
  dims <- dim(x)
  noise <- matrix(sample(c(1,0),
                         size = dims[1] * dims[2],
                         replace = TRUE,
                         prob=c(learning_rate,1-learning_rate)),
                  nrow = dims[1],
                  ncol = dims[2])
  if (merge == FALSE) return(noise)
  if (merge == TRUE) return(x*noise)
}
