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
#'
#' @return matrix
#' @export
#'
#' @examples
#' make_event_vectors()
make_event_vectors <- function(event_names = LETTERS[1:10],
                               vector_length = 10) {
  event_vectors <- matrix(1,
                          ncol = vector_length,
                          nrow = length(event_names))

  row.names(event_vectors) <- event_names

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
    all_events_matrix <- matrix(nrow = timesteps)
    for (e in 1:length(current_trial$event)){
      e_matrix <- matrix(0,
                         nrow = timesteps,
                         ncol = length(event_vectors[current_trial[e,'event'],]))

      e_matrix[(current_trial[e,'onset']:current_trial[e,'offset']),] <- event_vectors[current_trial[e,'event'],]

      #append matrix for each event to whole trial matrix
      all_events_matrix <-  cbind(all_events_matrix,e_matrix)
    }
    all_events_matrix <- all_events_matrix[,-1] # delete column of NAs

    # append temporal vectors
    all_events_matrix <- cbind(all_events_matrix,temporal_vectors)

    # save to list
    timeline_list[[t]] <- all_events_matrix
  }

  return(timeline_list)
}
