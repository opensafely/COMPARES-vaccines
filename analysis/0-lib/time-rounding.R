# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Core functions for rounding event times and event counts in survival data
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# taken from 
# https://github.com/opensafely-actions/kaplan-meier-function/blob/main/analysis/time-rounding.R

library("tidyverse")

times2counts <- function(times, max_time=max(times)){
  # converts a (person-level) vector of event times into a table of [time, event_count] pairs
  dplyr::count(tibble(time=times), time) |>
    tidyr::complete(
      time = seq_len(max_time),
      fill = list(n = 0L)
    )
}

counts2times <- function(time, count){
  # converts a [time, event_count] pair into a vector of event times
  rep(time, count)
}

# main rounding functions -----

# In summary: take a vector of event times and round them to a given precision.
# TODO: take a survival table of [time, number_of_events_per_period] and return rounded times
# this essentially involves replacing `time` argument with outputs of `rle(time)`
# and outputting counts across whole time period
# TODO: create test that ensures precision 1 (if integer time) or 0 (if fully continuous time) returns the same values
round_event_times <- function(time_to_event, min_count, origin=0L, group.method="conjoin", value.method="ceiling", min_increment=1L){
  
  # Purpose:
  # For a set of event times (or any ordinal / ordered categorical variable that can be encoded as 1,2,3,...), put event times into groups such that no group contains fewer than `min_count` events.
  # This is a way to satisfy "no cells less than `min_count`" criteria when reporting summaries of survival data,
  # such as those from Kaplan-Meier estimators and other cumulative incidence type estimates / quantities.
  
  # There are a few different methods to do this deterministically.
  # All methods considered here DO NOT split event times that occur at the same time into different groups (except even spacing? TODO: check this)
  # Otherwise it's vulnerable to data-ordering (and solutions become non-deterministic if not accounting for this ordering)
  
  # Grouping methods:
  # "forward"     - Moving from time zero forwards, set the next group boundary to be an (existing) event time such that the number of events less than or equal to that time is at least `min_count`.
  # "backward"    - Same as forwards, but start at the final event time and move backwards. TODO: not yet implemented
  # "conjoin"     -   * Initialise group boundaries at all unique event times;
  #                   * merge the group with minimum count into the adjacent group with the smallest count;
  #                   * repeat until there are no groups with count < `min_count`
  #                 In other words, merge smallest group to the smallest adjacent group until there are no small groups left.
  #                 Ties are dealt with deterministically, by choosing the earliest group.
  # "cumulative"  -   * Calculate the cumulative incidence of event times;
  #                   * round cumulative incidence to `min_count`;
  #                   * return to event count scale (eg with `diff()`).
  #               - Note this does not result in a total count that is equal to the original total count, as it rounds it up to `min_count`.
  
  
  # Once groups have been set, the actual modified event times can be chosen in different ways:
  
  # Output value options:
  # floor   - round all values to lower boundary of the group
  # ceiling - round all values to upper boundary of the group
  # mid     - round all values to midpoint of group
  # spaced  - assume equal spacing of event times within each group from lower to upper bound (similar to linear interpolation between group boundaries)
  #           usually optionally apply rounding to these spaced values (to ensure they're integers, for example)
  
  # The floor, ceiling, and mid methods result in "steppy" outputs that suggest event counts are `min_count` times fewer than the actual counts.
  # This might look unnatural to those who are familiar with eg KM curves.
  # The spacing method look more natural.
  # However (??), this method also smooths out steppiness even if the counts are non-disclosive (i think? but check). TODO: Is there a way to deal with this?
  
  
  ## argument tests
  stopifnot("origin should be strictly less than minimum observation time_to_event" = origin< min(time_to_event))
  stopifnot("time_to_event must be numeric" = is.numeric(time_to_event))
  stopifnot("min_count must be numeric" = is.numeric(min_count))
  stopifnot("min_increment must be numeric" = is.numeric(min_increment))
  stopifnot("group.method should be one of conjoin, forward, or cumulative" = group.method %in% c("conjoin", "forward", "cumulative"))
  stopifnot("value.method should be one of ceiling, floor, mid, or spaced" = value.method %in% c("ceiling", "floor", "mid", "spaced"))
  
  # order event times
  order_tte <- order(time_to_event)
  sort_tte <- sort(time_to_event)
  
  # if number of event times is less than `min_count`, then avoid rounding methods completely and output single group -- TODO: THIS WILL NEED TO BE REDACTED!
  if (length(time_to_event) < min_count){
    cuts <- c(origin, max(time_to_event))
    counts <- length(time_to_event)
  } else {
    
    #run length encoding
    rleobj <- rle(sort_tte)
    
    # all times where at least one event occurred, in order
    event_times <- rleobj$values
    
    # number of events at each event time (excluding zero-counts)
    event_n <- rleobj$lengths
    
    
    # define new event time boundaries to ensure at least 6 events per group
    switch(
      group.method,
      
      # start by defining one group per event time
      # the size of the group is defined by the number of events at that time
      # merge the smallest group with the smallest of the two adjacent groups (if not a boundary)
      # if the smallest group is tied, take the first (lowest in time) group
      # stop when all groups are at least `min_count` in size
      "conjoin" = {
        
        # initialise group boundaries and group counts
        # cuts are closed on the right, i.e., (a,b]
        cuts <- c(origin, event_times)
        counts <- event_n
        
        while(any(counts<min_count)){
          
          # find smallest group and its size
          minindex <- which.min(counts)
          len <- length(counts)
          
          # merge group with smallest adjacent group (if equal select left/lower group)
          if(c(Inf,counts)[minindex] <= c(counts, Inf)[minindex+1]) {
            cuts <- cuts[-(minindex)]
          } else {
            cuts <- cuts[-(minindex+1)]
          }
          
          # recalculate group sizes
          counts <- c(as.integer(table(cut(sort_tte, cuts))))
        }
        print(cuts)
      },
      
      # create groups by moving forward from time-zero.
      # a group is created when at least `min_count` events are in the group.
      # then start again at the next even time
      "forward" = {
        
        # one possibly more efficient way to derive the groups, but not using for now
        # event_cml <- cumsum(event_n)
        # cml_floor <- floor_any(event_cml, min_count)
        # cml_ceiling <- lag(cml_floor, 1, 0L) + min_count
        
        # function to define partial cumulative sum which restarts each time the sum exceeds rounding width
        # when the sum restarts, this defines the cut points for the rounded groups
        partialsum <- accumulate(
          event_n,
          function(a,b){
            sumab = a+b
            if(sumab>=min_count) 0 else sumab
          }
        )
        
        # choose cut points
        cuts_prelim <- unique(c(origin, event_times[partialsum==0], max(event_times)))
        
        # group
        groups_prelim <- cut(sort_tte, cuts_prelim, include.lowest=TRUE)
        
        counts_prelim <- as.integer(table(groups_prelim, useNA="ifany"))
        
        # the last group might not be large enough, if not merge with penultimate group
        if(counts_prelim[length(counts_prelim)] < min_count){
          cuts <- cuts_prelim[-(length(cuts_prelim)-1)]
          groups <- cut(sort_tte, cuts, include.lowest=TRUE)
        } else {
          cuts <- cuts_prelim
          groups <- groups_prelim
        }
        print(cuts)
        counts <- as.integer(table(groups, useNA="ifany"))
      },
      # round events to the nearest `min_count` on the cumulative event count scale
      # then convert back to the events_per_time scale, using `diff()`
      "cumulative" = {
        event_cml <- cumsum(event_n)
        
        ceiling_event_cml <- ceiling_any(event_cml, min_count)
        
        counts <- diff(c(0,ceiling_event_cml))
        cuts <- c(0,event_times)[counts>0]
        counts <- counts[counts>0]
        print(cuts)
      }
    )
  }
  
  switch(value.method,
         "ceiling" = {
           obstime <- rep(cuts[-1], times=counts)
         },
         "floor" = {
           obstime <- rep(cuts[-length(cuts)]+min_increment, times=counts)
         },
         "mid" = {
           obstime <- rep((cuts[-length(cuts)] + cuts[-1])/2, times=counts)
         },
         "spaced" = {
           # duration of interval between cut points
           duration <- diff(cuts)
           
           # starting value within each group
           a <- rep(cuts[-(length(cuts))], times=counts)
           
           # within-group index
           b <- unlist(lapply(counts, seq_len))
           
           # within-group incremental duration
           c <- rep(duration/counts, times=counts)
           
           obstime <- a + b*c
         },
         
  )
  
  # round event times to be integer
  # could leave this as a post-function option, but temptation is do use `as.integer` instead of `ceiling`.
  # Ceiling more appropriate because `as.integer` (or `floor`) will put events times at the origin, which is not allowed.
  obstime <- ceiling(obstime/min_increment)*min_increment
  
  # recover original order from input vector
  obstime_original_order <- obstime[order(order_tte)]
  
  return(obstime_original_order)
}



round_cmlcount <- function(x, time, min_count, method="linear", integer.times=TRUE, integer.counts=TRUE) {
  
  # Purpose:
  # For a vector of cumulative event counts `x`, indexed at `time`, this function
  # rounds `x` so that each time `x` increases, there are at least `min_count` events occurring at each step
  # This is different to the "group event times" approach the `round_event_times` function above.
  # In this function, there _may_ be grouped times where the underlying count is less than `min_count`, but
  # these are rounded in a way that is undetectable.
  
  stopifnot("x must be non-descreasing" = all(diff(x)>=0))
  if(integer.counts) stopifnot("x must be integer" = all(x %% 1 ==0))
  
  # round events such that they are no fewer than min_count events per step
  # steps are then shifted by ` - floor(min_count/2)`, to give  to remove bias
  if(method=="constant") {
    rounded_counts <- roundmid_any(x, min_count)
  }
  
  # as above, but will linearly-interpolate event times between rounded steps,
  # so that events are spaced equally within periods between each `min_count`th event
  # this will also linearly interpolate event if _true_ counts are safe but "steppy" -- can we avoid this by not over-interpolating?
  if(method=="linear") {
    
    # group follow-up time such that the `min_count`th event occurs on the last day within each group
    # eg for:
    # min_count = 3;
    # event_time              = c(0,0,1,0,0,1,0,1,0,1,1,0,0,1,0,1,0,0,0,1,0, 2, 0, 1)... ;
    # cumulative_event_count  = c(0,0,1,1,1,2,2,3,3,4,5,5,5,6,6,7,7,7,7,8,8,10,10,11)... .
    # then:
    # grouped_time            = c(0,0,0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,2,2, 2, 3, 3)...;
    
    x_floor <- floor_any(x, min_count)
    x_ceiling <- lag(x_floor, 1, 0L) + min_count
    
    # an alternative is the `min_count`th event occurs on the first day within each group,so that:
    # grouped_time            = c(0,0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,2,2,2, 3, 3, 3)... .
    # using ceiling in it's pure form.
    # this will start incrementing incidence at the first event time, rather than from time-zero onwards
    
    # x_ceiling <- ceiling_any(x, min_count)
    
    #naturally_steppy <- which((x - x_mid) == 0)
    x_rle <- rle(x_ceiling)
    
    # get index locations of step increases
    steptime <- c(0,time[cumsum(x_rle$lengths)])
    
    # get cumulative count at each step
    stepheight <- c(0,x_rle$values)
    
    rounded_counts <- approx(x=steptime, y=stepheight, xout = time, method="linear", rule=1)$y
    if(integer.times) rounded_counts <- floor(rounded_counts)
  }
  
  return (rounded_counts)
}