# Function to append a list into an existing list
list_append <- function(lst, ...){
  lst <- c(lst, list(...))
  return(lst)
}
