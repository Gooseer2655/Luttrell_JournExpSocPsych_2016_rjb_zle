display.lm <- function(model, rnd = 2, rnd.p.two = FALSE) {
  
  if(rnd <= 2 & rnd.p.two == FALSE) {
    cutoff <- .001
    p.display <- str_c("< ", cutoff)
    
    display.obj <- tidy(model) %>% bind_cols(confint_tidy(model)) %>% 
      as.data.frame() %>% mutate(p.value.char = ifelse(p.value < cutoff, p.display, p.value),
                                 ` ` = ifelse(p.value < .001, "***", ifelse(p.value < .01, "**", 
                                                                            ifelse(p.value < .05, "*", ifelse(p.value < .10, "+", " "))))) %>%
      mutate_if(is.numeric, round, 3) %>%
      mutate(p.value = ifelse(p.value < cutoff, p.display, p.value)) %>%
      select(-p.value.char) %>%
      mutate_if(is.numeric, round, 2)
  } else {
    cutoff <- str_c("0.", paste(rep(0,rnd-1), collapse = ""), "1") %>% as.numeric()
    p.display <- str_c("< ", cutoff)
    
    display.obj <- tidy(model) %>% bind_cols(confint_tidy(model)) %>% 
      as.data.frame() %>% mutate(p.value.char = ifelse(p.value < cutoff, p.display, p.value),
                                 ` ` = ifelse(p.value < .001, "***", ifelse(p.value < .01, "**", 
                                                                            ifelse(p.value < .05, "*", ifelse(p.value < .10, "+", " "))))) %>%
      mutate_if(is.numeric, round, rnd) %>%
      mutate(p.value = ifelse(p.value < cutoff, p.display, p.value)) %>%
      select(-p.value.char)
  }
  return(display.obj)
}
