funcShaded <- function(x, m, sd, lower_bound, upper_bound) {
  y = dnorm(x, mean = m, sd = sd)
  y[x < lower_bound] <- NA
  y[x > upper_bound] <- NA
  return(y)
}

ggplot(data = data.frame(x = c(-20, 20)), aes(x)) +
  stat_function(fun = dnorm, n = 201, args = list(mean = -5, sd = 5^0.5)) + ylab("") +
  stat_function(fun = funcShaded, args = list(m = -5, sd = 5^0.5, lower_bound = qnorm(0.025,-5,5^0.5), upper_bound = qnorm(0.975,-5,5^0.5)), 
                geom = "area", fill = "#882255", alpha = .2) + 
  stat_function(fun = dnorm, n = 201, args = list(mean = 5, sd = 3^0.5)) + ylab("") +
  stat_function(fun = funcShaded, args = list(m =  5, sd = 3^0.5,
                                              lower_bound = qnorm(0.025, 5,3^0.5),
                                              upper_bound = qnorm(0.975, 5,3^0.5)), 
                geom = "area", fill = "#44aa99", alpha = .2) + 
  scale_y_continuous(breaks = NULL) + 
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank())