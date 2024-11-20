grid_sort_d <- grid %>% 
  mutate(prob = ifelse(dhalflife > 0, 1-prob, prob),
         abs_dh = abs(dhalflife)) %>% 
  filter(combo == "7")
grid_sort_d = grid_sort_d[order(grid_sort_d$abs_dh),]

x = grid_sort_d$abs_dh[2:nrow(grid_sort_d)]
y = grid_sort_d$prob
sfun_d  <- stepfun(x, y, f=0)

tt_d1 <- seq(min(grid_sort_d$abs_dh), max(grid_sort_d$abs_dh), by = 0.01)
tt_d2 <- seq(min(grid_sort_d$abs_dh), max(grid_sort_d$abs_dh), by = 0.05)
tt_d3 <- seq(min(grid_sort_d$abs_dh), max(grid_sort_d$abs_dh), by = 0.1)


grid_sort_r <- grid %>% 
  mutate(prob = ifelse(rhalflife > 1, 1-prob, prob),
         h_times_gt = ifelse(rhalflife < 1, 1/rhalflife-1, rhalflife-1)) %>% 
  filter(combo == "7")
grid_sort_r = grid_sort_r[order(grid_sort_r$h_times_gt),]

x = grid_sort_r$h_times_gt[2:nrow(grid_sort_r)]
y = grid_sort_r$prob
sfun_r  <- stepfun(x, y, f=0)

tt_r1 <- seq(min(grid_sort_r$h_times_gt), max(grid_sort_r$h_times_gt), by = 0.1)
tt_r2 <- seq(min(grid_sort_r$h_times_gt), max(grid_sort_r$h_times_gt), by = 0.2)
tt_r3 <- seq(min(grid_sort_r$h_times_gt), max(grid_sort_r$h_times_gt), by = 0.5)


op <- par(mfrow = c(3,2))
plot(sfun_r, xval = tt_r1, col.hor = "red", main="small steps", xlab="", ylab="");
plot(sfun_d, xval = tt_d1, col.hor = "red", main="", xlab="", ylab="");
plot(sfun_r, xval = tt_r2, col.hor = "orange", main="medium steps", xlab="", ylab="P(halflife_i > halflife_j); h_i > h_j");
plot(sfun_d, xval = tt_d2, col.hor = "orange", main="", xlab="", ylab="");
plot(sfun_r, xval = tt_r3, col.hor = "blue", main="large steps", xlab="n times h_i greater than h_j", ylab="");
plot(sfun_d, xval = tt_d3, col.hor = "blue", main="", xlab="h_i - h_j", ylab="")






