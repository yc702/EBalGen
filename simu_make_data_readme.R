## Help functions for examples in the readme demo.
make_data <- function(n,
                      p = 5,
                      setting_x = 1,
                      setting_src = 1,
                      setting_trt = 1, 
                      setting_main = 1,
                      setting_cate = 1,
                      err_sd = 1)
{

   x_fn <- function(n) {
     switch(setting_x,
            `1` = matrix(4 * runif(n * p) - 2, n, p),
            `2` = matrix(8 * runif(n * p) - 2, n, p))
   }
   
   prob_src_fn <- function(x) {
      switch(setting_src,
             `1` = plogis(.4 * x[,1] + .3 * x[,2] - .2 * x[,4]))
   }
   
   prob_trt_fn <- function(x) {
      z <- switch(setting_trt,
                  `1` = (.7 * x[,2] + .5 * x[,3]))
      plogis(z)
   }
   
   eff_cate_fn <- function(x) {
      
      z <- switch(setting_cate,
                  `1` = (x[,1] - .6 * x[,2] - .4 * x[,3]))
      z

   }
   
   eff_main_fn <- function(x) {
      switch(setting_main,
             `1` = .5 * x[,1] + .3 * x[,2] + .3 * x[,3]- .4 * x[,4] - .7 * x[,5])
   }
   
   err_fn <- function(n) rnorm(n, 0, err_sd)
   
   .make_data_inner(n, x_fn, err_fn, prob_src_fn, prob_trt_fn, eff_main_fn, eff_cate_fn)
   
}

.make_data_inner <- function(n, x_fn, err_fn, prob_src_fn, prob_trt_fn, eff_main_fn, eff_cate_fn)
{
   x <- x_fn(n)
   prob_src <- prob_src_fn(x)
   prob_trt <- prob_trt_fn(x)
   eff_main <- eff_main_fn(x)
   eff_cate <- eff_cate_fn(x)
   err <- err_fn(n)
   
   s <- runif(n) <= prob_src
   trt <- runif(n) <= prob_trt
   
   y0 <- eff_main - eff_cate / 2
   y1 <- eff_main + eff_cate / 2
   
   y <-  (trt * y1 + (1 - trt) * y0 + err)
   
   list(x   = x,
        s   = as.integer(s),
        trt = trt,
        y   = y,
        prob_src = prob_src,
        prob_trt = prob_trt,
        y0  = y0,
        y1  = y1)
}



