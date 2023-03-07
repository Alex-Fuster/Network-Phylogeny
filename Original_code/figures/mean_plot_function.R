MeanPlotPosNeg <- function(measure, thing_we_measure){
	mean_measure_pos <- apply(measure[,1:nsim/2],1, mean, na.rm = TRUE)
	quantiles_pos <- apply(measure[,1:nsim/2], 1, quantile, probs = c(0.05, 0.95), na.rm = T)

	y_min_pos <- min(quantiles_pos, na.rm = T)
	y_max_pos <- max(quantiles_pos, na.rm = T)

	if(y_min_pos < 0 | y_max_pos > 1){
		y_range_pos = c(y_min_pos, y_max_pos)
	}else{
		y_range_pos = c(0, 1)
	}

	mean_measure_neg <- apply(measure[,((nsim/2)+1):nsim],1, mean, na.rm = TRUE)
	quantiles_neg <- apply(measure[,((nsim/2)+1):nsim], 1, quantile, probs = c(0.05, 0.95), na.rm = T)

	y_min_neg <- min(quantiles_neg, na.rm = T)
	y_max_neg <- max(quantiles_neg, na.rm = T)

	if(y_min_neg < 0 | y_max_neg > 1){
		y_range_neg = c(y_min_neg, y_max_neg)
	}else{
		y_range_neg = c(0, 1)
	}

	## Plot for Positive interactions
	par(mar = c(3,3,0.5,0.5), mgp = c(1.5, 0.3, 0), tck = -.008, las = 1)
	plot(mean_measure_pos[1:150], pch=19, xlab = "Time", xlim = c(0,150), ylab = paste("Mean", thing_we_measure), ylim = y_range_pos, col = c(rgb(0, 158, 115, max = 255)))
	polygon(c(1:150,rev(1:150)), c(quantiles_pos[1,1:150], rev(quantiles_pos[2,1:150])), col = adjustcolor(rgb(161,215,106, max = 255), alpha = 0.5), border = NA)
	points(mean_measure_neg[1:150], pch=19, col = c(rgb(230, 159, 0, max = 255)))
	polygon(c(1:150,rev(1:150)), c(quantiles_neg[1,1:150], rev(quantiles_neg[2,1:150])), col = adjustcolor(rgb(230, 159, 0, max = 255), alpha = 0.3), border = NA)

}
