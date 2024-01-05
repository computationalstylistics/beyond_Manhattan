# load the corpora for 100 EN novels, and 100 PL novels

load("freqs_100_novels_EN.RData")
load("freqs_100_novels_PL.RData")

# load the corpora with Southern authors, and fantasy authors

library(stylo)
data(lee)
data(galbraith)




#############################################################
#
# define functions
#
#############################################################


# function to compute a distance measure with custom norm l
# (please keep in mind that the function is not optimized for speed!)
dist_l_norm = function(frequencies, l = 1, scale = TRUE) {

	no_of_txts = dim(frequencies)[1]
	no_of_dims = dim(frequencies)[2]
	dist_matrix = matrix(nrow = no_of_txts, ncol = no_of_txts, 
		dimnames = list(rownames(frequencies), rownames(frequencies)))
	if(scale == TRUE) {
		frequencies = scale(frequencies)
	}
	for(source_txt in 1:(no_of_txts -1)) {
		for(target_txt in (source_txt +1):no_of_txts) {
			raw_diff = frequencies[source_txt, ] - frequencies[target_txt, ]
			### get rid of 0 differences
			### this will not harm for powers >0, while sanitizing powers <=0
			raw_diff = raw_diff[raw_diff != 0]
			delta = sum(abs(raw_diff) ^l) ^(1/l)
			dist_matrix[target_txt, source_txt] = delta
		}
	}
	# optionally divide by the number of dimensions
	dist_matrix = dist_matrix / no_of_dims
	return(as.dist(dist_matrix))
}




#############################################################




# function to perform leave-one-out Delta classification
leave_one_out_classification = function(freq, l_norm = 1) {

	distance_table = dist_l_norm(freq, l = l_norm, scale = TRUE)
	distance_table = as.matrix(distance_table)

	predicted_classes = c()
	expected_classes = c()

	# iterate over the rows of the table
	for(i in 1:length(distance_table[,1]) ) {
	    # for each row, sort the distances, get the second closest
	    # (the second, because the first is 0, for itself)
	    nearest_meighbors = sort(distance_table[i,])
	    # get its name, and from the name get the class ID
	    predicted_class = gsub("_.*", "", names(nearest_meighbors[2]))
	    expected_class = gsub("_.*", "", names(nearest_meighbors[1]))
	    predicted_classes = c(predicted_classes, predicted_class)
	    expected_classes = c(expected_classes, expected_class)
	}

	results = performance.measures(predicted_classes, expected_classes)$avg.f
	return(results)

}



#############################################################



# perform the above classification x 100 with randomly selected feature subspace
# (the default value is 0.5, to randomly pick 50% of the input features)
bootstrapped_leave_one_out = function(input_features, features = 0.5, l_norm = 1) {

	# how many variables (features) we have:
	feature_space = dim(input_features)[2]
	f_scores = c()

	for(n in 1:100) {
		# get randomly a subset, namely a proportion provided by 'features'
		features_subset = sample(feature_space, floor(feature_space * features))
		# select a table of frequencies for the current feature subset
		table_of_frequencies = input_features[, features_subset]
		f = leave_one_out_classification(table_of_frequencies, l_norm = l_norm)
		f_scores = c(f_scores, f)
	}
	f_scores = mean(f_scores)
	#
	return(f_scores)
}



#############################################################



# function to perform bootstrapped Delta (see above) iteratively,
# in order to cover different values of the norm l
assess_all_l_norms = function(frequencies, grid_l = seq(0.1, 10, 0.1), ...){

	f_scores_all = c()
	for(l in grid_l) {
		current_f = bootstrapped_leave_one_out(frequencies, l_norm = l, ...)
		f_scores_all = c(f_scores_all, current_f)
	}

	names(f_scores_all) = grid_l
	return(f_scores_all)
}



#############################################################



# function to run the above procedure iteratively over different MFW ranges
iterate_over_mfw_ranges = function(frequencies, mfw = seq(100, 1000, 100), ...){

	results_all = c()
	for(current_mfw in mfw) {
		message(current_mfw)
		results = assess_all_l_norms(frequencies[, 1:current_mfw], ...)
		results_all = rbind(results_all, results)
	}
	rownames(results_all) = mfw
	return(results_all)

}







#############################################################
#
# run the procedure
#
#############################################################

results_EN = iterate_over_mfw_ranges(freqs_100_novels_EN, mfw = seq(100, 1000, 100), grid_l = seq(-10, 10, 0.1))
save(results_EN, file = "freqs_100_novels_EN_RESULTS.RData")

results_PL = iterate_over_mfw_ranges(freqs_100_novels_PL, mfw = seq(100, 1000, 100), grid_l = seq(-10, 10, 0.1))
save(results_PL, file = "freqs_100_novels_PL_RESULTS.RData")

results_LEE = iterate_over_mfw_ranges(lee, mfw = seq(100, 1000, 100), grid_l = seq(-10, 10, 0.1))
save(results_LEE, file = "freqs_100_novels_LEE_RESULTS.RData")

results_GALBRAITH = iterate_over_mfw_ranges(galbraith, mfw = seq(100, 1000, 100), grid_l = seq(-10, 10, 0.1))
save(results_GALBRAITH, file = "freqs_100_novels_GALBRAITH_RESULTS.RData")






#############################################################
#
# visualize
#
#############################################################



library(viridis)
load("freqs_100_novels_PL_RESULTS.RData")
load("freqs_100_novels_EN_RESULTS.RData")
load("freqs_100_novels_GALBRAITH_RESULTS.RData")
load("freqs_100_novels_LEE_RESULTS.RData")



#png(file = "results_GALBRAITH.png", width = 7, height = 5, units = "in", res = 300)

current_results = results_PL
no_of_results = dim(current_results)[1]
f1_min = min(current_results)
f1_max = max(current_results)
axis_x = as.numeric(colnames(current_results))
color_palette = viridis(no_of_results)

#f1_min = 0.5
#f1_max = 0.8

plot(1:length(axis_x) ~ axis_x, type = "n", ylim = c(f1_min, f1_max), xlim = c(-2,7),
	ylab = "performance (f1 value)", xlab = "L norm (power)")
segments(1, 0.77, 1, 1, lty = 3, lwd = 2, col = "darkgray")
segments(2, 0.77, 2, 1, lty = 3, lwd = 2, col = "darkgray")
text(1, 0.77, "L1 norm  ", srt = 90, adj = 1, col = "darkgray")
text(2, 0.77, "L2 norm  ", srt = 90, adj = 1, col = "darkgray")
for(i in 1:no_of_results) {
	#lines(lowess(current_results[i,])$y ~ axis_x)
	lines(current_results[i,] ~ axis_x, col = color_palette[i], lwd = 3)
}
legend("bottomleft",
       legend = c("100", "500", "1000"), 
       col = c(color_palette[c(1, 5, 10)]),
       lwd = 3,
       title = "number of MFWs",
       bty = "n" # no box around the legend
       )





# fitting a linear model on max values across different MFW ranges
#
#max_f1_values = apply(current_results, 1, max)
#find_max_location = apply(current_results, 2, function(x) x == max_f1_values)
#find_max_loc = apply(find_max_location, 1, which)
#find_max = as.numeric(colnames(current_results)[find_max_loc])
#abline(lm(max_f1_values ~ find_max))





