# Function to determine elbow
# Source: https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve
################################################################################
elbow_finder <- function(x_values, y_values) {
        # Max values to create line
        max_x_x <- max(x_values)
        max_x_y <- y_values[which.max(x_values)]
        max_y_y <- max(y_values)
        max_y_x <- x_values[which.max(y_values)]
        max_df <- data.frame(x = c(max_y_x, max_x_x), y = c(max_y_y, max_x_y))
        
        # Creating straight line between the max values
        fit <- lm(max_df$y ~ max_df$x)
        
        # Distance from point to line
        distances <- c()
        for(i in 1:length(x_values)) {
                distances <- c(distances, abs(coef(fit)[2]*x_values[i] - y_values[i] + coef(fit)[1]) / sqrt(coef(fit)[2]^2 + 1^2))
        }
        
        # Max distance point
        x_max_dist <- x_values[which.max(distances)]
        y_max_dist <- y_values[which.max(distances)]
        
        c(x_max_dist, y_max_dist)
}
################################################################################
