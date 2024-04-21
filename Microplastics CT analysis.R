##### use a neural net trained on large clusters (i.e. not the edges of materials) 
##### of pixels classified as of a certain density and see how that improves
##### classification

# Megan M. Trusler, Craig J. Sturrock, Christopher H. Vane, Sarah Cook, Barry H. Lomax,
# X-ray computed tomography: A novel non-invasive approach for the detection of microplastics in sediments?,
# Marine Pollution Bulletin,
# Volume 194, Part A,
# 2023,
# 115350,
# ISSN 0025-326X,
# https://doi.org/10.1016/j.marpolbul.2023.115350.
# (https://www.sciencedirect.com/science/article/pii/S0025326X23007841)
# Abstract: As a non-invasive imaging technique, this study explores the application of Computed Tomography (CT) in microplastics research, assessing its potential to distinguish different types and sizes of microplastics (polypropylene, polyethylene terephthalate, polyethylene, and polyvinyl chloride) from homogenised river-estuarine sediment. When examined in layers within artificial cores, all microplastic types could be observed by CT imagery, with good contrast in X-ray attenuation (based on image gray level intensity) against background sediments. Large microplastics (4 mm diameter) were also detectable when distributed randomly amongst the sediment. These spiked cores had sufficient difference in attenuation to allow segmentation between type, and therefore isolate individual microplastics. Due to limitations on scan resolution, smaller microplastics (≤125 μm diameter) could not be detected in spiked cores. Scans of two sediment cores from a Thames River tributary (UK) revealed two distinctive sediment structures which could influence microplastic accumulation. This information would be lost using conventional recovery procedures.
# Keywords: Particle imaging; Plastic polymer; Sediment core; Microplastic pollution; Image analysis; River Thames

# loading packages
library(dplyr)
library(pixmap)
library(parallel)
library(purrr)
library(reshape2)
library(scico)

### the scan files are quite large (2764048 pixels)
(original_pixel_count <- read.pnm("/Users/mukappa/Documents/Megan Trusler microplastics/Microplastics core 3.ppm")@size %>% prod())

### to make this script run faster, we'll reduce the resolution using these functions
# this just rotates the image matrix so it's in the correct orientation
rotate <- function(x) t(apply(x, 2, rev))

# reduce the resolution by averaging (1 + spread)^2 pixels, if the amounts of rows and columns are not divisible by spread + 1 (capture in the function) then they are first trimmed from the image
mat_cell_avg <- function(mat, spread = 1, parallelise = T){
  capture <- 1 + spread
  diffrows <- nrow(mat) %% capture
  if(diffrows != 0){
    mat <- mat[-(rev(1:nrow(mat))[1:diffrows]), ]
    warning(paste0("There were ", 
                   diffrows, 
                   " excess rows which were trimmed from the end of the data"))
  }
  diffcols <- ncol(mat) %% capture
  if(diffcols != 0){
    mat <- mat[, -(rev(1:ncol(mat))[1:diffcols])]
    warning(paste0("There were ", 
                   diffcols, 
                   " excess columns which were trimmed from the end of the data"))
  }
  matseqc <- seq(capture, ncol(mat), by = capture) - spread
  matseqr <- seq(capture, nrow(mat), by = capture) - spread
  p0 <- expand.grid(x0 = matseqc, y0 = matseqr)
  p1 <- p0 + spread
  names(p1) <- c("x1", "y1")
  matseqtab <- cbind(p0, p1)
  matseqtab$xnew <- matseqtab$x1/capture
  matseqtab$ynew <- matseqtab$y1/capture
  matseqtab$value <- NA
  if(parallelise){
    matseqlist <- apply(matseqtab[, 1:4], 1, as.list)
    matseqlist <- lapply(matseqlist, as.numeric)
    matseqtab$value <- as.numeric(mclapply(matseqlist, function(x){mean(mat[x[2]:x[4], x[1]:x[3]])}))
  } 
  else{
    for(i in 1:nrow(matseqtab)){
      coords <- as.numeric(matseqtab[i, 1:4])
      matseqtab$value[i] <- mean(mat[coords[2]:coords[4], coords[1]:coords[3]])
    }
  }
  matrix(matseqtab$value, ncol = max(matseqtab$xnew), byrow = T)
}

# read.pnm reads in 3 colour channels, but as the ct image is in grayscale we can choose any channel (they are the same)
# the function then rotates the input file and averages the cells according to the spread
# the resulting image will be approx (spread + 1)^2 fewer pixels
reduce_resolution <- function(ppm, spread = 1){
  mp <- read.pnm(ppm)
  mp <- mp@red
  mp <- rotate(mp)
  mat_cell_avg(mp, spread = spread)
}

mp <- reduce_resolution("/Users/mukappa/Documents/Megan Trusler microplastics/Microplastics core 3.ppm", spread = 3)

# number of pixels in the matrix
(new_res_pixels <- length(mp))
original_pixel_count /  new_res_pixels

# for plotting, we'll create another image by melting the image into x, y, and value
mpmelt <- reshape2::melt(mp)

# and using a simple function to plot the image
# it can convert the image into a negative as some ct images can be dark
# it converts the decimal greyscale value into a hexidecimal one
plot_ct <- function(mpmelt, blank = F, image_negative = F){
  if(image_negative){mpmelt$value <- 1 - mpmelt$value}
  mpmelt$value <- paste0("#000000", toupper(as.hexmode(round(mpmelt$value * 255))))
  plot(mpmelt$Var1, mpmelt$Var2, type = "n", xlab = "x", ylab = "y", las = 1)
  if(!blank){
    rect(mpmelt$Var1, mpmelt$Var2, mpmelt$Var1 + 1, mpmelt$Var2 + 1, col = mpmelt$value, border = NA)
  }
}

par(mfrow = c(1, 2))
# here's the ct scan of the core containing the microplastics that appear as light grey blobs
plot_ct(mpmelt)
# it looks like we can pick them out based on colour
xleft = 50; xright = 95; ybottom = 450; ytop = 475
# so let's have a closer look at one of the particles
rect(xleft, ybottom, xright, ytop, lwd = 2, border = "#0000FF")
grey_values <- mpmelt %>% filter(Var1 > xleft & Var1 < xright & Var2 > ybottom & Var2 < ytop)
# there's a range of gray values
grey_values %>% pull(value) %>% hist(breaks = 20, las = 1, xlim = c(0, 1), col = "#0000FF55")
check_grey <- mpmelt[mpmelt$value > min(grey_values$value) & mpmelt$value < max(grey_values$value), ]

# we can't take those grey values as they occur in other parts of the scan which do not pertain to microplastics
par(mfrow = c(1, 2))
plot_ct(mpmelt)
rect(check_grey$Var1, check_grey$Var2, check_grey$Var1 + 1, check_grey$Var2 + 1, col = "#0000FF", border = NA)
# when the background is removed we can see that it's not picking out just the microplastics particles, or that pixels are missing inside some of them
plot_ct(mpmelt, blank = T)
rect(check_grey$Var1, check_grey$Var2, check_grey$Var1 + 1, check_grey$Var2 + 1, col = "#000000", border = NA)

# to our eyes we can see continuous areas of shades of grey and we don't need prior knowledge of what those shades are
# we work from the point of those 

# to improve things we can think of the pixels of the microplastic particles as being related to one another;
# pixels pertaining to microplastics particles should have neighbouring pixels within the range of the greyscale we discovered above
# to get the values of a pixel's nearest neighbours


# show a pca!
read_pixels <- function(mat, n, spread = 2){
  init_row <- c(n, n + (1:spread * nrow(mat)), n - (1:spread * nrow(mat)))
  c(
    init_row,
    sapply(1:spread, function(x){init_row + x}),
    sapply(1:spread, function(x){init_row - x})
  ) %>% sort()
}

ct_pca <- function(mp, spread = 1){
  pixel_ids <- matrix(1:length(mp), nrow = nrow(mp), ncol = ncol(mp))
  xminexc <- unique(mpmelt$Var1)[1:spread]
  xmaxexc <- rev(unique(mpmelt$Var1))[1:spread]
  xexc <- sort(c(xminexc, xmaxexc))
  yminexc <- unique(mpmelt$Var2)[1:spread]
  ymaxexc <- rev(unique(mpmelt$Var2))[1:spread]
  yexc <- sort(c(yminexc, ymaxexc))
  cropped_pixels <- pixel_ids[-xexc, -yexc]
  pixel_seq <- sort(as.numeric(cropped_pixels))
  pixel_groups <- mclapply(pixel_seq, function(n){c(n, mp[read_pixels(mp, n, spread = spread)])})
  pixel_groups <- do.call(rbind, pixel_groups)
  names(pixel_groups) <- c("pixel_ID", paste0("pixel", 1:(ncol(pixel_groups) + 1)))
  pixel_pca <- prcomp(pixel_groups[, 2:ncol(pixel_groups)], center = T)
  list(pixel_ids = pixel_ids, cropped_pixels = cropped_pixels, pixel_pca = pixel_pca)
}

ct_cluster <- function(pixel_pca, n_clusters = 3, n_start = 5){
  kp <- kmeans(pixel_pca$x, centers = n_clusters, nstart = n_start)
  print(table(kp$cluster))
  kp
}

ct_cluster_elbow <- function(ct_pca, k_vals = 1:10){
  tot_withinss <- purrr::map_dbl(k_vals, function(k){
    model <- kmeans(ct_pca$x, centers = k)
    model$tot.withinss
  })
  elbow_df <- data.frame(k = 1:10, tot_withinss = tot_withinss)
  ggplot(elbow_df, aes(k, tot_withinss)) + geom_line() + geom_point() + scale_x_continuous(breaks = 1:10)
  elbow_df
}

ct_cluster_overlay <- function(mpmelt, n_clusters = 3, alpha = 1){
  mpmelt_col <- mpmelt
  mpmelt_col$id <- sapply(1:nrow(mpmelt), function(x){pixel_ids[mpmelt$Var1[x], mpmelt$Var1[x]]})
  mpmelt_col <- mpmelt_col
  mpmelt_col$cluster <- NA
  mpmelt_col$cluster[pixel_ids %in% cropped_pixels] <- kp$cluster
  mpmelt_col$col <- NA
  mpmelt_col$col[pixel_ids %in% cropped_pixels] <- scico::scico(n_clusters + 1, palette = "roma")[1:n_clusters][kp$cluster]
  mpmelt_col$col[which(!is.na(mpmelt_col$col))] <- paste0(mpmelt_col$col[which(!is.na(mpmelt_col$col))], as.hexmode(round(255 * alpha)))
  rect(mpmelt$Var1, mpmelt$Var2, mpmelt$Var1 + 1, mpmelt$Var2 + 1, col = mpmelt_col$col, border = NA)
}

###############################

ct_pca_dat <- ct_pca(mp, spread = 2)
pixel_pca <- ct_pca_dat$pixel_pca
cropped_pixels <- ct_pca_dat$cropped_pixels
pixel_ids <- ct_pca_dat$pixel_ids

ct_cluster_elbow(pixel_pca)
kp <- ct_cluster(pixel_pca)

par(mfrow = c(1, 2))
plot_ct(mpmelt)
plot_ct(mpmelt)
ct_cluster_overlay(mpmelt)

mp <- reduce_resolution("/Users/mukappa/Documents/Megan Trusler microplastics/Microplastics core 1.ppm", spread = 3)


### messing about with colour reduction
n <- 6
reduced_palette <-seq(0, 1, length.out = 16)
new_colours <- sapply(mpmelt$value, function(x){reduced_palette[which.min(abs(x - reduced_palette))]})
new_hex_values <- as.hexmode(round(new_colours * 255))
col_df <- data.frame(table(round(mp$image_matrix, 3)))
names(col_df) <- c("value", "count")
col_df$value <- as.numeric(col_df$value)
par(mfrow = c(1, 1))
plot(count ~ value, col_df, type = "l")
diff_vals <- diff(col_df$count, differences = 2)
abundant_cols <- col_df[which(diff_vals < -100), ]
abundant_cols
plot(diff_vals, type = "l")
abline(h = -125, col = "red")
abline(v = 555, col = "blue")
abline(v = 700, col = "blue")
col_peaks <- abundant_cols$value[abundant_cols$value > 575 & abundant_cols$value < 700]
(max(col_peaks) - min(col_peaks)) / length(col_peaks)

mpmelt <- mp$image_melted
mpmelt$value <- round(mpmelt$value, 1)
table(mpmelt$value)
table(mpmelt$value) %>% length()
mp <- mp$image_matrix
mp <- round(mp, 1)

ct_pca_dat <- ct_pca(mp, spread = 2)
pixel_pca <- ct_pca_dat$pixel_pca
cropped_pixels <- ct_pca_dat$cropped_pixels
pixel_ids <- ct_pca_dat$pixel_ids

ct_elbow <- ct_cluster_elbow(pixel_pca)
ggplot(ct_elbow, aes(k, tot_withinss)) + geom_line() + geom_point() + scale_x_continuous(breaks = 1:10)

n_clusters <- 5

kp <- ct_cluster(pixel_pca, n_clusters = n_clusters)

par(mfrow = c(1, 2))
plot_ct(mpmelt)
plot_ct(mpmelt, blank = T)
ct_cluster_overlay(mpmelt, n_clusters = n_clusters, alpha = 1)

mp <- reduce_resolution("/Users/mukappa/Documents/Megan Trusler microplastics/Microplastics core 4.ppm", spread = 3)



