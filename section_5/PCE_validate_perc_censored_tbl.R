censoring_proportion_table <- function(seeds, data_file = "PCE.RData") {
  load(data_file)  # loads data.PCE.censored
  dat_us <- data.PCE.censored
  
  states <- spData::us_states
  pts <- sf::st_as_sf(
    dat_us, coords = c("lon", "lat"),
    crs = 4326, remove = FALSE
  )
  texas <- states[states$NAME == "Texas", ]
  in_tx <- lengths(sf::st_intersects(
    pts, sf::st_transform(texas, 4326)
  )) > 0
  dat_tx <- dat_us[in_tx, ]
  
  heldout_prop <- function(dat, seed) {
    group <- paste(dat$lon, dat$lat, sep = ":")
    set.seed(seed)
    held <- sample(
      unique(group),
      max(1, round(0.2 * length(unique(group))))
    )
    mean(dat$left_censored[group %in% held])
  }
  
  proportions <- rbind(
    TX = vapply(seeds, heldout_prop, numeric(1), dat = dat_tx),
    US = vapply(seeds, heldout_prop, numeric(1), dat = dat_us)
  )
  colnames(proportions) <- paste("Seed", seeds)
  
  table <- data.frame(
    Region = rownames(proportions),
    proportions,
    row.names = NULL,
    check.names = FALSE
  )
  
  knitr::kable(
    table,
    format = "latex",
    digits = 3,
    row.names = FALSE,
    align = c("l", rep("c", length(seeds)))
  )
}

latex_table <- censoring_proportion_table(1:5)
cat(latex_table, sep = "\n")