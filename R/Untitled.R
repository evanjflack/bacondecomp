library(optparse)
library(yaml)


arg_config <- list( make_option("--config", type = "character",
                                default = "year2_mort"))
arg_parser <- OptionParser(option_list = arg_config)
arg_list   <- parse_args(arg_parser)

message("Loading data...")
config <- read_yaml("~/Desktop/config.yaml")

tuning_design_raw <- do.call(crossing, config$gbm)

library(dplyr)
library()

crossing()

arg_list$config

file.path('', "Desktop")

config
