loadlib = function(x) {
    if (!require(x, character.only = TRUE)) {
        install.packages(x)
    } else {
        message(paste0('Loaded required library ', x))
    }
}
