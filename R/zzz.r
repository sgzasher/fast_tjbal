utils::globalVariables(c("id", "xmin", "xmax", "ymin", "ymax"))

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("fast.tjbal: Trajectory Balancing via KMMD")
}
