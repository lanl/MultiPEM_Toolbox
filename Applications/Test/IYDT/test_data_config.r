data_dir <- file.path(tdir, "data")
schemas <- list(
  seismic_cal.csv = list(
    columns = c("Y1", "Y2", "Source", "Path", "Station", "Type",
                "lRange", "W", "C2N", "HOB"),
    model = c("Y1", "Y2", "Type", "lRange", "W", "C2N", "HOB")
  ),
  seismic_new.csv = list(
    columns = c("Y1", "Y2", "Source", "Path", "Station", "Type",
                "lRange", "W", "HOB"),
    model = c("Y1", "Y2", "Type", "lRange", "W", "HOB")
  ),
  acoustic_cal.csv = list(
    columns = c("Y1", "Y2", "Source", "Path", "Station", "Type",
                "logTempSc", "logPressureSc", "lRange", "W", "C2N", "HOB"),
    model = c("Y1", "Y2", "Type", "logTempSc", "logPressureSc",
              "lRange", "W", "C2N", "HOB")
  ),
  acoustic_new.csv = list(
    columns = c("Y1", "Y2", "Source", "Path", "Station", "Type",
                "logTempSc", "logPressureSc", "lRange", "W", "HOB"),
    model = c("Y1", "Y2", "Type", "logTempSc", "logPressureSc",
              "lRange", "W", "HOB")
  ),
  optical_cal.csv = list(
    columns = c("Y1", "Y2", "Source", "logTempSc", "logPressureSc", "W", "HOB"),
    model = c("Y1", "Y2", "logTempSc", "logPressureSc", "W", "HOB")
  ),
  optical_new.csv = list(
    columns = c("Y1", "Y2", "Source", "logTempSc", "logPressureSc", "W", "HOB"),
    model = c("Y1", "Y2", "logTempSc", "logPressureSc", "W", "HOB")
  ),
  crater_cal.csv = list(
    columns = c("Y1", "Y2", "Source", "W", "HOB"),
    model = c("Y1", "Y2", "W", "HOB")
  ),
  crater_new.csv = list(
    columns = c("Y1", "Y2", "Source", "W", "HOB"),
    model = c("Y1", "Y2", "W", "HOB")
  )
)

read_fixture <- function(name, schema)
{
  value <- read.csv(
    file.path(data_dir, name), check.names = FALSE, stringsAsFactors = FALSE
  )
  if (!identical(names(value), schema$columns)) {
    stop(sprintf("%s has an unexpected column contract", name), call. = FALSE)
  }
  if (!all(vapply(value[schema$model], is.numeric, logical(1)))) {
    stop(sprintf("%s has a nonnumeric model column", name), call. = FALSE)
  }
  if (!all(is.finite(as.matrix(value[schema$model])))) {
    stop(sprintf("%s has a non-finite model value", name), call. = FALSE)
  }
  if (!nrow(value)) stop(sprintf("%s is empty", name), call. = FALSE)
  value
}

datasets <- Map(read_fixture, names(schemas), schemas)
names(datasets) <- names(schemas)

seismic_beta <- c(1.1, -0.65, 0.4, 0.08, -0.3)
seismic_x <- as.matrix(datasets[["seismic_cal.csv"]][
  , c("lRange", "W", "C2N", "HOB"), drop = FALSE
])
seismic_params <- list(
  pbeta = 5, iresp = TRUE, yield_scaling = 1 / 3, X = seismic_x,
  cal = FALSE, cal_par_names = character(), ncalp = 0,
  theta_names = character(), notExp = notExp, dnotExp = dnotExp
)
check_jacobian(
  "IYDT seismic calibration fixture", f_s, g_s, seismic_beta, seismic_params
)

seismic_new <- datasets[["seismic_new.csv"]]
seismic_event <- list(
  beta = seismic_beta, theta_names = c("W", "HOB"), iresp = TRUE,
  yield_scaling = 1 / 3,
  X = as.matrix(seismic_new[, "lRange", drop = FALSE]), notExp = notExp
)
check_jacobian(
  "IYDT seismic event fixture", f0_s, g0_s,
  unlist(seismic_new[1, c("W", "HOB")], use.names = FALSE), seismic_event
)

acoustic_beta <- c(0.7, -0.45, -0.8)
acoustic_x <- as.matrix(datasets[["acoustic_cal.csv"]][
  , c("logTempSc", "logPressureSc", "lRange", "W", "C2N", "HOB"),
  drop = FALSE
])
acoustic_params <- list(
  pbeta = 3, iresp = TRUE, yield_scaling = 1 / 3,
  pressure_scaling = 1 / 3, temp_scaling = 1 / 2, X = acoustic_x,
  cal = FALSE, cal_par_names = character(), ncalp = 0,
  theta_names = character()
)
check_jacobian(
  "IYDT acoustic calibration fixture", f_a, g_a, acoustic_beta, acoustic_params
)

acoustic_new <- datasets[["acoustic_new.csv"]]
acoustic_event <- list(
  beta = acoustic_beta, theta_names = c("W", "HOB"), iresp = TRUE,
  yield_scaling = 1 / 3, pressure_scaling = 1 / 3, temp_scaling = 1 / 2,
  X = as.matrix(acoustic_new[
    , c("logTempSc", "logPressureSc", "lRange"), drop = FALSE
  ])
)
check_jacobian(
  "IYDT acoustic event fixture", f0_a, g0_a,
  unlist(acoustic_new[1, c("W", "HOB")], use.names = FALSE), acoustic_event
)

optical_beta <- c(-1, 0.25, 0.4, -0.2)
optical_x <- as.matrix(datasets[["optical_cal.csv"]][
  , c("W", "HOB"), drop = FALSE
])
optical_params <- list(
  pbeta = 4, yield_scaling = 1 / 3, X = optical_x, cal = FALSE,
  cal_par_names = character(), ncalp = 0, theta_names = character(),
  notExp = notExp, dnotExp = dnotExp
)
check_jacobian(
  "IYDT optical calibration fixture", f_o, g_o, optical_beta, optical_params
)

optical_new <- datasets[["optical_new.csv"]]
optical_event <- list(
  beta = optical_beta, theta_names = c("W", "HOB"),
  yield_scaling = 1 / 3,
  X = as.matrix(optical_new[, c("logTempSc", "logPressureSc"), drop = FALSE]),
  notExp = notExp
)
check_jacobian(
  "IYDT optical event fixture", f0_o, g0_o,
  unlist(optical_new[1, c("W", "HOB")], use.names = FALSE), optical_event
)

crater_beta <- c(0.5, 1 / 3)
crater_x <- as.matrix(datasets[["crater_cal.csv"]][
  , c("W", "HOB"), drop = FALSE
])
crater_params <- list(
  pbeta = 2, X = crater_x, cal = FALSE, cal_par_names = character(),
  ncalp = 0, theta_names = character()
)
check_jacobian(
  "IYDT crater calibration fixture", f_c, g_c, crater_beta, crater_params
)

crater_new <- datasets[["crater_new.csv"]]
crater_event <- list(
  beta = crater_beta[1], theta_names = "W", yield_scaling = 1 / 3,
  X = as.matrix(crater_new[, "HOB", drop = FALSE])
)
check_jacobian(
  "IYDT crater event fixture", f0_c, g0_c, crater_new$W[[1]], crater_event
)

stopifnot(identical(
  c(yield = 1 / 3, pressure = 1 / 3, temperature = 1 / 2),
  c(yield = seismic_params$yield_scaling,
    pressure = acoustic_params$pressure_scaling,
    temperature = acoustic_params$temp_scaling)
))

assert_close(
  f_s(seismic_beta, seismic_params),
  c(-0.39724995674087382, -0.40470098604087412),
  "IYDT seismic golden calibration output", 1e-12, 1e-12
)
assert_close(
  f0_s(unlist(seismic_new[1, c("W", "HOB")], use.names = FALSE), seismic_event),
  c(4.4575219949866209, 4.1821584038866204),
  "IYDT seismic golden event output", 1e-12, 1e-12
)
assert_close(
  f_a(acoustic_beta, acoustic_params),
  c(-0.45964296561385226, 0.033000719186147554),
  "IYDT acoustic golden calibration output", 1e-12, 1e-12
)
assert_close(
  f0_a(unlist(acoustic_new[1, c("W", "HOB")], use.names = FALSE), acoustic_event),
  c(4.2604876958481226, 4.1244613034481228),
  "IYDT acoustic golden event output", 1e-12, 1e-12
)
assert_close(
  f_o(optical_beta, optical_params),
  c(3.8525405118094564, 3.8231319507982717),
  "IYDT optical golden calibration output", 1e-12, 1e-12
)
assert_close(
  f0_o(unlist(optical_new[1, c("W", "HOB")], use.names = FALSE), optical_event),
  3.4123832725916197,
  "IYDT optical golden event output", 1e-12, 1e-12
)
assert_close(
  f_c(crater_beta, crater_params),
  c(6.0686274386201333, 6.3306037432904327),
  "IYDT crater golden calibration output", 1e-12, 1e-12
)
assert_close(
  f0_c(crater_new$W[[1]], crater_event), 5.1659440382527331,
  "IYDT crater golden event output", 1e-12, 1e-12
)

message("IYDT golden forward outputs passed")
