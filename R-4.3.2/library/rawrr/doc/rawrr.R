## ----style, echo = FALSE, results = 'asis'------------------------------------
BiocStyle::markdown()
knitr::opts_chunk$set(fig.wide = TRUE, fig.retina = 3, error=FALSE)

## ----HexSticker, echo=FALSE, out.width="50%", eval=TRUE-----------------------
  knitr::include_graphics("rawrr_logo.png")

## ----installAssemblies, echo=TRUE---------------------------------------------
if (isFALSE(rawrr::.checkDllInMonoPath())){
  rawrr::installRawFileReaderDLLs()
}
rawrr::installRawrrExe()

## ----tartareEH4547, warning=FALSE, message=FALSE, eval=TRUE-------------------
# fetch via ExperimentHub
library(ExperimentHub)
eh <- ExperimentHub::ExperimentHub()
EH4547 <- normalizePath(eh[["EH4547"]])

(rawfile <- paste0(EH4547, ".raw"))
if (!file.exists(rawfile)){
  file.copy(EH4547, rawfile)
}

## ----check, eval=TRUE, echo=FALSE---------------------------------------------
stopifnot(file.exists(rawfile))

## ----readFileHeader, eval=TRUE------------------------------------------------
H <- rawrr::readFileHeader(rawfile = rawfile)

## ----echo=FALSE---------------------------------------------------------------

H <- list("Instrument model"=NA, "Number of scans"=NA, "RAW file"=NA, "Time range"=c(NA,NA))


## -----------------------------------------------------------------------------
i <- c(9594, 11113, 11884, 12788, 12677, 13204, 13868, 14551, 16136, 17193, 17612)
S <- rawrr::readSpectrum(rawfile = rawfile, scan = i)
class(S)
class(S[[1]])
summary(S[[1]])
plot(S[[1]], centroid=TRUE)

## ----checkScan, echo=FALSE, message=FALSE, eval=FALSE-------------------------
#  # AGC-related
#  S[[1]]$`AGC:`
#  S[[1]]$`AGC Target:`
#  # Injection time-related
#  S[[1]]$`Ion Injection Time (ms)`
#  S[[1]]$`Max. Ion Time (ms):`
#  # Resolving power-related
#  S[[1]]$`FT Resolution:`
#  # HCD-related
#  S[[1]]$`HCD Energy:`
#  S[[1]]$`HCD Energy eV:`

## ----plotSN,  eval = TRUE, fig.cap = "Spectrum plot using signal-to-noise option. The vertical grey lines indicate the *in-silico* computed y-ions of the peptide precusor LGGNEQVTR++ as calculated by the [protViz]( https://CRAN.R-project.org/package=protViz) package [@protViz]. The horizontal blue line indicates a signal to noise ratio of 5. Diagnostic information related to AGC is shown in grey."----
plot(S[[1]], centroid=TRUE, SN = TRUE, diagnostic = TRUE)
# S/N threshold indicator
abline(h = 5, lty = 2, col = "blue")
# decorate plot with y-ion series of target peptide LGGNEQVTR++
(yIonSeries <- protViz::fragmentIon("LGGNEQVTR")[[1]]$y[1:8])
names(yIonSeries) <- paste0("y", seq(1, length(yIonSeries)))
abline(v = yIonSeries, col='#DDDDDD88', lwd=5)
axis(3, yIonSeries, names(yIonSeries))

## -----------------------------------------------------------------------------
# value casting
(x <- S[[1]]$scan)
class(x)
(y <- rawrr::scanNumber(S[[1]]))
class(y)
# error handling
S[[1]]$`FAIMS Voltage On:`
try(rawrr::faimsVoltageOn(S[[1]]))
# generating a well-behaved accessor
maxIonTime <- rawrr::makeAccessor(key = "Max. Ion Time (ms):", returnType = "double")
maxIonTime(S[[1]])

## ----TIC, fig.cap="Total ion chromatogram (TIC) calculated from all MS1-level scans contained in 20181113_010_autoQC01.raw."----
message(rawfile)
rawrr::readChromatogram(rawfile = rawfile, type = "tic") |>
  plot()

## ----BPC, fig.cap="Base peak chromatogram (BPC) calculated from all MS1-level scans contained in 20181113_010_autoQC01.raw."----
rawrr::readChromatogram(rawfile = rawfile, type = "bpc") |>
  plot()

## ----plotrawrrchromatogram, fig.cap="Extracted ion chromatograms (XIC) for iRT peptide precursors. Each XIC was calculated using a tolerance of 10 ppm around the target mass and using only MS1-level scans.", error=TRUE----
iRTmz <- c(487.2571, 547.2984, 622.8539, 636.8695, 644.8230, 669.8384, 683.8282,
            683.8541, 699.3388, 726.8361, 776.9301)

names(iRTmz) <- c("LGGNEQVTR", "YILAGVENSK", "GTFIIDPGGVIR", "GTFIIDPAAVIR",
                 "GAGSSEPVTGLDAK", "TPVISGGPYEYR", "VEATFGVDESNAK",
                 "TPVITGAPYEYR", "DGLDAASYYAPVR", "ADVTPADFSEWSK",
                 "LFLQFGAQGSPFLK")

C <- rawrr::readChromatogram(rawfile, mass = iRTmz, tol = 10, type = "xic", filter = "ms")
class(C)
plot(C, diagnostic = TRUE)

## ----fragmentIonTraces, fig.cap="Extracted ion chromatograms (XICs) for LGGNEQVTR++ fragment ions y6, y7, and y8. Each target ion (see legend) was extracted with a tolerance of 10 ppm from all scans matching the provided filter.", fig.retina=1----
rawrr::readChromatogram(rawfile = rawfile,
      mass = yIonSeries[c("y6", "y7", "y8")],
      type = 'xic',
      tol = 10,
      filter = "FTMS + c NSI Full ms2 487.2567@hcd27.00 [100.0000-1015.0000]") |>
  plot(diagnostic = TRUE)

## ----readIndex----------------------------------------------------------------
rawrr::readIndex(rawfile = rawfile) |>
  subset(scanType == "FTMS + c NSI Full ms2 487.2567@hcd27.00 [100.0000-1015.0000]") |>
  head()

## ----fittedAPEX, fig.cap="The plots show the fitted peak curves in red."------
par(mfrow = c(3, 4), mar=c(4,4,4,1))
rtFittedAPEX <- C |>
  rawrr:::pickPeak.rawrrChromatogram() |>
  rawrr:::fitPeak.rawrrChromatogram() |>
  lapply(function(x){
    plot(x$times, x$intensities, type='p',
         ylim=range(c(x$intensities,x$yp)),
         main=x$mass); lines(x$xx, x$yp,
                             col='red'); x}) |>
  sapply(function(x){x$xx[which.max(x$yp)[1]]})

## a simple alternative to derive rtFittedAPEX could be
rt <- sapply(C, function(x) x$times[which.max(x$intensities)[1]])

## ----iRTscoreFit, error=TRUE--------------------------------------------------
iRTscore <- c(-24.92, 19.79, 70.52, 87.23, 0, 28.71, 12.39, 33.38, 42.26, 54.62, 100)
fit <- lm(rtFittedAPEX ~ iRTscore)

## ----iRTscoreFitPlot, fig.small=TRUE, echo=FALSE, fig.cap="iRT regression. Plot shows observed peptide RTs as a function of iRT scores and fitted regression line of corresponding linear model obtained by ordinary least squares (OLS) regression."----
# iRTscoreFitPlot
plot(rtFittedAPEX ~ iRTscore,
     ylab = 'Retention time [min]',
     xlab = "iRT score",
     pch=16, frame.plot = FALSE)
abline(fit, col = 'grey')
abline(v = 0, col = "grey", lty = 2)
legend("topleft", legend = paste("Regression line: ", "rt =",
                                 format(coef(fit)[1], digits = 4), " + ",
                                 format(coef(fit)[2], digits = 2), "score",
                                 "\nR2: ", format(summary(fit)$r.squared, digits = 4)),
       bty = "n", cex = 0.75)
text(iRTscore, rt, iRTmz, pos=1,cex=0.5)

## ----faq1---------------------------------------------------------------------
# use sample.raw file
f <- rawrr::sampleFilePath()
(rawrr::readIndex(f) |>
  subset(MSOrder == "Ms"))$scan |>
  sapply(FUN = rawrr::readSpectrum, rawfile=f) |>
  parallel::mclapply(FUN = function(x){
    idx <- x$intensity == max(x$intensity);
    data.frame(scan=x$scan,
               StartTime=x$StartTime,
               mZ=x$mZ[idx],
               intensity=x$intensity[idx]
               )}) |>
  base::Reduce(f=rbind)

## ----sessioninfo, echo=FALSE--------------------------------------------------
sessionInfo()

