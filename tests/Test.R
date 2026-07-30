require(pkr)
tblNCA(Theoph, "Subject", "Time", "conc", dose=320, concUnit="mg/L")
tblNCA(Theoph, "Subject", "Time", "conc", dose=320, down="Log", concUnit="mg/L")
tblNCA(Indometh, "Subject", "time", "conc", dose=25, adm="Bolus", concUnit="mg/L")
tblNCA(Indometh, "Subject", "time", "conc", dose=25, adm="Bolus", down="Log", concUnit="mg/L")
tblNCA(Indometh, "Subject", "time", "conc", dose=25, adm="Infusion", dur=0.25, concUnit="mg/L")
tblNCA(Indometh, "Subject", "time", "conc", dose=25, adm="Infusion", dur=0.25, down="Log", concUnit="mg/L")
tblNCA(Indometh, "Subject", "time", "conc", dose=25, concUnit="mg/L")
tblNCA(Indometh, "Subject", "time", "conc", dose=25, down="Log", concUnit="mg/L")

# The core computation is delegated to NonCompart; check that it is passed through
# unchanged rather than re-implemented here.
stopifnot(all.equal(
  tblNCA(Theoph, "Subject", "Time", "conc", dose=320, concUnit="mg/L"),
  NonCompart::tblNCA(Theoph, "Subject", "Time", "conc", dose=320, concUnit="mg/L"),
  check.attributes = FALSE))

# Figure generation. The cases below are the ones that used to fail or to lose
# data silently: more series than the old 18-colour vector, more treatments than
# fit in one row of panels, and a missing concentration.
Out = file.path(tempdir(), "pkrTest")
Sim = function(nSubj, nTrt) {
  set.seed(1)
  Time = c(0, 0.25, 0.5, 1, 2, 3, 4, 6, 8, 12, 24)
  do.call(rbind, lapply(seq_len(nSubj), function(i) {
    do.call(rbind, lapply(seq_len(nTrt), function(j) {
      ka = 1.2*exp(rnorm(1, 0, 0.3)) ; ke = 0.12*exp(rnorm(1, 0, 0.3)) ; V = 30*exp(rnorm(1, 0, 0.2))
      data.frame(ID = sprintf("S%03d", i), TRT = paste0("TRT", j), Time = Time,
                 conc = pmax(100/V*ka/(ka - ke)*(exp(-ke*Time) - exp(-ka*Time)), 0),
                 stringsAsFactors = FALSE)
    }))
  }))
}

nDev = length(dev.list())
Wide = Sim(25, 1)
stopifnot(length(plotPK(Wide, "ID", "Time", "conc", outdir=Out, name="Wide")) == 5)

Cross = Sim(6, 8)
stopifnot(length(plotPK(Cross, "ID", "Time", "conc", trt="TRT", outdir=Out, name="Cross")) == 5)

Miss = Sim(6, 1) ; Miss[3, "conc"] = NA
stopifnot(length(plotPK(Miss, "ID", "Time", "conc", outdir=Out, name="Miss")) == 5)

# No device may be left open, and none of the figures may be empty.
stopifnot(length(dev.list()) == nDev)
stopifnot(all(file.size(list.files(Out, full.names=TRUE)) > 0))

pdfNCA(file.path(Out, "Theoph.pdf"), Theoph, "Subject", "Time", "conc", dose=320, concUnit="mg/L")
stopifnot(length(dev.list()) == nDev)

unlink(Out, recursive=TRUE)

# ---------------------------------------------------------------- rNCA (SDTM)
# rNCA had no test coverage at all, which is how it went from 0.1.0 to 0.1.5
# unable to run on the data.frames loadEXPC() returns.
DTC = function(base, h) format(as.POSIXct(base, tz="UTC") + h*3600, "%Y-%m-%dT%H:%M:%S")
Conc = function(tt, D=100, ka=1.1, ke=0.15, V=30) round(D/V*ka/(ka - ke)*(exp(-ke*tt) - exp(-ka*tt)), 3)
Base = "2024-01-01 08:00:00"
Tp = c(0, 0.5, 1, 2, 4, 6, 8, 12, 24)

EXrec = function(subj, start, route="ORAL", end=start, trt="DRUGA")
  data.frame(STUDYID="ST01", USUBJID=subj, EXTRT=trt, EXDOSE="100", EXDOSU="mg",
             EXROUTE=route, EXSTDTC=start, EXENDTC=end, stringsAsFactors=FALSE)
PCrec = function(subj, dtc, conc, testcd="DRUGA")
  data.frame(STUDYID="ST01", USUBJID=subj, PCTESTCD=testcd, PCDTC=dtc,
             PCSTRESN=as.character(conc), PCSTRESC=as.character(conc), PCSTRESU="ug/L",
             PCSPEC="PLASMA", PCLLOQ="0", stringsAsFactors=FALSE)

# A plain single-dose study must simply work, from data.frames.
EX = EXrec("S-01", DTC(Base, 0))
PC = PCrec("S-01", DTC(Base, Tp), Conc(Tp))
R1 = rNCA(EX, PC)
stopifnot(is.data.frame(R1), nrow(R1) == 1, R1$AdmMethod == "Extravascular")
stopifnot(all(c("STUDYID","USUBJID","EXTRT","PCTESTCD","PCRFTDTC","Dose","CMAX","AUCLST") %in% names(R1)))
# It must agree with a direct single-profile analysis of the same numbers.
Direct = sNCA(Tp, Conc(Tp), dose=100, concUnit="ug/L")
stopifnot(isTRUE(all.equal(R1$AUCLST, as.numeric(Direct["AUCLST"]), tolerance=1e-6)))

# A matrix must give the same answer as a data.frame.
stopifnot(isTRUE(all.equal(R1, rNCA(as.matrix(EX), as.matrix(PC)), check.attributes=FALSE)))

# The result must not depend on the row order of either domain.
Shuf = function(d, s) { set.seed(s); d[sample(nrow(d)), , drop=FALSE] }
stopifnot(isTRUE(all.equal(R1, rNCA(Shuf(EX,1), Shuf(PC,2)), check.attributes=FALSE)))

# Multiple dosing: one row per interval, and shuffling still changes nothing.
EX2 = rbind(EXrec("S-01", DTC(Base, 0)), EXrec("S-01", DTC(Base, 24)))
PC2 = PCrec("S-01", DTC(Base, c(Tp, 24 + Tp[-1])), c(Conc(Tp), Conc(Tp[-1]) + 3))
R2 = rNCA(EX2, PC2)
stopifnot(nrow(R2) == 2)
stopifnot(isTRUE(all.equal(R2, rNCA(Shuf(EX2,3), Shuf(PC2,4)), check.attributes=FALSE)))

# Crossover: period 1 must stop at the period-2 dose, not absorb it.
EXx = rbind(EXrec("S-01", DTC(Base, 0), trt="DRUGA"),
            EXrec("S-01", DTC(Base, 168), trt="DRUGB"))
PCx = PCrec("S-01", DTC(Base, c(Tp, 168 + Tp)), c(Conc(Tp), Conc(Tp)*1.3))
Rx = rNCA(EXx, PCx, trt="DRUGA")
stopifnot(nrow(Rx) == 1, Rx$TLST == 24)
stopifnot(isTRUE(all.equal(Rx$AUCLST, as.numeric(Direct["AUCLST"]), tolerance=1e-6)))

# Mixing bolus and infusion must not shift values into the wrong columns.
Iv = c(0, 0.25, 0.5, 1, 2, 4, 6, 8, 12, 24)
EXm = rbind(EXrec("S-01", DTC(Base,0), route="INTRAVENOUS"),
            EXrec("S-02", DTC(Base,0), route="INTRAVENOUS", end=DTC(Base,1)))
PCm = rbind(PCrec("S-01", DTC(Base,Iv), round(4*exp(-0.2*Iv),3)),
            PCrec("S-02", DTC(Base,Iv), round(4*exp(-0.2*Iv),3)))
Rm = rNCA(EXm, PCm)
stopifnot(nrow(Rm) == 2, identical(Rm$AdmMethod, c("Bolus","Infusion")))
stopifnot(is.na(Rm$C0[2]), all(!is.na(Rm$VSSO)))   # C0 exists for bolus only

# An oral record carrying an end date/time stays extravascular.
Ro = rNCA(EXrec("S-01", DTC(Base,0), end=DTC(Base,1)), PC)
stopifnot(Ro$AdmMethod == "Extravascular", "VZFO" %in% names(Ro), !("VSSO" %in% names(Ro)))

# Case-insensitive analyte selection, and the label comes from the data.
Rc = rNCA(EX, PC, analyte="druga")
stopifnot(nrow(Rc) == 1, Rc$PCTESTCD == "DRUGA")

# One unreadable concentration must not abort the run.
PCbad = PC ; PCbad$PCSTRESN[3] = "<0.1"
stopifnot(nrow(suppressWarnings(rNCA(EX, PCbad))) == 1)

# Nothing qualifying is a warning and NULL, not an error.
stopifnot(is.null(suppressWarnings(rNCA(EX, PC[1:4, , drop=FALSE]))))

# A placebo dose delivers no drug, so it must not close an active interval.
EXdd = rbind(EXrec("S-01", DTC(Base,0), trt="DRUGA"),
             data.frame(STUDYID="ST01", USUBJID="S-01", EXTRT="PLACEBO", EXDOSE="0", EXDOSU="mg",
                        EXROUTE="ORAL", EXSTDTC=DTC(Base,12), EXENDTC=DTC(Base,12),
                        stringsAsFactors=FALSE))
stopifnot(isTRUE(all.equal(suppressWarnings(rNCA(EXdd, PC, trt="DRUGA"))$AUCLST, R1$AUCLST)))

# Neither must a dose recorded under a study the caller excluded.
EXs2 = EXrec("S-01", DTC(Base,6)) ; EXs2$STUDYID = "ST02"
stopifnot(isTRUE(all.equal(suppressWarnings(rNCA(rbind(EX, EXs2), PC, study="ST01"))$AUCLST, R1$AUCLST)))

# A blank EXROUTE is extravascular, so mixing it with intravenous must be refused.
EXbr = EXrec("S-01", DTC(Base,0), route="INTRAVENOUS", end=DTC(Base,1))
EXbr2 = EXrec("S-02", DTC(Base,0)) ; EXbr2$EXROUTE = ""
PCbr = rbind(PCrec("S-01", DTC(Base,Tp), Conc(Tp)), PCrec("S-02", DTC(Base,Tp), Conc(Tp)))
stopifnot(inherits(try(rNCA(rbind(EXbr, EXbr2), PCbr), silent=TRUE), "try-error"))

# Duplicate sampling and dosing times must still give an order-independent result.
PCdup = rbind(PC, PCrec("S-01", DTC(Base,0), 9.99))
Rdup = suppressWarnings(rNCA(EX, PCdup))
for (s in 1:5) {
  set.seed(s)
  stopifnot(isTRUE(all.equal(Rdup,
    suppressWarnings(rNCA(EX, PCdup[sample(nrow(PCdup)), , drop=FALSE])), check.attributes=FALSE)))
}

# Selections match case-insensitively for every selector, not just some.
stopifnot(nrow(suppressWarnings(rNCA(rbind(EX, EXrec("S-02", DTC(Base,0))),
                                    rbind(PC, PCrec("S-02", DTC(Base,Tp), Conc(Tp))),
                                    id=c("S-01","s-02")))) == 2)

# excludeInfusion: the terminal slope must not be taken from samples drawn while
# the infusion is still running, and that is the DEFAULT. The profile below is
# already log-linear from 8 h, two hours before the infusion ends, so the
# unrestricted search reaches back into it and reports 8 h as the terminal phase.
Tinf = c(0, 0.5, 1, 2, 4, 6, 8, 9.5, 10.5, 11, 12, 14, 16, 24)
Cinf = c(0, 1200, 1600, 1400, 1200, 1100, round(1000*exp(-0.30*(Tinf[Tinf >= 8] - 8)), 2))
EXinf = EXrec("S-01", DTC(Base, 0), route="INTRAVENOUS", end=DTC(Base, 10))
PCinf = PCrec("S-01", DTC(Base, Tinf), Cinf)
Roff = suppressWarnings(rNCA(EXinf, PCinf, MinPoints=3, excludeInfusion=FALSE))
Ron  = suppressWarnings(rNCA(EXinf, PCinf, MinPoints=3, excludeInfusion=TRUE))
Rdef = suppressWarnings(rNCA(EXinf, PCinf, MinPoints=3))
stopifnot(Roff$AdmMethod == "Infusion", Roff$Dur == 10)
stopifnot(Roff$LAMZLL < 10)      # what the unrestricted search does
stopifnot(Ron$LAMZLL  > 10)      # what the restriction gives
stopifnot(Ron$LAMZNPT < Roff$LAMZNPT)
# The default must be the restricted behaviour.
stopifnot(isTRUE(all.equal(Rdef, Ron, check.attributes=FALSE)))
stopifnot(Rdef$LAMZLL > Rdef$Dur)
# The option may only move the terminal slope. Every observed parameter - one not
# computed from lambda z - must be bit-identical, and nothing outside the
# lambda-z-derived set may move at all. AUCIFO and the other extrapolated-to-
# infinity quantities are lambda-z derived by definition and so are expected to
# differ; AUClast, AUCall and AUMClast are observed and are not.
Observed = c("CMAX","CMAXD","TMAX","TLAG","CLST","TLST","AUCLST","AUCALL","AUMCLST",
             "C0","MRTIVLST","MRTEVLST","Dose","Dur")
LamzDerived = c("b0","CLSTP","LAMZ","LAMZHL","LAMZLL","LAMZUL","LAMZNPT","CORRXY","R2","R2ADJ",
                "AUCIFO","AUCIFOD","AUCIFP","AUCIFPD","AUCPEO","AUCPEP","AUCPBEO","AUCPBEP",
                "AUMCIFO","AUMCIFP","AUMCPEO","AUMCPEP","VZO","VZP","VZFO","VZFP",
                "CLO","CLP","CLFO","CLFP","VSSO","VSSP",
                "MRTIVIFO","MRTIVIFP","MRTEVIFO","MRTEVIFP")
Num = intersect(names(Roff), names(Ron))
Num = Num[vapply(Num, function(n) is.numeric(Roff[[n]]) && is.numeric(Ron[[n]]), logical(1))]
Moved = Num[!vapply(Num, function(n) isTRUE(all.equal(Roff[[n]], Ron[[n]], tolerance=1e-12)), logical(1))]
stopifnot(length(intersect(Moved, Observed)) == 0)   # no observed parameter moved
stopifnot(length(setdiff(Moved, LamzDerived)) == 0)  # nothing unexpected moved
stopifnot("LAMZLL" %in% Moved)                       # and the slope window really did

# It is a no-op where there is no infusion to exclude, so the default cannot
# change an extravascular or a bolus result.
EXev = EXrec("S-01", DTC(Base, 0))
stopifnot(isTRUE(all.equal(suppressWarnings(rNCA(EXev, PC, excludeInfusion=FALSE)),
                           suppressWarnings(rNCA(EXev, PC, excludeInfusion=TRUE)),
                           check.attributes=FALSE)))
EXbo = EXrec("S-01", DTC(Base, 0), route="INTRAVENOUS")
stopifnot(isTRUE(all.equal(suppressWarnings(rNCA(EXbo, PC, excludeInfusion=FALSE)),
                           suppressWarnings(rNCA(EXbo, PC, excludeInfusion=TRUE)),
                           check.attributes=FALSE)))
# sNCA stays a faithful NonCompart wrapper - the restriction lives in NCA0/rNCA.
stopifnot(isTRUE(all.equal(
  sNCA(Tinf, Cinf, dose=1000, adm="Infusion", dur=10, concUnit="ng/mL", returnNA=FALSE),
  NonCompart::sNCA(Tinf, Cinf, dose=1000, adm="Infusion", dur=10, concUnit="ng/mL"),
  check.attributes=FALSE)))

# Too few post-infusion samples: warn and fall back rather than silently honour
# nothing, and never drop the profile.
Wi = character(0)
Rfew = withCallingHandlers(rNCA(EXinf, PCinf[Tinf <= 12, , drop=FALSE], MinPoints=3, excludeInfusion=TRUE),
        warning=function(w) { Wi <<- c(Wi, conditionMessage(w)); invokeRestart("muffleWarning") })
stopifnot(any(grepl("post-infusion", Wi)), nrow(Rfew) == 1)

# A pre-dose sample is taken from the dosing day, not from the previous period.
EXcr = rbind(EXrec("S-01", DTC(Base,0), trt="DRUGA"), EXrec("S-01", DTC(Base,72), trt="DRUGB"))
Late = c(6,8,10,12,16,24)
PCcr = rbind(PCrec("S-01", DTC(Base,c(0,0.5,1,1.5,2)), Conc(c(0,0.5,1,1.5,2))),
             PCrec("S-01", DTC(Base,72+Late), Conc(Late)))
Rcr = suppressWarnings(rNCA(EXcr, PCcr, MinPoints=3))
Rb = Rcr[Rcr$EXTRT == "DRUGB", ]
stopifnot(nrow(Rb) == 1, Rb$TMAX == 6)   # not 0, which a period-1 sample would give
stopifnot(isTRUE(all.equal(Rb$AUCLST,
          as.numeric(sNCA(Late, Conc(Late), dose=100, concUnit="ug/L")["AUCLST"]), tolerance=1e-6)))
