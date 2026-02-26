library(dplyr)
library(stringr)
library(MatchIt)
library(SummarizedExperiment)
library(tibble)
library(phyloseq)
library(curatedMetagenomicData)

###### 1.Read all samples ######
meta_all <- curatedMetagenomicData::sampleMetadata

###### 2. filter ######
mk_region_df <- function(meta_all, cols, keep_only_stool = TRUE) {
  df <- meta_all
  
  if (keep_only_stool && !is.na(cols$body) && cols$body %in% names(df)) {
    bs <- tolower(df[[cols$body]])
    df <- df[grepl("^stool$|^fec(es)?$|^faec(es)?$", bs), , drop = FALSE]
  }
  
  
  pick <- function(x) if (!is.na(x) && x %in% names(df)) df[[x]] else NULL
  pick_num <- function(x) if (!is.na(x) && x %in% names(df)) suppressWarnings(as.numeric(df[[x]])) else NA_real_
  
  
  sampleID   <- pick(cols$id)
  subject_id <- pick(cols$subject)
  study_name <- pick(cols$study)
  
  country_std <- if (!is.na(cols$country) && cols$country %in% names(df)) df[[cols$country]] else NA
  
  BMI    <- pick_num(cols$bmi)
  age    <- pick_num(cols$age)
  sex    <- pick(cols$sex)
  sex_std <- tolower(trimws(as.character(if (!is.null(sex)) sex else NA)))
  sex_std[grepl("^m(ale)?$",  sex_std)] <- "Male"
  sex_std[grepl("^f(emale)?$", sex_std)] <- "Female"
  sex_std[sex_std %in% c("", "na", "nan", "missing")] <- NA
  sex_std[!is.na(sex_std) & !(sex_std %in% c("Male","Female"))] <- "Other"
  sex_std[is.na(sex_std)] <- "Unknown"
  sex_std <- factor(sex_std, levels = c("Male","Female","Other","Unknown"))
  
  height <- pick_num(cols$height)
  weight <- pick_num(cols$weight)
  
  antibiotics         <- pick(cols$antibiotics)
  antibiotics_family  <- pick(cols$antibiotics_family)
  alcohol             <- pick(cols$alcohol)
  smoking             <- pick(cols$smoking)
  diet                <- pick(cols$diet)
  
  disease_vec <- NULL
  for (cand in c(cols$disease, cols$study_condition, cols$health_status, "disease","study_condition","health_status","status")) {
    if (!is.na(cand) && cand %in% names(df)) { disease_vec <- df[[cand]]; break }
  }
  
  birth_weight     <- pick_num(cols$birth_weight)
  gestational_age  <- pick_num(cols$gestational_age)
  premature        <- pick(cols$premature)
  feeding_practice <- pick(cols$feeding_practice)
  born_method      <- pick(cols$born_method)
  pregnant         <- pick(cols$pregnant)
  lactating        <- pick(cols$lactating)
  
  travel_destination <- pick(cols$travel_destination)
  
  non_westernized <- pick(cols$non_western)
  
  hscrp      <- pick_num(cols$hscrp)
  esr        <- pick_num(cols$esr)
  creatinine <- pick_num(cols$creatinine)
  albumin    <- pick_num(cols$albumin)
  
  out <- dplyr::tibble(
    sampleID   = if (!is.null(sampleID))   sampleID   else NA,
    subject_id = if (!is.null(subject_id)) subject_id else NA,
    study_name = if (!is.null(study_name)) study_name else NA,
    country_std = country_std,   
    
    BMI   = BMI,
    age   = age,
    sex   = sex_std,
    height = height,
    weight = weight,
    
    antibiotics        = if (!is.null(antibiotics))        antibiotics        else NA,
    antibiotics_family = if (!is.null(antibiotics_family)) antibiotics_family else NA,
    alcohol            = if (!is.null(alcohol))            alcohol            else NA,
    smoking            = if (!is.null(smoking))            smoking            else NA,
    diet               = if (!is.null(diet))               diet               else NA,
    
    disease_status = if (!is.null(disease_vec)) disease_vec else NA,
    
    birth_weight      = birth_weight,
    gestational_age   = gestational_age,
    premature         = if (!is.null(premature))        premature        else NA,
    feeding_practice  = if (!is.null(feeding_practice)) feeding_practice else NA,
    born_method       = if (!is.null(born_method))      born_method      else NA,
    pregnant          = if (!is.null(pregnant))         pregnant         else NA,
    lactating         = if (!is.null(lactating))        lactating        else NA,
    
    travel_destination = if (!is.null(travel_destination)) travel_destination else NA,
    
    non_westernized = if (!is.null(non_westernized)) non_westernized else NA,
    
    hsCRP      = hscrp,
    ESR        = esr,
    creatinine = creatinine,
    albumin    = albumin
  ) %>%
    dplyr::filter(!is.na(sampleID)) %>%
    dplyr::distinct(sampleID, .keep_all = TRUE)%>%   #######
  dplyr::filter(is.na(pregnant)|pregnant!="yes")%>%
    dplyr::filter(is.na(alcohol)|alcohol!="yes") %>%
    dplyr::filter(is.na(smoking)|smoking!="yes")
  
  out
}


detect_cols <- function(df) {
  nm <- names(df)
  
  find1 <- function(patterns) {
    hits <- unique(unlist(lapply(patterns, function(p) grep(p, nm, ignore.case = TRUE, value = TRUE))))
    if (length(hits) == 0) return(NA_character_)
    hits[1]
  }
  
  find_all <- function(patterns) {
    unique(unlist(lapply(patterns, function(p) grep(p, nm, ignore.case = TRUE, value = TRUE))))
  }
  
  list(
    id            = find1(c("^sample[_ ]?id$", "^sample$", "^SampleID$")),
    subject       = find1(c("^subject[_ ]?id$", "^subject$")),
    study         = find1(c("^study[_ ]?name$", "^study$", "cohort")),
    bmi           = find1(c("^bmi$", "bmi_?kg", "body.?mass")),
    age           = find1(c("^age$", "host[_ ]?age")),
    age_category  = find1(c("^age[_ ]?category$")),
    sex           = find1(c("^sex$", "^gender$")),
    height        = find1(c("^height$")),
    weight        = find1(c("^weight$")),
    body          = find1(c("^body[_ ]?site$", "^Body[_ ]?site$")),
    country       = find1(c("^country$", "host[_ ]?country", "study[_ ]?country")),
    location      = find1(c("^location$", "geo[_ ]?loc", "geo[_ ]?loc[_ ]?name", "city", "state|province|region|continent")),
    non_western   = find1(c("^non[_ ]?westernized$", "westernization", "urban|rural")),
    disease       = find1(c("^disease$", "diagnosis", "phenotype")),
    study_condition = find1(c("^study[_ ]?condition$", "^condition$")),
    health_status = find1(c("^health[_ ]?status$", "^status$")),
    antibiotics   = find1(c("antibiotic", "antibiotics_current_use", "recent[_ ]?antibiotics", "antibiotics")),
    antibiotics_family = find1(c("antibiotics[_ ]?family")),
    alcohol       = find1(c("^alcohol(_numeric)?$", "alcohol[_ ]?use", "drinks?")),
    smoking       = find1(c("^smok(er|ing)(_status)?$", "brinkman", "^smoker$")),
    diet          = find1(c("^diet$", "dietary", "food", "fiber")),
    birth_weight      = find1(c("^birth[_ ]?weight$")),
    gestational_age   = find1(c("^gestational[_ ]?age$", "^ga$")),
    premature         = find1(c("^premature$", "preterm")),
    feeding_practice  = find1(c("^feeding[_ ]?practice$", "breast|formula|weaning")),
    born_method       = find1(c("born[_ ]?method", "delivery|caesarean|cesarean|vaginal")),
    pregnant          = find1(c("^pregnan")),
    lactating         = find1(c("^lactat")),
    travel_destination = find1(c("^travel[_ ]?destination$", "travel")),
    hscrp          = find1(c("^hs[_-]?crp$", "^crp$", "high[_ ]?sensitivity[_ ]?crp")),
    esr            = find1(c("^esr$")),
    creatinine     = find1(c("^creatin(e|ine)$", "^creat$")),
    albumin        = find1(c("^albumin(e)?$")),
    
    
    covariates_available = unique(c(
      find_all(c("^age$", "host[_ ]?age", "^age[_ ]?category$")),
      find_all(c("^sex$", "^gender$")),
      find_all(c("^height$", "^weight$", "^bmi$", "body.?mass")),
      find_all(c("antibiotic", "antibiotics_current_use", "recent[_ ]?antibiotics", "antibiotics[_ ]?family")),
      find_all(c("^alcohol(_numeric)?$", "alcohol[_ ]?use", "drinks?")),
      find_all(c("^smok(er|ing)(_status)?$", "brinkman", "^smoker$")),
      find_all(c("^diet$", "dietary", "food", "fiber")),
      find_all(c("^health[_ ]?status$", "^status$", "^study[_ ]?condition$", "^condition$", "^disease$")),
      find_all(c("^country$", "host[_ ]?country", "study[_ ]?country", "^location$", "geo[_ ]?loc", "geo[_ ]?loc[_ ]?name",
                 "city", "state|province|region|continent", "^non[_ ]?westernized$", "westernization", "urban|rural")),
      find_all(c("^birth[_ ]?weight$", "^gestational[_ ]?age$", "^ga$", "^premature$", "preterm",
                 "^feeding[_ ]?practice$", "breast|formula|weaning", "born[_ ]?method", "delivery|caesarean|cesarean|vaginal",
                 "^pregnan", "^lactat")),
      find_all(c("^travel[_ ]?destination$", "travel")),
      find_all(c("^hs[_-]?crp$", "^crp$", "^esr$", "^creatin(e|ine)$", "^creat$", "^albumin(e)?$"))
    ))
  )
}

cols <- detect_cols(meta_all)

df0 <- mk_region_df(meta_all, cols)


###### 3.PSM ######
do_psm_binary <- function(
    df, treat_var, treat_ref, treat_cmp,
    y_var = "BMI",
    covars = NULL,                
    id_var = "sampleID",
    use_overlap_trim = TRUE,
    trim_q = c(0.05, 0.95),
    exact_vars = "sex",
    caliper_sd = 0.1,
    verbose = TRUE
) {
  stopifnot(all(c(treat_var, id_var) %in% names(df)))
  msg <- function(...) if (verbose) message(...)
  
  dat <- df %>%
    dplyr::mutate(.treat = dplyr::case_when(
      .data[[treat_var]] %in% treat_cmp ~ 1L,
      .data[[treat_var]] %in% treat_ref ~ 0L,
      TRUE ~ NA_integer_
    )) %>%
    dplyr::filter(!is.na(.treat), !is.na(.data[[id_var]]))
  
  tab0 <- table(dat$.treat)
  if (length(tab0) < 2 || any(tab0 < 2)) {
    msg(sprintf("stop",
                unname(tab0[as.character(0)]), unname(tab0[as.character(1)])))
    return(NULL)
  }
  msg("base:", nrow(dat), ";0/1=", paste(tab0, collapse = "/"))
  
  covars <- covars[covars %in% names(dat)]
  
  use_cov <- covars[vapply(covars, function(v) {
    x <- dat[[v]]
    any(!is.na(x)) && (length(unique(na.omit(x))) > 1)
  }, logical(1))]
  
  for (v in use_cov) {
    if (is.character(dat[[v]]) || is.logical(dat[[v]])) dat[[v]] <- as.factor(dat[[v]])
  }
  
  vars_for_ps <- c(".treat", use_cov)
  mf <- dat[, vars_for_ps, drop = FALSE]
  cc <- stats::complete.cases(mf)
  mf_cc <- mf[cc, , drop = FALSE]
  
  is_bad <- vapply(names(mf_cc), function(v) {
    if (v == ".treat") return(FALSE)
    x <- mf_cc[[v]]
    if (is.factor(x)) {
      nlevels(droplevels(x)) < 2
    } else if (is.numeric(x)) {
      sd(x, na.rm = TRUE) == 0 || length(unique(na.omit(x))) < 2
    } else {
      FALSE
    }
  }, logical(1))
  
  use_cov2 <- setdiff(use_cov, names(mf_cc)[is_bad])
  
  
  form <- if (length(use_cov2)) {
    as.formula(paste(".treat ~", paste(use_cov2, collapse = " + ")))
  } else {
    .treat ~ 1
  }
  
  ps_fit <- glm(form, data = dat, family = binomial(), na.action = na.exclude,
                control = list(maxit = 100))
  ps_hat <- predict(ps_fit, newdata = dat, type = "response")
  
  eps <- 1e-6
  ps_hat <- pmin(pmax(ps_hat, eps), 1 - eps)
  
  if (use_overlap_trim) {
    keep_finite <- is.finite(ps_hat) & !is.na(ps_hat)
    if (!all(keep_finite)) msg(sprintf("remove infinite PS:%d", sum(!keep_finite)))
    dat    <- dat[keep_finite, , drop = FALSE]
    ps_hat <- ps_hat[keep_finite]
    
    tr <- dat$.treat
    if (length(unique(tr)) < 2L) { msg(" "); return(NULL) }
    
    ps_T <- ps_hat[tr == 1L]; ps_C <- ps_hat[tr == 0L]
    lo <- max(min(ps_T, na.rm = TRUE), min(ps_C, na.rm = TRUE))
    hi <- min(max(ps_T, na.rm = TRUE), max(ps_C, na.rm = TRUE))
    if (!is.finite(lo) || !is.finite(hi) || lo >= hi) {
      msg(" "); return(NULL)
    }
    keep_overlap <- (ps_hat >= lo) & (ps_hat <= hi)
    
    qT <- stats::quantile(ps_T, probs = trim_q, names = FALSE, na.rm = TRUE)
    qC <- stats::quantile(ps_C, probs = trim_q, names = FALSE, na.rm = TRUE)
    lo2 <- max(qT[1], qC[1]); hi2 <- min(qT[2])
    keep_trim <- if (!is.finite(lo2) || !is.finite(hi2) || lo2 >= hi2) {
      rep(TRUE, length(ps_hat))
    } else (ps_hat >= lo2) & (ps_hat <= hi2)
    
    keep <- keep_overlap & keep_trim & is.finite(ps_hat)
    if (sum(!keep) > 0) {
      # msg(sprintf("Remove %d samples (overlap:[%.3f,%.3f]; trim:[%.3f,%.3f])",
      #             sum(!keep), lo, hi, lo2, hi2))
    }
    dat    <- dat[keep, , drop = FALSE]
    ps_hat <- ps_hat[keep]
  }
  
  
  tab1 <- table(dat$.treat)
  if (length(tab1) < 2 || any(tab1 < 2)) {
    msg(sprintf("",
                unname(tab1[as.character(0)]), unname(tab1[as.character(1)])))
    return(NULL)
  }
  
  
  dat$ps_logit <- qlogis(ps_hat)
  sd_logit <- stats::sd(dat$ps_logit, na.rm = TRUE)
  cal_val <- if (is.finite(sd_logit) && sd_logit > 0) caliper_sd * sd_logit else {
    msg(""); 0.2
  }
  
  
  exact_use <- exact_vars[exact_vars %in% names(dat)]
  if (length(exact_use)) {
    for (v in exact_use) {
      if (!is.factor(dat[[v]])) dat[[v]] <- as.factor(dat[[v]])
      if (nlevels(droplevels(dat[[v]])) < 2) {
        msg(sprintf("exact one level", v))
      }
    }
    msg("exact=", paste(exact_use, collapse = ", "))
  }
  
  
  m.out <- MatchIt::matchit(
    .treat ~ 1, data = dat,
    method   = "nearest",
    distance = dat$ps_logit,
    caliper  = cal_val,
    ratio    = 1,
    replace  = FALSE,
    exact    = if (length(exact_use)) exact_use else NULL
  )
  matched <- MatchIt::match.data(m.out) %>% tibble::as_tibble()
  
  
  matched$..y <- if (y_var %in% names(dat)) matched[[y_var]] else NA_real_
  
  
  tab2 <- table(matched$.treat)
  msg(sprintf("complete : matched N=%d；0/1=%s/%s；caliper=%.3f",
              nrow(matched),
              unname(tab2[as.character(0)]), unname(tab2[as.character(1)]),
              cal_val))
  
  list(
    matched_df = matched %>%
      dplyr::transmute(
        !!id_var := .data[[id_var]],
        treat = as.integer(.treat),
        y = ..y
      ),
    matchit_object = m.out,
    settings = list(
      covars_used = use_cov2,
      caliper = cal_val,
      exact_vars = exact_use,
      overlap_trim = use_overlap_trim,
      trim_q = trim_q
    )
  )
}


mk_compare <- function(df0, ref, cmp, y_var = "BMI",
                       age_band_width = 1L,
                       caliper_sd = 0.1,
                       use_overlap_trim = TRUE,
                       trim_q = c(0.05, 0.95),
                       verbose = TRUE) {
  df_use <- df0 %>%
    dplyr::mutate(country_std = ifelse(.data$country_std %in% c(ref, cmp), .data$country_std, NA)) %>%
    dplyr::filter(!is.na(.data$country_std))
  
  if ("sex" %in% names(df_use)) {
    df_use$sex <- as.character(df_use$sex)
    df_use$sex[is.na(df_use$sex) | df_use$sex %in% c("", "NA", "nan", "missing")] <- "Unknown"
    df_use$sex <- factor(df_use$sex, levels = c("Male","Female","Other","Unknown"))
  }
  

  age_num <- suppressWarnings(as.numeric(df_use$age))
  
  low  <- floor(age_num / age_band_width) * age_band_width
  high <- low + (age_band_width - 1)
  age_band <- ifelse(is.finite(low), paste0(low, "-", high), NA)
  df_use$age_band <- factor(age_band)
  
  
  covars_ps <- c("antibiotics")   
  if ("age" %in% names(df_use)) covars_ps <- union("age", covars_ps)
  
  
  do_psm_binary(
    df = df_use,
    treat_var = "country_std",
    treat_ref = ref,
    treat_cmp = cmp,
    y_var = y_var,
    covars = covars_ps,
    id_var = "sampleID",
    exact_vars = c("sex","age_band","antibiotics"), #######for exact match
    caliper_sd = caliper_sd,
    use_overlap_trim = use_overlap_trim,
    trim_q = trim_q,
    verbose = verbose
  )
}


df_CHN <- df0 %>% filter(country_std == "CHN")
df_US <- df0 %>% filter(country_std == "USA")

res_cu <- mk_compare(df0, ref = "CHN", caliper_sd = 0.2,cmp = "USA", y_var = "BMI")

data_china_usa <- if (!is.null(res_cu)) list(
  treat = res_cu$matched_df$treat,
  y = res_cu$matched_df$y,
  sample_id = res_cu$matched_df$sampleID
) else NULL


data_china_usa <- if (!is.null(res_cu)) list(
  treat = res_cu$matched_df$treat,
  y = res_cu$matched_df$y,
  sample_id = res_cu$matched_df$sampleID
) else NULL


###### 4. get OTU ######

## China vs USA

matched_ids <- data_china_usa$sample_id   # or: data_china_uk$sample_id
id_col <- {
    cand <- c("sample_id","sampleID","SampleID","sample")
    cand[which(cand %in% names(meta_all))][1]
  }
stopifnot(!is.na(id_col))
  
  
meta_matched <- meta_all %>%
  semi_join(tibble(!!id_col := matched_ids), by = id_col)
  
  
tse <- returnSamples(
    sampleMetadata = meta_matched,
    dataType       = "relative_abundance",
    counts         = TRUE,      
    rownames       = "short"    
  )
  
current_ids <- colnames(tse) 
  
idx <- match(matched_ids, current_ids)
  
missing <- is.na(idx)
  
  
tse_matched <- tse[, idx, drop = FALSE]
stopifnot(all(colnames(tse_matched) == matched_ids))
  
count_mat <- as.matrix(assay(tse_matched))
colnames(count_mat) <- matched_ids   
  
  
keep <- match(matched_ids, data_china_usa$sample_id)
treat_vec <- data_china_usa$treat[keep]
y_vec     <- data_china_usa$y[keep]


count1 <- t(count_mat)
y1 <- y_vec
treat1 <-  treat_vec
count1 <- as.matrix(count1)
#count1 <- count1[,c(1:1000)]
  
all_data <- data.frame(count1, y1 = y1, treat1 = treat1)
  
keep_rows <- complete.cases(all_data)
# keep_rows <- keep_rows[c(c(1:200),c(1000:1200))]
  
count1 <- count1[keep_rows, ]
y1     <- y1[keep_rows]
treat1 <- treat1[keep_rows]
  
zero_prop <- colMeans(count1 == 0)  
count1 <- count1[, zero_prop <= 0.90]  

save(count1, y1, treat1, file = "real_data.RData")
