#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# CuratePurity.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Curate the purity of glioma samples using clonal alterations (IDHmt, 1p19q and optionally TP53)
#
# Author: Jurriaan Janssen (j.janssen.1@erasmusmc.nl)
#
# condaenv: PureCN
# Usage: 
#
# TODO:
# 1) 
#
# History:
#  13-07-2026: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
suppressMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))
suppressPackageStartupMessages(library(PureCN))
suppressMessages(library(dplyr))


#-------------------------------------------------------------------------------
# 0.2  Parse command line arguments
#-------------------------------------------------------------------------------
option_list <- list(
    make_option(c("--variants"), type = "character", default = NULL,
                help = "Input to old PureCN variants file"),
    make_option(c("--ACE_results"), type = "character", default = NULL,
        help = "Input to ACE results file"),
    make_option(c("--seg-file"), type = "character", default = NULL,
        help = "QDNAseq/CBS segmentation file"),
    make_option(c("--margin"), type = "double", default = 0.1,
        help = "+/- margin (in purity units) placed around the consensus estimate [default %default]"),
    make_option(c("--fallback-min-purity"), type = "double", default = 0.15,
        help = "Used verbatim when zero anchors are available [default %default]"),
    make_option(c("--fallback-max-purity"), type = "double", default = 0.95,
        help = "Used verbatim when zero anchors are available [default %default]"),
    make_option(c("--out"), type = "character", default = NULL,
        help = "Output CSV path")
)

opt <- parse_args(OptionParser(option_list = option_list), convert_hyphens_to_underscores = TRUE)


#opt$variants <- "/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/WES/PureCN/500kbp/SG_061/SG_061_tumor1_variants.csv"
#opt$ACE_results <- "/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/WES/ACE/500kbp/SG_061/ACE_fits.txt"
#opt$seg_file <- "/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/WES/QDNAseq/500kbp/SG_061/data/QDNAseq_Segments.txt"
#-------------------------------------------------------------------------------
# 1.1  Read files
#-------------------------------------------------------------------------------
variants <- read.csv(opt$variants)
ACE_results <- read.delim(opt$ACE_results)
segments <- read.delim(opt$seg_file)

#-------------------------------------------------------------------------------
# 1.2  Fetch purity measures
#-------------------------------------------------------------------------------
get_seg_mean <- function(chrom_q, pos, segments) {
    row <- segments %>%
        filter( chrom == chrom_q, loc.start <= pos, loc.end >= pos)
    if (nrow(row) == 0) return(NA_real_)
    row$seg.mean[1]
}

# 1. IDH1/2 VAF based purity
# Do not use recovered variants, allow different copy number states (usually gain)
IDH_LOCUS <- list(chr = "chr2", pos = 208248388)  # hg38

IDH_candidates <- variants %>%
    filter(gene.symbol %in% c('IDH1', 'IDH2'), ML.SOMATIC == TRUE, recovered == FALSE)

if (nrow(IDH_candidates) > 1) {
     flog.warn("Multiple IDH1/2 somatic calls found (%d); prioritizing IDH1 over IDH2, then highest depth.",
              nrow(IDH_candidates))
    IDH_candidates <- IDH_candidates %>%
        mutate(gene_priority = ifelse(gene.symbol == "IDH1", 1, 2)) %>%
        arrange(gene_priority, desc(depth)) %>%
        slice(1) %>%
        select(-gene_priority)
}
IDH_purity <- NA_real_
IDH_copy_ratio <- NA_real_


if (nrow(IDH_candidates) == 1) {
    IDH_copy_ratio <- get_seg_mean(IDH_candidates$chr[1], IDH_candidates$start[1], segments)
    diploid_band <- 0.08

    if (is.na(IDH_copy_ratio)) {
        flog.warn("No segmentation coverage at IDH locus; skipping as anchor.")
    } else if (abs(IDH_copy_ratio - 1) > diploid_band) {
        flog.warn("IDH locus not diploid (copy_ratio=%.3f); skipping as anchor -- likely a local/chromosome gain or loss.",
                   IDH_copy_ratio)
    } else {
        IDH_purity <- IDH_candidates$AR[1] * 2
    }
}


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# 2. TP53 VAF based purity
# Retrieve CN status: neutral, loss or cnLOH
TP53_variants <- variants %>%
    filter(gene.symbol == 'TP53', ML.SOMATIC == TRUE)

purity_from_row <- function(vaf, loh, cn_total = NULL, copy_ratio = NULL) {
    if (!loh) {
        return(2 * vaf)  # diploid het
    }
    # LOH == TRUE: distinguish deletion vs copy-neutral
    is_deletion <- (!is.null(cn_total) && cn_total == 1) ||
                    (!is.null(copy_ratio) && copy_ratio < 0.85)
    if (is_deletion) {
        return(2 * vaf / (1 + vaf)) # LOH via deletion
    }
    vaf  # copy-neutral LOH
}



if (nrow(TP53_variants) != 0) {
    TP53_annotated <- TP53_variants %>%
    rowwise() %>%
    mutate(copy_ratio = get_seg_mean(chr, start, segments)) %>%
    ungroup() %>%
    mutate(
        loh = case_when(
            is.na(copy_ratio)  ~ NA,
            copy_ratio < 0.85  ~ TRUE,    # deletion
            copy_ratio > 1.15  ~ NA,      # gain -- multiplicity ambiguous, exclude as anchor
            AR > 0.5           ~ TRUE,    # copy-neutral LOH-like state
            TRUE               ~ FALSE
        )
    ) %>%
    filter(!is.na(loh))
    
    purity_from_row_vec <- Vectorize(purity_from_row)

    TP53_scored <- TP53_annotated %>%
        mutate(purity = purity_from_row_vec(vaf = AR, loh = loh, copy_ratio = copy_ratio))

    if (nrow(TP53_scored) > 1) {
        flog.warn("Multiple TP53 somatic calls found (%d); using highest-depth VAF.", nrow(TP53_scored))
        TP53_scored <- TP53_scored %>% slice_max(AR, n = 1, with_ties = FALSE)
    }

    TP53_purity <- if (nrow(TP53_scored) == 1) TP53_scored$purity[1] else NA_real_
}else{
    TP53_purity <- NA
}


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# 3. 1p19q based purity
# Mean linear copy ratio across a chromosome arm, width-weighted.
# region: list(chrom=, start=, end=)
.arm_copy_ratio <- function(segments, region) {
    seg <- segments[segments$chrom == region$chrom &
                     segments$loc.start < region$end &
                     segments$loc.end   > region$start, ]
    if (nrow(seg) == 0) return(NA_real_)
    # clip segment widths to the arm boundary
    w_start <- pmax(seg$loc.start, region$start)
    w_end   <- pmin(seg$loc.end, region$end)
    width   <- w_end - w_start
    sum(seg$seg.mean * width) / sum(width)
}

# What fraction of an arm's width is covered by segments below a loss threshold
.arm_loss_fraction <- function(segments, region, threshold = 0.85) {
    seg <- segments[segments$chrom == region$chrom &
                     segments$loc.start < region$end &
                     segments$loc.end   > region$start, ]
    if (nrow(seg) == 0) return(NA_real_)
    w_start <- pmax(seg$loc.start, region$start)
    w_end   <- pmin(seg$loc.end, region$end)
    width   <- w_end - w_start
    is_loss <- seg$seg.mean < threshold
    sum(width[is_loss]) / sum(width)
}

# hg38 approximate centromere-based arm boundaries
ARM_1P  <- list(chrom = "chr1",  start = 0,        end = 123400000)
ARM_1Q  <- list(chrom = "chr1",  start = 123400000, end = 248956422)
ARM_19P <- list(chrom = "chr19", start = 0,        end = 26200000)
ARM_19Q <- list(chrom = "chr19", start = 26200000, end = 58617616)


purity_from_hemizygous_loss <- function(copy_ratio) {
    2 * (1 - copy_ratio)
}
detect_1p19q_codel <- function(segments,
                                min_loss_depth = 0.114,          # ~2x empirical noise SD -- minimum depth to count as real loss
                                retained_depth_ceiling = 0.114,  # below this, partner arm is indistinguishable from noise
                                magnitude_ratio_threshold = 0.5) {

    seg_s <- segments
    cr_1p  <- .arm_copy_ratio(seg_s, ARM_1P)
    cr_1q  <- .arm_copy_ratio(seg_s, ARM_1Q)
    cr_19p <- .arm_copy_ratio(seg_s, ARM_19P)
    cr_19q <- .arm_copy_ratio(seg_s, ARM_19Q)

    if (any(is.na(c(cr_1p, cr_1q, cr_19p, cr_19q)))) {
        return(list(codel = NA, cr_1p = cr_1p, cr_1q = cr_1q, cr_19p = cr_19p, cr_19q = cr_19q,
                    reason = "insufficient segmentation coverage"))
    }

    depth_of_loss <- function(cr) pmax(0, 1 - cr)
    depth_1p  <- depth_of_loss(cr_1p);  depth_1q  <- depth_of_loss(cr_1q)
    depth_19p <- depth_of_loss(cr_19p); depth_19q <- depth_of_loss(cr_19q)

    # An arm is "lost" once its loss depth clears the noise floor -- no longer
    # gated on any single segment crossing a fixed absolute copy-ratio cutoff
    p_arm_lost <- depth_1p  >= min_loss_depth
    q_arm_lost <- depth_19q >= min_loss_depth

    # Partner arm is "retained" if its own loss is within noise, OR shallow
    # relative to the codeleted arm (distinct/secondary event, as before)
    q1_retained  <- (depth_1q  < retained_depth_ceiling) || (depth_1q  < magnitude_ratio_threshold * depth_1p)
    p19_retained <- (depth_19p < retained_depth_ceiling) || (depth_19p < magnitude_ratio_threshold * depth_19q)

    codel <- p_arm_lost && q_arm_lost && q1_retained && p19_retained

    reason <- if (codel) {
        sprintf("true 1p/19q codeletion (depth 1p=%.3f, 19q=%.3f vs partner depths 1q=%.3f, 19p=%.3f)",
                depth_1p, depth_19q, depth_1q, depth_19p)
    } else if (!p_arm_lost || !q_arm_lost) {
        "1p/19q loss criteria not met (loss depth below noise floor)"
    } else {
        "1p and 19q both lost, but 1q/19p also lost at comparable depth -- likely whole-chr1/chr19 loss, not true codel"
    }

    list(codel = codel, cr_1p = cr_1p, cr_1q = cr_1q, cr_19p = cr_19p, cr_19q = cr_19q,
         depth_1p = depth_1p, depth_1q = depth_1q, depth_19p = depth_19p, depth_19q = depth_19q,
         reason = reason)
}


purity_from_1p19q <- function(codel_result, segments = NULL,
                               discrepancy_threshold = 0.15) {
    if (is.null(codel_result) || is.na(codel_result$codel) || !codel_result$codel) {
        return(list(purity = NA_real_, flag = NA_character_))
    }

    p_1p  <- purity_from_hemizygous_loss(codel_result$cr_1p)
    p_19q <- purity_from_hemizygous_loss(codel_result$cr_19q)
    disagreement <- abs(p_1p - p_19q)

    if (disagreement > discrepancy_threshold) {
        return(list(
            purity = p_1p,
            flag = sprintf(
                "1p (%.3f) and 19q (%.3f) purity estimates disagree by %.3f; using 1p (larger, more robust arm). Possible additional subclonal loss on 19q -- worth manual review.",
                p_1p, p_19q, disagreement)
        ))
    }

    if (is.null(segments)) {
        stop("purity_from_1p19q: segments is required when arms agree (needed to pool copy ratios). Did you forget to pass it?")
    }

    seg_1p  <- segments[segments$chrom == ARM_1P$chrom  & segments$loc.start < ARM_1P$end  & segments$loc.end > ARM_1P$start, ]
    seg_19q <- segments[segments$chrom == ARM_19Q$chrom & segments$loc.start < ARM_19Q$end & segments$loc.end > ARM_19Q$start, ]
    w_1p  <- pmin(seg_1p$loc.end, ARM_1P$end)   - pmax(seg_1p$loc.start, ARM_1P$start)
    w_19q <- pmin(seg_19q$loc.end, ARM_19Q$end) - pmax(seg_19q$loc.start, ARM_19Q$start)
    cr_pooled <- (sum(seg_1p$seg.mean * w_1p) + sum(seg_19q$seg.mean * w_19q)) / (sum(w_1p) + sum(w_19q))
    
    list(purity = purity_from_hemizygous_loss(cr_pooled), flag = NA_character_)
}


codel_result <- detect_1p19q_codel(segments)
purity_1p19q <- purity_from_1p19q(codel_result, segments = segments)

#-------------------------------------------------------------------------------
# 1.3  Select ACE model using purity anchors
#-------------------------------------------------------------------------------
select_purity_from_ACE <- function(ACE_results, anchor_purities,
                                   ploidy_target = 2,
                                   min_ploidy = 1.5,
                                   margin = 0.05) {

    anchor_purities <- anchor_purities[!is.na(anchor_purities)]
    candidates <- ACE_results %>% filter(minimum == TRUE) %>% filter(ploidy >= min_ploidy)

    if (length(anchor_purities) == 0) {
        flog.warn("No anchor purity available; falling back to lowest-error near-diploid ACE solution.")
        chosen <- candidates %>%
            mutate(ploidy_dist = abs(ploidy - ploidy_target)) %>%
            arrange(ploidy_dist, error) %>%
            dplyr::slice(1)
        return(list(chosen = chosen, anchor_consensus = NA_real_,
                    anchor_bounds = c(NA_real_, NA_real_), candidates_in_bounds = NA))
    }

    anchor_consensus <- median(anchor_purities)
    lower <- min(anchor_purities) - margin
    upper <- max(anchor_purities) + margin
    
    in_bounds <- candidates %>% filter(cellularity >= lower, cellularity <= upper)

    if (nrow(in_bounds) == 0) {
        flog.warn("No ACE local minima fall within anchor-derived bounds [%.3f, %.3f]; using closest candidate instead.",
                  lower, upper)
        chosen <- candidates %>%
            mutate(dist_to_consensus = abs(cellularity - anchor_consensus),
                   ploidy_dist = abs(ploidy - ploidy_target)) %>%
            arrange(dist_to_consensus, ploidy_dist, error) %>%
            dplyr::slice(1)
        return(list(chosen = chosen, anchor_consensus = anchor_consensus,
                    anchor_bounds = c(lower, upper), candidates_in_bounds = in_bounds))
    }

    chosen <- in_bounds %>%
        mutate(ploidy_dist = abs(ploidy - ploidy_target)) %>%
        arrange(ploidy_dist, error) %>%
        dplyr::slice(1)

    list(chosen = chosen, anchor_consensus = anchor_consensus,
         anchor_bounds = c(lower, upper), candidates_in_bounds = in_bounds)
}

use_TP53_anchor <- is.na(codel_result$codel) || !isTRUE(codel_result$codel)

TP53_purity_final <- if (use_TP53_anchor) TP53_purity else NA_real_

if (!use_TP53_anchor && !is.na(TP53_purity)) {
    flog.info("1p/19q codeletion detected; excluding TP53 (%.3f) from anchor consensus.", TP53_purity)
}
resolve_final_purity_ploidy <- function(ace_selection, anchor_consensus,
                                         discordance_threshold = 0.19,
                                         ploidy_target = 2) {

    ace_purity <- ace_selection$chosen$cellularity
    ace_ploidy <- ace_selection$chosen$ploidy

    if (length(ace_purity) == 0 || length(ace_ploidy) == 0) {
        flog.warn("No ACE candidate solution available (0 rows after filtering); falling back to anchor consensus with target ploidy %.1f.",
                   ploidy_target)
        if (is.na(anchor_consensus)) {
            stop("resolve_final_purity_ploidy: no ACE candidates AND no anchor consensus -- cannot determine a purity estimate for this sample.")
        }
        return(list(final_purity = anchor_consensus, final_ploidy = ploidy_target,
                    purity_source = "anchor_consensus_no_ace_candidates", discordance = NA_real_))
    }

    if (is.na(anchor_consensus)) {
        return(list(final_purity = ace_purity, final_ploidy = ace_ploidy,
                    purity_source = "ace_selection", discordance = NA_real_))
    }

    discordance <- abs(ace_purity - anchor_consensus)

    if (discordance > discordance_threshold) {
        flog.warn(
            "ACE purity (%.2f) and anchor consensus (%.2f) disagree by %.3f (> %.2f); using anchor consensus, ploidy defaulted to %.1f.",
            ace_purity, anchor_consensus, discordance, discordance_threshold, ploidy_target
        )
        return(list(final_purity = anchor_consensus, final_ploidy = ploidy_target,
                    purity_source = "anchor_consensus", discordance = discordance))
    }

    list(final_purity = ace_purity, final_ploidy = ace_ploidy,
         purity_source = "ace_selection", discordance = discordance)
}

anchor_purities <- c(IDH_purity, TP53_purity_final, purity_1p19q$purity)


ace_selection <- select_purity_from_ACE(ACE_results, anchor_purities,
                                        ploidy_target = 2, margin = opt$margin)

resolved <- resolve_final_purity_ploidy(ace_selection, ace_selection$anchor_consensus,
                                         discordance_threshold = 0.15, ploidy_target = 2)
#-------------------------------------------------------------------------------
# 2.1  Prepare output files
#-------------------------------------------------------------------------------
final_purity <- resolved$final_purity
final_ploidy <- resolved$final_ploidy

anchor_names <- c("IDH", "TP53", "1p19q")
anchor_values <- c(
    if (length(IDH_purity) > 0) IDH_purity[1] else NA_real_,
    if (length(TP53_purity) > 0) TP53_purity[1] else NA_real_,
    purity_1p19q$purity
)
names(anchor_values) <- anchor_names

anchors_used_str  <- paste(anchor_names[!is.na(anchor_values)], collapse = ";")
anchors_value_str <- paste(sprintf("%s=%.4f", anchor_names[!is.na(anchor_values)],
                                    anchor_values[!is.na(anchor_values)]), collapse = ";")

pureCN_defaults <- list(
    min_purity = formals(PureCN::runAbsoluteCN)$test.purity[[2]],
    max_purity = formals(PureCN::runAbsoluteCN)$test.purity[[3]],
    min_ploidy = formals(PureCN::runAbsoluteCN)$min.ploidy,
    max_ploidy = formals(PureCN::runAbsoluteCN)$max.ploidy
)
build_purity_ploidy_bounds <- function(resolved, n_anchors,
                                        margin = 0.05,
                                        ploidy_margin = 0.3,
                                        low_confidence_ploidy_target = 2,
                                        low_confidence_threshold = 2,
                                        extreme_ploidy_lower = 1.5,
                                        extreme_ploidy_upper = 3.5,
                                        defaults = pureCN_defaults) {
    if (n_anchors == 0) {
        flog.warn("No anchors available; using PureCN's own defaults (purity %.2f-%.2f, ploidy %.2f-%.2f).",
                   defaults$min_purity, defaults$max_purity, defaults$min_ploidy, defaults$max_ploidy)
        return(list(
            min_purity = defaults$min_purity, max_purity = defaults$max_purity,
            min_ploidy = defaults$min_ploidy, max_ploidy = defaults$max_ploidy,
            confidence_flag = "no_anchors_pureCN_default"
        ))
    }

    final_purity <- resolved$final_purity
    final_ploidy <- resolved$final_ploidy

    low_confidence <- n_anchors < low_confidence_threshold
    extreme_ploidy <- final_ploidy < extreme_ploidy_lower || final_ploidy > extreme_ploidy_upper

    if (isTRUE((low_confidence || extreme_ploidy) &&
        abs(final_ploidy - low_confidence_ploidy_target) > ploidy_margin)) {

        ploidy_lo <- min(final_ploidy, low_confidence_ploidy_target) - ploidy_margin
        ploidy_hi <- max(final_ploidy, low_confidence_ploidy_target) + ploidy_margin
        confidence_flag <- if (extreme_ploidy) "extreme_ploidy_widened" else "low_confidence_ploidy_widened"
    } else {
        ploidy_lo <- final_ploidy - ploidy_margin
        ploidy_hi <- final_ploidy + ploidy_margin
        confidence_flag <- ifelse(low_confidence, "low_confidence_single_anchor", "high_confidence")
    }

    if (resolved$purity_source %in% c("anchor_consensus", "anchor_consensus_no_ace_candidates")) {
        confidence_flag <- "anchor_override"
    }

    list(
        min_purity = round(max(0.05, final_purity - margin), 4),
        max_purity = round(min(0.99, final_purity + margin), 4),
        min_ploidy = round(max(1.0, ploidy_lo), 4),
        max_ploidy = round(ploidy_hi, 4),
        confidence_flag = confidence_flag
    )
}

bounds <- build_purity_ploidy_bounds(resolved, n_anchors = sum(!is.na(anchor_values)))

bounds

result <- data.frame(
    n_anchors = sum(!is.na(anchor_values)),
    anchors_used = anchors_used_str,
    anchors_detail = anchors_value_str,
    IDH_purity = ifelse(is.na(anchor_values["IDH"]), NA_real_, round(anchor_values["IDH"], 4)),
    TP53_purity = ifelse(is.na(anchor_values["TP53"]), NA_real_, round(anchor_values["TP53"], 4)),
    purity_1p19q = ifelse(is.na(anchor_values["1p19q"]), NA_real_, round(anchor_values["1p19q"], 4)),
    codel_1p19q = codel_result$codel,
    codel_reason = codel_result$reason,
    anchor_consensus = round(ace_selection$anchor_consensus, 4),
    anchor_bounds_lower = round(ace_selection$anchor_bounds[1], 4),
    anchor_bounds_upper = round(ace_selection$anchor_bounds[2], 4),
    ace_candidates_in_bounds = ifelse(is.data.frame(ace_selection$candidates_in_bounds),
                                       nrow(ace_selection$candidates_in_bounds), NA_integer_),
    ace_ploidy = if (length(ace_selection$chosen$ploidy) == 0) NA_real_ else ace_selection$chosen$ploidy,
    ace_purity = if (length(ace_selection$chosen$cellularity) == 0) NA_real_ else ace_selection$chosen$cellularity,
    purity_discordance = round(resolved$discordance, 4),
    purity_source = resolved$purity_source,
    final_purity = round(final_purity, 4),
    final_ploidy = round(final_ploidy, 4),
    min_purity = bounds$min_purity,
    max_purity = bounds$max_purity,
    min_ploidy = bounds$min_ploidy,
    max_ploidy = bounds$max_ploidy,
    confidence = bounds$confidence_flag
)

write.csv(result, file = opt$out, row.names = FALSE)

# Also write a simple, shell-sourceable bounds file (no quoting/CSV parsing needed downstream)
bounds_env_file <- sub("\\.csv$", ".env", opt$out)
writeLines(c(
    sprintf("min_purity=%s", result$min_purity),
    sprintf("max_purity=%s", result$max_purity),
    sprintf("min_ploidy=%s", result$min_ploidy),
    sprintf("max_ploidy=%s", result$max_ploidy)
), con = bounds_env_file)
