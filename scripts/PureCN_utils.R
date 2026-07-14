
Fetch_IDH_mut <- function(ann){
    hotspot_idx <- c()
    for (i in seq_along(ann)) {
        entries <- unlist(strsplit(as.character(ann[[i]]), ","))
        for (x in entries) {
            fields <- unlist(strsplit(x, "\\|"))
            if (length(fields) < 11) next
            gene <- fields[4]
            protein <- fields[11]
            # remove p.
            protein <- gsub("^p\\.", "", protein)

            idh1_hits <- c("Arg132His","Arg132Cys","Arg132Gly",
                           "Arg132Ser","Arg132Leu")
            idh2_hits <- c("Arg140Gln",
                           "Arg172Lys","Arg172Met","Arg172Trp",
                           "Arg172Gly","Arg172Ser")
            if ((gene == "IDH1" && protein %in% idh1_hits) ||
                (gene == "IDH2" && protein %in% idh2_hits)) {
                hotspot_idx <- c(hotspot_idx, i)
                break
            }
        }
    }

    hotspot_idx <- unique(hotspot_idx)
    return(hotspot_idx)
    
}


SetPriorVcf_IDH_mutant <- function(vcf,
                                   prior.somatic = c(0.5, 5e-04, 0.999, 1e-04, 0.995, 5e-04),
                                   ...) {
    # First let PureCN assign default priors
    vcf <- PureCN::setPriorVcf(vcf, prior.somatic = prior.somatic , ...)

    # Fetch annotations
    ann <- info(vcf)$ANN
    if (is.null(ann)) return(vcf)
    # Fetch IDH mutation
    hotspot_idx <- Fetch_IDH_mut(ann)

    if (length(hotspot_idx) > 0) {
        info(vcf)$PureCN.PR[hotspot_idx] <- 0.9999
        flog.info("Set high somatic prior for ", length(hotspot_idx),
                " IDH hotspot mutation(s).")
    }

    vcf
}

SetPriorPurity <- function(vcf_file, test.purity, stringency = 5,
                            min.likelihood.floor = 1e-4) {

    vcf <- VariantAnnotation::readVcf(vcf_file)
    ann <- info(vcf)$ANN
    hotspot_idx <- Fetch_IDH_mut(ann)

    if (length(hotspot_idx) == 0) {
        return(rep(1, length(test.purity)) / length(test.purity))
    }
    if (length(hotspot_idx) > 1) hotspot_idx <- hotspot_idx[1]

    ad <- geno(vcf)$AD[hotspot_idx, ][[1]]
    alt_reads   <- ad[2]
    total_reads <- ad[1] + ad[2]

    expected_vaf <- test.purity / 2
    expected_vaf <- pmin(pmax(expected_vaf, 1e-4), 1 - 1e-4)

    ll <- dbinom(alt_reads, total_reads, expected_vaf, log = TRUE)
    ll <- ll * stringency
    ll <- ll - max(ll)
    likelihood <- exp(ll)
    likelihood <- pmax(likelihood, min.likelihood.floor * max(likelihood))
    likelihood / sum(likelihood)
}

recover_filtered_variants <- function(vcf_raw, gene_name, ret, seg_file,
                                      min_af, min_dp, protein_patterns = NULL) {
    # Use hotspot-aware matching if patterns provided, else gene-level
    if (!is.null(protein_patterns)) {
        gene_idx <- which(annotates_to_hotspot(info(vcf_raw)$ANN, gene_name, 
                                               protein_patterns))
    } else {
        gene_idx <- which(annotates_to_gene(info(vcf_raw)$ANN, gene_name))
    }
    if (length(gene_idx) == 0) return(NULL)
    gene_vcf <- vcf_raw[gene_idx, ]

    # Filter PASS, min AF, min depth
    gene_vcf <- gene_vcf[rowRanges(gene_vcf)$FILTER == "PASS", ]
    if (nrow(gene_vcf) == 0) return(NULL)

    # Only apply AF/DP filters for gene-level recovery, not hotspots
    if (is.null(protein_patterns)) {
        keep <- as.numeric(geno(gene_vcf)$AF) >= min_af &
                as.numeric(geno(gene_vcf)$DP) >= min_dp
        gene_vcf <- gene_vcf[keep, ]
        if (nrow(gene_vcf) == 0) return(NULL)
    }

    AF           <- as.numeric(geno(gene_vcf)$AF)
    DP           <- as.numeric(geno(gene_vcf)$DP)
    POPAF        <- as.numeric(info(gene_vcf)$POPAF)
    popAF_linear <- 10^(-POPAF)
    var_chr      <- as.character(seqnames(gene_vcf))
    var_pos      <- start(gene_vcf)

    if (all(var_chr == "chrX")) {
        segments  <- read.delim(seg_file)
        seg_mean  <- mean(segments[segments$chrom == "chrX", "seg.mean"])
        log_ratio <- rep(log2(seg_mean), length(var_pos))
    } else {
        purecn_seg <- ret$results[[1]]$seg
        log_ratio  <- sapply(seq_along(var_pos), function(i) {
            seg_row <- purecn_seg[purecn_seg$chrom == var_chr[i] &
                                  purecn_seg$loc.start <= var_pos[i] &
                                  purecn_seg$loc.end   >= var_pos[i], ]
            if (nrow(seg_row) == 0) return(NA)
            seg_row$log.ratio[1]
        })
    }

    data.frame(
        chr               = var_chr,
        start             = var_pos,
        end               = end(gene_vcf),
        ID                = names(gene_vcf),
        REF               = as.character(ref(gene_vcf)),
        ALT               = sapply(alt(gene_vcf), function(x)
            paste(as.character(x), collapse = ",")),
        ML.SOMATIC        = popAF_linear < 0.0001,
        POSTERIOR.SOMATIC = ifelse(popAF_linear < 0.0001, 1 - popAF_linear,
                                   popAF_linear),
        AR                = AF,
        AR.ADJUSTED       = AF,
        log.ratio         = log_ratio,
        depth             = DP,
        prior.somatic     = 1 - popAF_linear,
        on.target         = 1L,
        pon.count         = as.integer(info(gene_vcf)$HMF_PON_SC),
        gene.symbol       = gene_name,
        recovered         = TRUE,
        stringsAsFactors  = FALSE
    )
}

.getFilePrefix <- function(out, sampleid) {
    isDir <- file.info(out)$isdir
    if (!is.na(isDir) && isDir) return(file.path(out, sampleid))
    out
}



.checkFileList <- function(file) {
    files <- read.delim(file, as.is = TRUE, header = FALSE)[, 1]
    numExists <- sum(file.exists(files), na.rm = TRUE)
    if (numExists < length(files)) {
        stop("File not exists in file ", file)
    }
    files
}


.getNormalCoverage <- function(normal.coverage.file) {
    if (!is.null(normalDB)) {
        if (is.null(normal.coverage.file)) {
            normal.coverage.file <- calculateTangentNormal(tumor.coverage.file,
                                                           normalDB)
        }
    } else if (is.null(normal.coverage.file) && is.null(seg.file) &&
               is.null(log.ratio)) {
        stop("Need either normalDB or normal.coverage.file")
    }
    normal.coverage.file
}

annotates_to_gene <- function(ann_list, gene_name) {
  sapply(ann_list, function(anns) {
    any(sapply(anns, function(ann) {
      fields <- strsplit(ann, "\\|")[[1]]
      length(fields) >= 4 && fields[4] == gene_name
    }))
  })
}

# Check if a variant annotates to a specific gene AND protein change
annotates_to_hotspot <- function(ann_list, gene_name, protein_patterns) {
  sapply(ann_list, function(anns) {
    any(sapply(anns, function(ann) {
      fields <- strsplit(ann, "\\|")[[1]]
      if (length(fields) < 11) return(FALSE)
      gene    <- fields[4]
      protein <- fields[11]
      gene == gene_name && any(sapply(protein_patterns, function(p) grepl(p, protein)))
    }))
  })
}
