############################################################
# Heuristic screen for potential chromothripsis cases
#
# Purpose:
# This script performs an initial heuristic search for samples/chromosomes
# with copy-number patterns suggestive of chromothripsis.
#
# The approach is based on FACETS-derived copy-number segments and looks for:
#   1. Chromosomes with multiple copy-number segments
#   2. Multiple total copy-number states
#   3. Oscillation between a baseline copy-number state and higher states
#   4. Sufficient segment count to prioritize cases for manual review
#
# Important:
# This is not intended to be a definitive chromothripsis caller.
# Candidate cases identified by this heuristic should be manually reviewed
# using copy-number plots, allele-specific copy number, pathology context,
# and other available genomic/clinical information.
############################################################


############################################################
# 1. Load FACETS copy-number segment file
############################################################

# Read FACETS annotated copy-number segment file.
# Expected columns include:
#   ID
#   chrom
#   loc.start
#   loc.end
#   tcn
#   lcn
#
# Here:
#   tcn = total copy number
#   lcn = lesser/minor copy number
cncf <- read.csv(
  "path/to/impact_facets_annotated.cncf.txt.gz",
  header = TRUE,
  as.is = TRUE,
  sep = "\t"
)

# Extract the DMP/sample identifier from the FACETS ID.
cncf$sample <- str_sub(cncf$ID, 1, 17)

# Restrict copy-number data to samples included in the study/sample list.
cncf <- cncf[cncf$sample %in% nsp$SAMPLE_ID,]


############################################################
# 2. Initialize chromosome-level summary table
############################################################

# Create an empty data frame to store per-sample, per-chromosome metrics.
#
# Columns:
#   dmp_sample          = sample ID
#   chrom               = chromosome
#   n_seg               = number of merged copy-number segments
#   n_tcn_state         = number of unique total copy-number states
#   n_oscil             = number of oscillating baseline-copy-number blocks
#   basecn              = inferred baseline copy-number state
#   n_seg_basecn        = number of segments at or below baseline copy number
#   length_seg_basecn   = total genomic length of baseline/low-copy segments
#   length_seg_total    = total genomic length covered by all segments
csum <- as.data.frame(matrix(NA, nrow = 0, ncol = 9))

colnames(csum) <- c(
  "dmp_sample",
  "chrom",
  "n_seg",
  "n_tcn_state",
  "n_oscil",
  "basecn",
  "n_seg_basecn",
  "length_seg_basecn",
  "length_seg_total"
)

# Explicitly set column classes.
# This helps avoid unintended factor/string behavior in older R workflows.
csum$dmp_sample <- as.character(csum$dmp_sample)
csum$chrom <- as.character(csum$chrom)
csum$n_seg <- as.numeric(as.character(csum$n_seg))
csum$n_oscil <- as.numeric(as.character(csum$n_oscil))
csum$basecn <- as.numeric(as.character(csum$basecn))
csum$n_seg_basecn <- as.numeric(as.character(csum$n_seg_basecn))
csum$length_seg_basecn <- as.numeric(as.character(csum$length_seg_basecn))
csum$length_seg_total <- as.numeric(as.character(csum$length_seg_total))


############################################################
# 3. Loop over pulmonary NET samples
############################################################

# net is assumed to be the main pulmonary NET annotation table.
# Required columns include:
#   dmp_sample
#   study_id
#   facets_id
#   dmp_pt
#   decision
#   facets_qc
for (i in 1:nrow(net)) {
  
  print(paste0("Processing ", net$dmp_sample[i], ": ", i))
  
  
  ############################################################
  # 3a. Retrieve FACETS copy-number segments for this sample
  ############################################################
  
  # If the sample is not present in the combined cncf object,
  # try to locate a sample-specific FACETS output file.
  if (!(net$dmp_sample[i] %in% unique(cncf$sample))) {
    
    # Search for FACETS files using the patient's study ID.
    # Only high-sensitivity FACETS output files are retained.
    infoline <- list.files(
      path = paste0("analysis/facets_files/Custom_Tischfield"),
      pattern = net$study_id[i],
      full.names = TRUE
    )
    
    infoline <- infoline[grepl("_hisens", infoline)]
    
    # If only one matching file exists, use it.
    if (length(infoline) == 1) {
      
      df <- read.csv(
        infoline,
        header = TRUE,
        as.is = TRUE,
        sep = "\t"
      )
      
    } else {
      
      # If multiple matching files exist, prioritize refit files when present.
      if (sum(grepl("refit", infoline)) != 0) {
        infoline <- infoline[grepl("refit", infoline)]
      }
      
      # If there are still multiple files, use the last matching file.
      df <- read.csv(
        infoline[length(infoline)],
        header = TRUE,
        as.is = TRUE,
        sep = "\t"
      )
    }
    
    # Assign the current DMP sample ID to the loaded segment table.
    df$sample <- net$dmp_sample[i]
    
  } else if (nrow(cncf[cncf$sample == net$dmp_sample[i],]) != 0) {
    
    # If the sample exists in the combined cncf object,
    # retrieve the corresponding FACETS high-sensitivity segment profile.
    #
    # facets$tumor_sample_id and facets$id are used to map from net$facets_id
    # to the FACETS sample ID.
    df <- cncf[
      cncf$ID == paste0(
        facets$tumor_sample_id[facets$id == net$facets_id[i]],
        "_hisens"
      ),
    ]
    
  } else {
    
    # If no copy-number information can be found, skip this sample.
    print(paste0("Sample ", net$dmp_sample[i], " has no copy number information"))
    next
  }
  
  
  ############################################################
  # 3b. Identify chromosomes with enough segmentation to evaluate
  ############################################################
  
  # Consider only chromosomes with at least 5 FACETS segments.
  # This threshold removes chromosomes with too few segments to evaluate
  # for oscillatory copy-number patterns.
  if (length(table(df$chrom)[table(df$chrom) >= 5]) != 0) {
    
    # Loop through chromosomes meeting the minimum segment-count criterion.
    for (j in names(table(df$chrom)[table(df$chrom) >= 5])) {
      
      # Add one row to the summary table for this sample/chromosome.
      csum[nrow(csum) + 1,] <- NA
      csum$dmp_sample[nrow(csum)] <- net$dmp_sample[i]
      csum$chrom[nrow(csum)] <- j
      
      
      ############################################################
      # 3c. Extract and simplify copy-number segments for chromosome j
      ############################################################
      
      # Keep only segments from the current chromosome.
      tempdf <- df[df$chrom == j,]
      
      # Keep minimal segment-level fields required for this heuristic.
      tempdf <- tempdf[, c(
        "sample",
        "chrom",
        "loc.start",
        "loc.end",
        "tcn",
        "lcn"
      )]
      
      # Identify adjacent segments with identical total and minor copy number.
      tempdf$contiguous <- ""
      
      for (k in 2:nrow(tempdf)) {
        
        if (
          tempdf$tcn[k] == tempdf$tcn[k - 1] &
          !is.na(tempdf$lcn[k]) &
          !is.na(tempdf$lcn[k - 1]) &
          tempdf$lcn[k] == tempdf$lcn[k - 1]
        ) {
          tempdf$contiguous[k] <- "yes"
        }
      }
      
      # Get starting positions of segments that should be merged into
      # the immediately preceding segment.
      infoline <- tempdf$loc.start[tempdf$contiguous == "yes"]
      
      # Merge adjacent segments with identical TCN and LCN.
      # For each contiguous segment:
      #   - extend the previous segment's end coordinate
      #   - remove the current segment
      for (k in infoline) {
        
        tempdf$loc.end[which(tempdf$loc.start == k) - 1] <-
          tempdf$loc.end[which(tempdf$loc.start == k)]
        
        tempdf <- tempdf[tempdf$loc.start != k,]
      }
      
      
      ############################################################
      # 3d. Calculate basic chromosome-level copy-number metrics
      ############################################################
      
      # Number of merged copy-number segments on this chromosome.
      csum$n_seg[nrow(csum)] <- nrow(tempdf)
      
      # Number of distinct total copy-number states.
      csum$n_tcn_state[nrow(csum)] <- length(unique(tempdf$tcn))
      
      
      ############################################################
      # 3e. Estimate baseline copy-number state and oscillation count
      ############################################################
      
      # Continue only if at least 3 merged segments remain.
      # Fewer than 3 segments are insufficient for this oscillation heuristic.
      if (nrow(tempdf) >= 3) {
        
        # Calculate segment length.
        tempdf$len <- tempdf$loc.end - tempdf$loc.start
        
        # Define a baseline copy-number state.
        #
        # Here, basecn is assigned as the lowest total copy-number state
        # observed more than once.
        csum$basecn[nrow(csum)] <-
          as.numeric(names(table(tempdf$tcn)[table(tempdf$tcn) > 1][1]))
        
        # Label segments at or below the inferred baseline copy number.
        # These are treated as the "baseline/low" state for oscillation counting.
        tempdf$contiguous[tempdf$tcn <= csum$basecn[nrow(csum)]] <- "base"
        
        # Count how many separate baseline/low-copy blocks exist.
        #
        # The idea is:
        #   chromothripsis-like profiles often show repeated oscillation
        #   between lower and higher copy-number states across a chromosome.
        #
        # v starts at 1 and increments each time a new baseline/low-copy block
        # begins after a non-baseline segment.
        v <- 1
        
        for (k in 2:nrow(tempdf)) {
          
          if (
            tempdf$contiguous[k] == "base" &
            tempdf$contiguous[k - 1] != "base"
          ) {
            v <- v + 1
          }
        }
        
        # Store the estimated number of oscillating baseline-copy blocks.
        csum$n_oscil[nrow(csum)] <- v
        
        # Store number of segments at or below baseline copy number.
        csum$n_seg_basecn[nrow(csum)] <-
          nrow(tempdf[tempdf$tcn <= csum$basecn[nrow(csum)],])
        
        # Store total genomic length of segments at or below baseline copy number.
        csum$length_seg_basecn[nrow(csum)] <-
          sum(tempdf$len[tempdf$tcn <= csum$basecn[nrow(csum)]])
        
        # Store total genomic length covered by all merged segments.
        csum$length_seg_total[nrow(csum)] <- sum(tempdf$len)
      }
    }
  }
}


############################################################
# 10. Export candidate tables for manual review
############################################################

# Export per-chromosome summary table for samples/chromosomes meeting
# a relaxed review threshold:
#   n_oscil >= 3
#   n_seg >= 5
#
# This table is intended for manual review of candidate chromosomes.
write.table(
  csum[
    csum$n_oscil >= 3 &
      csum$n_seg >= 5,
  ],
  "path/to/per_chromosome_summary.tsv",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = "\t"
)

# Export per-case annotation table for samples that have at least one
# chromosome meeting the relaxed review threshold.
#
# This is intended to support case-level manual review.
write.table(
  net[
    net$dmp_sample %in%
      csum$dmp_sample[
        csum$n_oscil >= 3 &
          csum$n_seg >= 5
      ],
  ],
  "path/to/per_case_annotation.tsv",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = "\t"
)