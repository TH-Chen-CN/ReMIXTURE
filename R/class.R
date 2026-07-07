#' ReMIXTURE
#'
#' @description
#' ReMIXTURE is a tool for visualising how genetic (or any) diversity in a population is distributed and shared between groups. Its canonical purpose is to describe the genetic diversity among individuals allocated to different geographic regions.
#'
#' It produces intuitive (and [pretty](https://github.com/mtrw/ReMIXTURE/blob/main/README.md)!) plots.
#' These plots are based on ReMIXTURE's diversity metric, which is designed to be maximally intuitive. They work like Venn diagrams: Diversity in a group is partitioned into diversity that is unique to that group, and that which is overlapped by the diversity in other groups.
#'
#' ReMIXTURE will work on any dataset given two pieces of data:
#' - (1) A matrix of pairwise distances between a collection of samples (say, at least 10?) from each of some collection of regions (or groups, more broadly). This may contain any metric you choose such as IBS distances for genetic data. ReMIXTURE can automatically create an IBS distance matrix given `.vcf` / `.vcf.gz` file, using the package `SNPRelate`, this will require an additional table mapping the sample names in the VCF file to regions.
#' - (2) A table describing where on Earth each region/group is located (latitude and longitude), which is used for plotting only.
#'
#' ## *GETTING STARTED?*
#' Scroll down and follow the short tutorial in the *examples* section.
#'
#' @examples
#'
#' ####################################
#' #### A SHORT TUTORIAL ####
#' ####################################
#'
#' # Load the library
#' library(ReMIXTURE)
#'
#' # ReMIXTURE is written using the [R6](https://r6.r-lib.org/index.html) Object Orientated Programming system. You can read more about OOP from the great Hadley Wickham [here](https://adv-r.hadley.nz/oo.html) or about R6 specifically [here](https://adv-r.hadley.nz/r6.html).
#'
#' # Create ReMIXTURE analysis object from an included dataset (Tripodi & Rabanus-Wallace et al. 2021. PNAS).
#' rm <- ReMIXTURE$new(
#'   distance_matrix = ReMIXTURE_example_distance_matrix,
#'   region_table = ReMIXTURE_example_region_table
#' )
#'
#' # An MDS plot is a quick way to get a feel for the structure of the data.
#' rm$plot_MDS()
#' #
#' # A run involves taking a subsample of all the samples in the distance matrix, such that every region has the same number of samples included. This subsample is hierarchically clustered (imagine creating a dendrogram and cutting it at a certain height). Clusters are then counted. The number of clusters in which a region appears is a proxy for that region's total genetic diversity. The number of clusters in which *only* one region appears is a proxy for the diversity unique to that region. The number of clusters in which members of a pair of clusters appears is a proxy for the diversity that is overlapping between those two regions. the concept is similar to the way a Venn diagram might work. This is done `iterations` times, and all these various counts are averaged at the end.
#' #
#' # How is it decided how high the tree is cut? The details are in `?ReMIXTURE` under '$run'.
#' # Ideally, a collection of parameter settings are chosen and the results assessed to choose a good one. By default, the run function tries to choose a good selection of values, and it does a good job.
#' rm$run(iterations = 25) # Iterations set low as we will be doing several runs and comparing results, so we can sacrifice a little accuracy for now.
#'
#' # For in an intuitive look at how the diversity metric changes with the parameter settings, view them simply as a grid, for each run (press <enter> to get through them). Notice how with low H values, unique diversity is over-emphasised and overlapping diversity is under-emphasised
#' rm$plot_results_grid()
#'
#' # - Circle sizes represent the diversity within a region (diagonals) or overlapped between regions (off-diagonals)
#' # - The white inner circle represents the diversity unique to the focal region. The black outer band therefore represents the amount of diversity in the focal region that is overlapped by diversity found in other regions.
#'
#' # Comparing the results between runs should make it clear why a ReMIXTURE analysis on one data set can not be straightforwardly compared. Certainly, any parameter settings will reliably give an indication of which regions are more or less diverse, which have more unique diversity, and which overlap more. But by altering the distribution of H, we can emphasise different features more or less.
#' # To create completely comparable analyses, therefore, you will need to ensure both datasets
#'
#' # A good balance can be found by summing region-unique/-overlapping clusters. At some point they are expected to cross over.
#' rm$plot_h_optimisation()
#'
#' # Note also that the number of multi-region clusters reaches a maximum, when the H is big enough to produce many such clusters, but before the clustersbecome so large that the total cluster number decline enough to drag this number down.
#' # H values between the crossover point and the multi-region maximum are usually good for showing a balanced representation. I prefer the crossover point, here about run 4.
#' # Some more diagnosis can be done, but for now let's produce our remixture maps!
#' # BEWARE: What we are doing is great for a demo but before you publish you should be sure to do the rest of this tutorial.
#' rm$plot_maps(run = 4, focalRegion = "Africa")
#'
#' # Note the overlapped diversity is now shown by the width/alpha of a line joining the focal region to each other region.
#'
#' # Let's make those diversity visualisations bigger on the plot and get some curved lines ...
#' set.seed(42.42)
#' rm$plot_maps(
#'   run=4,
#'   focalRegion = "Africa",
#'   width_max = 30,
#'   curvature_matrix = "random"
#' )
#'
#' # Let's crop the map and try another projection ...
#' rm$plot_maps(
#'   run=4,
#'   focalRegion = "Africa",
#'   width_max = 30,
#'   curvature_matrix = "random",
#'   range_lon = c(-120,140),
#'   range_lat = c(-30,60),
#'   projection = equirectangular
#' )
#'   # Finally let's do this as a plot with all maps using the standard R layout options. I will also show each region's diversity metrics only when it is the focal region, using `diversityCirclesFocalOnly = TRUE`.
#' par(mfrow=c(4,4),mar=c(0.2,0.2,1,0.2))
#' rm$plot_maps(
#'   run=4,
#'   #focalRegion = "Africa", # without this argument, will sequentially plot maps for all focal regions
#'   width_max = 30,
#'   curvature_matrix = "random",
#'   range_lon = c(-120,140),
#'   range_lat = c(-30,60),
#'   projection = equirectangular,
#'   diversityCirclesFocalOnly = TRUE
#' )
#' # clear plots and reset parameters
#' dev.off()
#'
#' ####################################
#' #### More functions and details ####
#' ####################################
#'
#' # Extract info from the object
#' rm$distance_matrix[1:10,1:10] # get back the distance matrix
#' rm$region_table # get back the region table
#' rm$mds # get back the mds data
#'
#' # Diagnostic tools -- the distribution of distances between/within different regions, and how they relate to the distribution of H
#' rm$plot_distance_densities(HdistFromRun=4) # will plot all sequentially including H distribution as used in run 4
#' rm$plot_distance_densities( # overlay a hypothetical H distribution of your choosing
#'   H = 0.01,
#'   H_truncNorm_sd = 0.02,
#'   H_truncNorm_lims = c(.005,.02)
#' )
#' rm$plot_distance_densities( # as above, but with a constant H
#'   H = 0.01
#' )
#'
#' # We can plot MDS plots as we perform a run, and draw hulls surrounding the clusters. This is great for giving an impression of how the subsampling level looks and whether there are any peculiarities interfering with the clustering steps.
#' # Note for this run we have chosen a constant H, as shown in the last plot
#' # Once you have seen enough iterations, hit 'f' to finish the run without plotting.
#' rm$run(
#'   iterations = 100,
#'   H = 0.01,
#'   diagnosticPlotMDSclusters = TRUE
#' )
#'
#' # That's essentially all the package's functionality, though there are some more details in the documentation for the specific methods.
#'
#'
#' @export
ReMIXTURE <- R6::R6Class("ReMIXTURE",

  #### PRIVATE ####
  private = list(

    m = matrix(), # a distance matrix with rownames and colnames giving regions
    rt = data.table(), # region table, an info table with columns "region", "lat" , "long" , and optionally "colour"

    results = NULL,
    mdsPlot = NULL,

    runflag = FALSE, #flags a run has been done, results computed and saved
    plotflag = FALSE, #flags a plot has been done, results computed and saved

    #### VALIDATORS ####
    validate_rt = function(in_rt){
      if( !(data.table::is.data.table(in_rt) || is.data.frame(in_rt)) ){
        stop("Region position table must be a data.frame or data.table")
      }
      in_rt <- data.table::as.data.table(in_rt)
      if( any(!c("region","lon","lat") %in% colnames(in_rt) ) ){
        stop("Region position table must include columns named \"lon\", \"lat\", and \"region\".")
      }
      #check types of cols
      if(!(is.numeric(in_rt$lon) | is.integer(in_rt$lon))){
        stop("Column `lon` must contain numeric or integer values.")
      }
      if(!(is.numeric(in_rt$lat) | is.integer(in_rt$lat))){
        stop("Column `lat` must contain numeric or integer values.")
      }
      if(!(is.character(in_rt$region))){
        stop("Column `region` must be a character vector.")
      }
      if( !all(in_rt$lon %between% c(-180,180)) ){
        stop("All lon (longitude) values must fall between +/- 180.")
      }
      if( !all(in_rt$lat %between% c(-85,85)) ){
        stop("All lat (latitude) values must fall between +/- 85.")
      }
      in_rt
    },

    validate_m = function(in_dm){
      if( !is.matrix(in_dm) ){
        stop( paste0("Argument to distance_matrix must be a matrix.") )
      }
      if( ncol(in_dm) != nrow(in_dm) ){
        stop( paste0("Argument to distance_matrix must be a square matrix.") )
      }
      if(any(!is.finite(in_dm))){
        stop("All distance matrix entries must be finite")
      }
      if( !all(diag(in_dm)==0) ){
        stop("Self-distance (i.e. distance matrix diagonals) should all be zero")
      }
      if ( !all(in_dm[upper.tri(in_dm)]==t(in_dm)[upper.tri(in_dm)]) ){
        stop("Distance matrix must be symmetric.")
      }
      if(!all(in_dm>=0)){
        stop("All distance matrix entries must be positive")
      }
      if (is.null(colnames(in_dm)) | is.null(rownames(in_dm))){
        stop( "Column and row names of input matrix must provide region information." )
      }
      if( !all(colnames(in_dm) == rownames(in_dm)) ) {
        stop( "Column and row names of input matrix must be the same." )
      }
      if( min(table(colnames(in_dm))) < 20 ) {
        warning( paste0("Some regions contain dangerously low numbers of samples (the smallest group has ",min(table(colnames(in_dm)))," members), which could lead to dodgy results. Careful!") )
      }
    },

    validate_dm_rt = function(in_dm,in_rt){
      if(! all(colnames(in_dm) %in% in_rt$region)){
        stop("All regions described in the column / row names of the input distance matrix must correspond to entries in the region information table.")
      }
    },

    validate_variant_sample_table = function(in_sample_table){
      if(is.null(in_sample_table)){
        stop("`sample_table` must be provided when using `vcf_file`.")
      }
      if(!data.table::is.data.table(in_sample_table) && !is.data.frame(in_sample_table)){
        stop("`sample_table` must be a data.frame or data.table.")
      }

      in_sample_table <- data.table::as.data.table(in_sample_table)

      if(!all(c("sample", "region") %in% colnames(in_sample_table))){
        stop("`sample_table` must contain columns named `sample` and `region`.")
      }
      if(!is.character(in_sample_table$sample)){
        stop("Column `sample` in `sample_table` must be character.")
      }
      if(!is.character(in_sample_table$region)){
        stop("Column `region` in `sample_table` must be character.")
      }
      if(any(is.na(in_sample_table$sample)) || any(in_sample_table$sample == "")){
        stop("Column `sample` in `sample_table` must not contain missing or empty values.")
      }
      if(any(is.na(in_sample_table$region)) || any(in_sample_table$region == "")){
        stop("Column `region` in `sample_table` must not contain missing or empty values.")
      }
      if(anyDuplicated(in_sample_table$sample)){
        stop("`sample_table$sample` must not contain duplicates.")
      }

      in_sample_table
    },

    filter_invalid_distance_samples = function(dm){
      invalid_message <- paste0(
        "Fewer than 2 samples remain after filtering invalid distances. ",
        "This often means the VCF has too much missing data. ",
        "Consider building the distance matrix yourself, performing missing-data handling or imputation first, or using data with fewer missing values."
      )

      if(nrow(dm) < 2){
        stop(invalid_message)
      }

      #browser()
      invalid <- !isBehaved(dm)
      drop_all_invalid <- apply(invalid, 1, all)
      if(any(drop_all_invalid)){
        dm <- dm[!drop_all_invalid, !drop_all_invalid, drop = FALSE]
      }
      if(nrow(dm) < 2){
        stop(invalid_message)
      }
      invalid <- !isBehaved(dm)
      #diag(invalid) <- FALSE
      drop_any_invalid <- apply(invalid, 1, any)
      if(any(drop_any_invalid)){
        dm <- dm[!drop_any_invalid, !drop_any_invalid, drop = FALSE]
      }

      if(nrow(dm) < 2){
        stop(invalid_message)
      }

      dm
    },

    distance_matrix_from_vcf_inputs = function(vcf_file, sample_table){
      elapsed_string <- function(start_time){
        sprintf("%.2f sec", proc.time()[["elapsed"]] - start_time)
      }

      total_start <- proc.time()[["elapsed"]]
      snprelate_verbose <- interactive()

      if(!requireNamespace("SNPRelate", quietly = TRUE)){
        stop("Reading VCF input requires the 'SNPRelate' package. Install it from Bioconductor and try again.")
      }
      if(!is.character(vcf_file) || length(vcf_file) != 1L || is.na(vcf_file) || vcf_file == ""){
        stop("`vcf_file` must be a length-1 character path.")
      }
      if(!file.exists(vcf_file)){
        stop("`vcf_file` does not exist: ", vcf_file)
      }
      if(grepl("\\.bcf$", vcf_file, ignore.case = TRUE)){
        stop("Only .vcf and .vcf.gz input files are supported. Please convert the BCF file to VCF or VCF.GZ before using ReMIXTURE.")
      }
      if(!grepl("\\.vcf(\\.gz)?$", vcf_file, ignore.case = TRUE)){
        stop("Only .vcf and .vcf.gz input files are supported.")
      }

      stage_start <- proc.time()[["elapsed"]]
      ce("\tValidating sample table and ...")
      st <- private$validate_variant_sample_table(sample_table)
      ce("\tValidating sample table ... done (", elapsed_string(stage_start), ")")

      gds_tmp_fname <- tempfile(fileext = ".gds")
      gds_handle <- NULL

      on.exit({
        if(!is.null(gds_handle)){
          tryCatch(
            SNPRelate::snpgdsClose(gds_handle),
            error = function(e) NULL
          )
        }
        if(file.exists(gds_tmp_fname)){
          unlink(gds_tmp_fname)
        }
      }, add = TRUE)

      stage_start <- proc.time()[["elapsed"]]
      ce("\tConverting VCF to temporary GDS ...")
      tryCatch(
        SNPRelate::snpgdsVCF2GDS(
          vcf_file,
          out.fn = gds_tmp_fname,
          method = "biallelic.only",
          verbose = snprelate_verbose
        ),
        error = function(e){
          stop("VCF could not be converted to GDS: ", conditionMessage(e))
        }
      )
      ce("\tConverting VCF to temporary GDS ... done (", elapsed_string(stage_start), ")")

      gds_handle <- tryCatch(
        SNPRelate::snpgdsOpen(gds_tmp_fname),
        error = function(e){
          stop("VCF was converted to GDS, but the temporary GDS file could not be opened: ", conditionMessage(e))
        }
      )

      stage_start <- proc.time()[["elapsed"]]
      ce("\tComputing IBS matrix ...")
      ibs <- tryCatch(
        SNPRelate::snpgdsIBS(
          gds_handle,
          sample.id = st$sample,
          missing.rate = NaN,
          remove.monosnp = FALSE,
          num.thread = NA_real_,
          verbose = snprelate_verbose
        ),
        error = function(e){
          stop("IBS calculation failed for the provided VCF: ", conditionMessage(e))
        }
      )
      ce("\tComputing IBS matrix ... done (", elapsed_string(stage_start), ")")

      stage_start <- proc.time()[["elapsed"]]
      ce("\tConverting IBS to distance matrix ...")
      sample_list <- ibs$sample.id
      missing_samples <- setdiff(st$sample, sample_list)
      if(length(missing_samples) > 0){
        stop(
          "`sample_table` contains sample IDs that could not be mapped to VCF samples. Examples: ",
          paste(utils::head(missing_samples, 5), collapse = ", ")
        )
      }

      region_labels <- st$region[match(sample_list, st$sample)]

      if(anyNA(region_labels)){
        stop(
          "Internal mapping error: VCF samples were returned by SNPRelate, but some could not be mapped back to `sample_table$region`."
        )
      }

      dm <- 1 - ibs$ibs
      colnames(dm) <- rownames(dm) <- region_labels
      ce("\tConverting IBS to distance matrix ... done (", elapsed_string(stage_start), ")")

      stage_start <- proc.time()[["elapsed"]]
      ce("\tFiltering invalid distances ...")
      dm <- private$filter_invalid_distance_samples(dm)
      ce("\tFiltering invalid distances ... done (", elapsed_string(stage_start), ")")

      ce("\tDone. Total elapsed time: ", elapsed_string(total_start))

      dm
    },

    normalize_main_panel = function(main_panel){
      if(!(is.numeric(main_panel) || is.integer(main_panel)) || length(main_panel) != 4){
        stop("`main_panel` must be a numeric vector of length 4.")
      }

      panel_vals <- if(is.null(names(main_panel))){
        as.numeric(main_panel)
      } else {
        required_names <- c("left", "right", "bottom", "top")
        if(!all(required_names %in% names(main_panel))){
          stop("Named `main_panel` must include `left`, `right`, `bottom`, and `top`.")
        }
        as.numeric(main_panel[required_names])
      }
      names(panel_vals) <- c("left", "right", "bottom", "top")

      if(any(!is.finite(panel_vals))){
        stop("`main_panel` values must be finite.")
      }
      if(any(panel_vals < 0 | panel_vals > 1)){
        stop("`main_panel` values must fall within [0,1].")
      }
      if(panel_vals["left"] >= panel_vals["right"]){
        stop("`main_panel` must satisfy `left < right`.")
      }
      if(panel_vals["bottom"] >= panel_vals["top"]){
        stop("`main_panel` must satisfy `bottom < top`.")
      }

      panel_vals
    },

    make_curvature_matrix = function(curvature_matrix, regions){
      n_regions <- length(regions)

      if(is.null(curvature_matrix)){
        return(matrix(0.0, nrow = n_regions, ncol = n_regions, dimnames = list(regions, regions)))
      }

      if(is.character(curvature_matrix) && length(curvature_matrix) >= 1L && curvature_matrix[1] == "random"){
        cm <- matrix(rnorm(n_regions^2, 0, 0.3), nrow = n_regions)
        rownames(cm) <- colnames(cm) <- regions
        return(cm)
      }

      if(!is.matrix(curvature_matrix)){
        stop("`curvature_matrix` must be NULL, \"random\", or a matrix.")
      }
      if(!(is.numeric(curvature_matrix) || is.integer(curvature_matrix))){
        stop("`curvature_matrix` must be numeric.")
      }

      rn <- rownames(curvature_matrix)
      cn <- colnames(curvature_matrix)
      has_row_names <- !is.null(rn)
      has_col_names <- !is.null(cn)

      if(has_row_names || has_col_names){
        if(!(has_row_names && has_col_names)){
          stop("Named `curvature_matrix` must have both rownames and colnames.")
        }
        if(!all(regions %in% rn) || !all(regions %in% cn)){
          stop("Named `curvature_matrix` must include all regions in both rownames and colnames.")
        }
        return(curvature_matrix[regions, regions, drop = FALSE])
      }

      if(nrow(curvature_matrix) != n_regions || ncol(curvature_matrix) != n_regions){
        stop("Unnamed `curvature_matrix` must have dimensions length(regions) x length(regions).")
      }

      cm <- curvature_matrix
      rownames(cm) <- colnames(cm) <- regions
      cm
    },

    get_map_plot_state = function(run, width_max, alpha_max, curvature_matrix){
      if(private$runflag == FALSE){
        stop("Analysis has not been run. Perform using `$run()`")
      }
      if(is.null(run) & length(private$results) > 1){
        stop("Please provide a run number to plot from (consider using `<ReMIXTURE Object>$plot_h_optimisation()`, `<ReMIXTURE Object>$plot_results_grid()`, and `<ReMIXTURE Object>$plot_clustercounts()` to assess which run parameters are appropriate).")
      }
      if(is.null(run) & length(private$results) == 1){
        run <- 1
      }

      st <- private$results[[run]]$overlap
      rt <- data.table(
        region = colnames(st),
        totDiv = private$results[[run]]$diversity
      )

      rt[, uniqueDiv := private$results[[run]]$overlap[r, r], by = .(r = region)]
      rt <- private$rt[rt, on = .(region)]

      width_scaling <- compute_plot_width_scaling(
        totDiv = rt$totDiv,
        uniqueDiv = rt$uniqueDiv,
        overlap = st,
        width_max = width_max
      )
      rt[, wTotDiv := width_scaling$wTotDiv]
      rt[, wUniqueDiv := width_scaling$wUniqueDiv]
      wst <- width_scaling$wst

      ct <- private$make_curvature_matrix(curvature_matrix, rt$region)

      at <- copy(st)
      diag(at) <- 0.0
      if(is.null(alpha_max)){
        at[,] <- 1.0
      } else {
        maxAlphaSource <- suppressWarnings(max(at, na.rm = TRUE))
        if(!is.finite(maxAlphaSource) || maxAlphaSource <= 0){
          at[,] <- 0.0
        } else {
          at[,] <- at / maxAlphaSource * alpha_max
        }
      }

      list(
        run = run,
        rt = rt,
        regions = rt$region,
        st = st,
        wst = wst,
        at = at,
        ct = ct
      )
    },

    validate_focal_regions = function(focalRegions, regions){
      if(is.null(focalRegions)){
        return(regions)
      }
      if(!is.character(focalRegions)){
        stop("`focalRegions` must be a character vector.")
      }

      focalRegions <- focalRegions[!duplicated(focalRegions)]
      missing_regions <- setdiff(focalRegions, regions)
      if(length(missing_regions) > 0){
        stop(
          "Unknown focal region(s): ",
          paste(missing_regions, collapse = ", ")
        )
      }

      focalRegions
    },

    draw_region_diversity_circle = function(trt, j, projection, shadow = FALSE, circle_scale = 1){
      cdt <- circle_seg(trt$lon[j], trt$lat[j], radius = trt$wTotDiv[j] / 2 * circle_scale) %>% mat2dtLL()
      udt <- circle_seg(trt$lon[j], trt$lat[j], radius = trt$wUniqueDiv[j] / 2 * circle_scale) %>% mat2dtLL()

      if(shadow){
        plotMapItem(cdt, projFun = projection, plotFun = polygon, col = "#00000055")
      }
      plotMapItem(cdt, projFun = projection, plotFun = polygon, col = "#000000FF")
      plotMapItem(udt, projFun = projection, plotFun = polygon, col = "#FFFFFF")
    },

    get_region_anchors = function(trt, projection){
      projected <- projection(
        dtLL = trt[, .(lon, lat)],
        projColNames = c("x", "y")
      )

      data.table(
        region = trt$region,
        x = grconvertX(projected$x, from = "user", to = "ndc"),
        y = grconvertY(projected$y, from = "user", to = "ndc")
      )
    },

    filter_visible_anchors = function(anchors, main_panel, tolerance = 0.005){
      main_panel <- private$normalize_main_panel(main_panel)
      anchors <- as.data.table(copy(anchors))

      anchors[
        x >= main_panel["left"] - tolerance &
          x <= main_panel["right"] + tolerance &
          y >= main_panel["bottom"] - tolerance &
          y <= main_panel["top"] + tolerance
      ][]
    },

    compute_inset_height = function(inset_width, inset_height, range_lon, range_lat, panel_aspect = NULL){
      if(!(is.numeric(inset_width) || is.integer(inset_width)) || length(inset_width) != 1 || !is.finite(inset_width) || inset_width <= 0){
        stop("`inset_width` must be a single positive finite number.")
      }

      if(!is.null(inset_height)){
        if(!(is.numeric(inset_height) || is.integer(inset_height)) || length(inset_height) != 1 || !is.finite(inset_height) || inset_height <= 0){
          stop("`inset_height` must be NULL or a single positive finite number.")
        }
        return(as.numeric(inset_height))
      }

      lon_span <- abs(diff(range_lon))
      lat_span <- abs(diff(range_lat))
      if(!is.finite(lon_span) || !is.finite(lat_span) || lon_span <= 0 || lat_span <= 0){
        stop("`range_lon` and `range_lat` must each define a non-zero finite range.")
      }

      if(!is.null(panel_aspect)){
        if(!(is.numeric(panel_aspect) || is.integer(panel_aspect)) || length(panel_aspect) != 1 || !is.finite(panel_aspect) || panel_aspect <= 0){
          stop("`panel_aspect` must be NULL or a single positive finite number.")
        }
        target_aspect <- as.numeric(panel_aspect)
      } else {
        raw_aspect <- lon_span / lat_span
        target_aspect <- max(2.2, min(2.8, raw_aspect))
      }
      device <- par("din")
      computed_height <- inset_width * device[1] / (target_aspect * device[2])
      # Keep automatic inset maps from collapsing into very flat strips on wide devices.
      computed_height <- max(0.13, min(0.24, computed_height))

      computed_height
    },

    validate_composite_geometry = function(main_panel, inset_width, gap, outer_margin){
      main_panel <- private$normalize_main_panel(main_panel)
      left_capacity <- unname(main_panel["left"] - gap - outer_margin)
      right_capacity <- unname(1 - outer_margin - main_panel["right"] - gap)

      list(
        left_capacity = left_capacity,
        right_capacity = right_capacity
      )
    },

    rectangles_overlap = function(a_left, a_right, a_bottom, a_top, b_left, b_right, b_bottom, b_top, tolerance = 1e-6){
      (a_right > b_left + tolerance) &&
        (a_left < b_right - tolerance) &&
        (a_top > b_bottom + tolerance) &&
        (a_bottom < b_top - tolerance)
    },

    render_map_panel = function(
      state,
      focalRegion,
      range_lon,
      range_lat,
      projection,
      mapData,
      overview = FALSE,
      diversityCirclesFocalOnly = FALSE,
      returnAnchors = FALSE,
      circle_scale = 1,
      line_scale = 1,
      targetRegions = NULL,
      clip_to_panel = FALSE,
      border_last = FALSE,
      tight_axes = FALSE,
      draw_grid = TRUE,
      draw_title = TRUE,
      title_cex = NULL,
      title_line = NULL,
      inset_style = "default"
    ){
      plotMiddle <- findCentreLL(range_lon, range_lat)
      trt <- copy(state$rt) %>% rotateLatLonDtLL(-plotMiddle[1], -plotMiddle[2], splitPlotGrps = FALSE)
      use_line_emphasis <- !overview && identical(inset_style, "line_emphasis")
      effective_circle_scale <- if(use_line_emphasis) circle_scale * 0.35 else circle_scale
      effective_line_scale <- if(use_line_emphasis) line_scale * 1.8 else line_scale
      draw_focal_circle_last <- TRUE
      curved_line_path <- function(x1, y1, x2, y2, curvature = 0, n = 200){
        if(curvature == 0){
          return(filled_line(x1, y1, x2, y2, n = n))
        }
        if(!(curvature > -pi & curvature < pi)){
          stop("Curvatures must be in [-pi,+pi]")
        }
        reflectX <- sign(curvature) == 1
        curvature <- abs(curvature)
        l <- euc_dist(x1, y1, x2, y2)
        x_c <- (l / 2) / tan(curvature)
        h_c <- (l / 2) / sin(pi - curvature)
        circle_seg(
          x_c,
          l / 2,
          h_c,
          (pi * (3 / 2)) - curvature,
          (pi * (3 / 2)) + curvature,
          n = n
        ) %>%
          reflect(about_x = reflectX, about_y = FALSE) %>%
          rotate(angle(x1, y1, x2, y2)) %>%
          translate(by_x = x1, by_y = y1)
      }

      pe <- plotEmptyMap(
        range_lon,
        range_lat,
        projFun = projection,
        xaxs = if(tight_axes) "i" else "r",
        yaxs = if(tight_axes) "i" else "r"
      )
      if(draw_grid){
        plotMapItem(makeMapDataLatLonLines(), range_lon, range_lat, projFun = projection, plotFun = lines, col = "#00000022", lwd = 0.4)
      }
      plotMapItem(mapData, range_lon, range_lat, projFun = projection, plotFun = polygon, col = "#f7bf25", lwd = 0.2)
      if(!border_last){
        plotMapBorder(range_lon, range_lat, projFun = projection, plotEdges = pe, lwd = 4)
      }
      if(isTRUE(clip_to_panel)){
        usr <- par("usr")
        clip(usr[1], usr[2], usr[3], usr[4])
      }

      if(overview){
        for(j in seq_len(nrow(trt))){
          private$draw_region_diversity_circle(trt, j, projection, shadow = TRUE, circle_scale = circle_scale)
        }
        if(border_last){
          plotMapBorder(range_lon, range_lat, projFun = projection, plotEdges = pe, lwd = 4)
        }
        if(returnAnchors){
          return(private$get_region_anchors(trt, projection))
        }
        return(invisible(NULL))
      }

      i <- match(focalRegion, trt$region)
      if(is.na(i)){
        stop("`focalRegion` must name a valid region.")
      }

      targetRegions <- if(is.null(targetRegions)) state$regions else targetRegions
      targetRegions <- targetRegions[targetRegions %in% trt$region]
      targetIdx <- which(trt$region %in% targetRegions)

      if(use_line_emphasis){
        private$draw_region_diversity_circle(trt, i, projection, shadow = FALSE, circle_scale = effective_circle_scale)
        draw_focal_circle_last <- FALSE
      }

      for(j in targetIdx){
        if(i == j){ next }
        raw_curve_width <- state$wst[trt$region[i], trt$region[j]] * effective_line_scale
        if(use_line_emphasis){
          if(!is.finite(raw_curve_width) || raw_curve_width <= 0){
            next
          }
          ldt <- curved_line_path(
            x1 = trt$lon[i], y1 = trt$lat[i],
            x2 = trt$lon[j], y2 = trt$lat[j],
            curvature = state$ct[trt$region[i], trt$region[j]]
          ) %>% mat2dtLL()
          plotMapItem(
            ldt,
            projFun = projection,
            plotFun = lines,
            col = "grey20",
            lty = 1,
            lwd = pmax(raw_curve_width * 2, 0.65)
          )
        } else {
          ldt <- curved_rounded_line(
            x1 = trt$lon[i], y1 = trt$lat[i],
            x2 = trt$lon[j], y2 = trt$lat[j],
            width = raw_curve_width,
            curvature = state$ct[trt$region[i], trt$region[j]]
          ) %>% mat2dtLL()
          plotMapItem(
            ldt,
            projFun = projection,
            plotFun = polygon,
            col = alpha("black", state$at[trt$region[i], trt$region[j]]),
            border = "#000000",
            lwd = 0.15
          )
        }
        if(diversityCirclesFocalOnly == FALSE){
          private$draw_region_diversity_circle(trt, j, projection, shadow = TRUE, circle_scale = effective_circle_scale)
        }
      }

      if(draw_focal_circle_last){
        private$draw_region_diversity_circle(trt, i, projection, shadow = FALSE, circle_scale = effective_circle_scale)
      }
      if(border_last){
        plotMapBorder(range_lon, range_lat, projFun = projection, plotEdges = pe, lwd = 4)
      }
      if(draw_title){
        if(is.null(title_cex) && is.null(title_line)){
          title(main = trt$region[i])
        } else {
          title_args <- list(main = trt$region[i])
          if(!is.null(title_cex)){
            title_args$cex.main <- title_cex
          }
          if(!is.null(title_line)){
            title_args$line <- title_line
          }
          do.call(title, title_args)
        }
      }

      if(returnAnchors){
        return(private$get_region_anchors(trt, projection))
      }
      invisible(NULL)
    },

    make_composite_layout = function(anchors, focalRegions, main_panel, inset_width, inset_height){
      main_panel <- private$normalize_main_panel(main_panel)
      anchors <- as.data.table(copy(anchors))
      anchors <- anchors[region %in% focalRegions]
      if(nrow(anchors) != length(focalRegions)){
        stop("Anchors are missing for one or more focal regions.")
      }

      anchors <- anchors[match(focalRegions, region)]
      cx <- mean(c(main_panel["left"], main_panel["right"]))
      cy <- mean(c(main_panel["bottom"], main_panel["top"]))
      outer_margin <- 0.03
      gap <- 0.045
      label_strip <- 0.025
      min_card_gap <- 0.025

      anchors[, side := fifelse(
        abs(x - cx) >= abs(y - cy),
        fifelse(x >= cx, "right", "left"),
        fifelse(y >= cy, "top", "bottom")
      )]

      distribute_intervals <- function(n, size, min_val, max_val, min_gap){
        span <- max_val - min_val
        if(n <= 0){
          return(list(starts = numeric(), size = size))
        }
        if(span <= 0){
          return(list(starts = rep(min_val, n), size = size))
        }
        adj_size <- size
        if(n == 1){
          start <- min_val + max(0, (span - adj_size) / 2)
          return(list(starts = start, size = adj_size))
        }
        effective_gap <- min_gap
        needed <- n * adj_size + (n - 1) * effective_gap
        if(needed > span){
          effective_gap <- max(0, (span - n * adj_size) / (n - 1))
          needed <- n * adj_size + (n - 1) * effective_gap
        }
        start0 <- min_val + max(0, (span - needed) / 2)
        starts <- start0 + (seq_len(n) - 1) * (adj_size + effective_gap)
        list(starts = starts, size = adj_size)
      }

      top_dt <- anchors[side == "top"][order(x)]
      bottom_dt <- anchors[side == "bottom"][order(x)]
      left_dt <- anchors[side == "left"][order(-y)]
      right_dt <- anchors[side == "right"][order(-y)]

      geometry <- private$validate_composite_geometry(
        main_panel = main_panel,
        inset_width = inset_width,
        gap = gap,
        outer_margin = outer_margin
      )
      if(nrow(left_dt) > 0 && inset_width > geometry$left_capacity){
        stop(
          "Insufficient horizontal space for left-side inset panels. ",
          "Current inset_width = ", signif(inset_width, 4),
          ", available left capacity = ", signif(geometry$left_capacity, 4),
          ". Reduce inset_width or move main_panel['left'] to the right."
        )
      }
      if(nrow(right_dt) > 0 && inset_width > geometry$right_capacity){
        stop(
          "Insufficient horizontal space for right-side inset panels. ",
          "Current inset_width = ", signif(inset_width, 4),
          ", available right capacity = ", signif(geometry$right_capacity, 4),
          ". Reduce inset_width or move main_panel['right'] to the left."
        )
      }

      map_width <- inset_width
      map_height <- inset_height
      card_width <- map_width
      card_height <- map_height + label_strip
      side_map_width <- inset_width
      side_map_height <- inset_height
      side_card_width <- side_map_width
      side_card_height <- side_map_height + label_strip
      horizontal_capacity <- floor(((1 - 2 * outer_margin) + min_card_gap) / (card_width + min_card_gap))
      if(nrow(top_dt) > horizontal_capacity){
        stop(
          "Insufficient horizontal space for top-side inset panels. ",
          "Need to place ", nrow(top_dt), " panel(s) with inset_width = ", signif(inset_width, 4),
          " and label strip height = ", signif(label_strip, 4),
          ". Reduce focalRegions, reduce inset_width, or edit returned panels manually."
        )
      }
      if(nrow(bottom_dt) > horizontal_capacity){
        stop(
          "Insufficient horizontal space for bottom-side inset panels. ",
          "Need to place ", nrow(bottom_dt), " panel(s) with inset_width = ", signif(inset_width, 4),
          " and label strip height = ", signif(label_strip, 4),
          ". Reduce focalRegions, reduce inset_width, or edit returned panels manually."
        )
      }

      top_card_bottom <- max(
        outer_margin,
        min(1 - outer_margin - card_height, main_panel["top"] + gap)
      )
      bottom_card_bottom <- max(
        outer_margin,
        min(1 - outer_margin - card_height, main_panel["bottom"] - gap - card_height)
      )
      bottom_card_top <- bottom_card_bottom + card_height

      left_right_y_min <- outer_margin
      left_right_y_max <- 1 - outer_margin
      if(nrow(top_dt) > 0){
        left_right_y_max <- min(left_right_y_max, top_card_bottom - gap)
      }
      if(nrow(bottom_dt) > 0){
        left_right_y_min <- max(left_right_y_min, bottom_card_top + gap)
      }
      if(left_right_y_max <= left_right_y_min){
        stop(
          "Insufficient vertical space for side inset panels after reserving top/bottom panel bands. ",
          "Reduce focalRegions, reduce inset_height, or adjust main_panel."
        )
      }
      vertical_capacity <- floor(((left_right_y_max - left_right_y_min) + min_card_gap) / (side_card_height + min_card_gap))
      if(nrow(left_dt) > vertical_capacity){
        stop(
          "Insufficient vertical space for left-side inset panels. ",
          "Need to place ", nrow(left_dt), " panel(s) with inset_height = ", signif(inset_height, 4),
          " and label strip height = ", signif(label_strip, 4),
          ". Reduce focalRegions, reduce inset_height, or edit returned panels manually."
        )
      }
      if(nrow(right_dt) > vertical_capacity){
        stop(
          "Insufficient vertical space for right-side inset panels. ",
          "Need to place ", nrow(right_dt), " panel(s) with inset_height = ", signif(inset_height, 4),
          " and label strip height = ", signif(label_strip, 4),
          ". Reduce focalRegions, reduce inset_height, or edit returned panels manually."
        )
      }

      build_side_layout <- function(side_name){
        dt <- switch(
          side_name,
          top = copy(top_dt),
          bottom = copy(bottom_dt),
          left = copy(left_dt),
          right = copy(right_dt)
        )
        if(nrow(dt) == 0){
          return(NULL)
        }

        if(side_name %in% c("top", "bottom")){
          x_layout <- distribute_intervals(
            n = nrow(dt),
            size = card_width,
            min_val = outer_margin,
            max_val = 1 - outer_margin,
            min_gap = min_card_gap
          )
          card_left <- x_layout$starts
          card_right <- card_left + x_layout$size
          map_width_side <- x_layout$size
          if(side_name == "top"){
            card_bottom <- rep(top_card_bottom, nrow(dt))
            card_top <- card_bottom + card_height
            map_left <- card_left
            map_right <- card_right
            map_bottom <- card_bottom
            map_top <- map_bottom + map_height
          } else {
            card_bottom <- rep(bottom_card_bottom, nrow(dt))
            card_top <- card_bottom + card_height
            map_left <- card_left
            map_right <- card_right
            map_top <- card_top
            map_bottom <- map_top - map_height
          }
          out <- data.table(
            region = dt$region,
            side = side_name,
            card_left = card_left,
            card_right = card_right,
            card_bottom = card_bottom,
            card_top = card_top,
            map_left = map_left,
            map_right = map_right,
            map_bottom = map_bottom,
            map_top = map_top
          )
        } else {
          y_layout <- distribute_intervals(
            n = nrow(dt),
            size = side_card_height,
            min_val = left_right_y_min,
            max_val = left_right_y_max,
            min_gap = min_card_gap
          )
          card_bottom <- y_layout$starts
          card_top <- card_bottom + y_layout$size
          map_height_side <- side_map_height
          card_left <- rep(
            if(side_name == "right"){
              main_panel["right"] + gap
            } else {
              main_panel["left"] - gap - side_card_width
            },
            nrow(dt)
          )
          card_right <- card_left + side_card_width
          map_left <- card_left
          map_right <- card_right
          map_bottom <- card_bottom
          map_top <- map_bottom + map_height_side
          out <- data.table(
            region = dt$region,
            side = side_name,
            card_left = card_left,
            card_right = card_right,
            card_bottom = card_bottom,
            card_top = card_top,
            map_left = map_left,
            map_right = map_right,
            map_bottom = map_bottom,
            map_top = map_top
          )
        }

        labels <- private$compute_panel_label_positions(out)
        out <- labels[out, on = .(region)]
        out[, `:=`(
          left = map_left,
          right = map_right,
          bottom = map_bottom,
          top = map_top
        )]
        out
      }

      panels <- rbindlist(
        lapply(c("top", "right", "bottom", "left"), build_side_layout),
        use.names = TRUE,
        fill = TRUE
      )

      if(nrow(panels) > 0){
        overlaps_main_idx <- which(
          panels$map_right > main_panel["left"] + 1e-6 &
            panels$map_left < main_panel["right"] - 1e-6 &
            panels$map_top > main_panel["bottom"] + 1e-6 &
            panels$map_bottom < main_panel["top"] - 1e-6
        )
        if(length(overlaps_main_idx) > 0){
          stop(
            paste0(
              "Automatic layout overlaps the central overview for panel(s): ",
              paste(panels$region[overlaps_main_idx], collapse = ", "),
              ". Reduce `inset_width` or narrow `main_panel`."
            )
          )
        }
      }

      if(nrow(panels) > 1){
        for(i in seq_len(nrow(panels) - 1)){
          for(j in (i + 1):nrow(panels)){
            if(private$rectangles_overlap(
              panels$map_left[i], panels$map_right[i], panels$map_bottom[i], panels$map_top[i],
              panels$map_left[j], panels$map_right[j], panels$map_bottom[j], panels$map_top[j]
            )){
              stop(
                paste0(
                  "Automatic layout produced overlapping panels: ",
                  panels$region[i], " and ", panels$region[j],
                  ". Consider passing fewer focalRegions or editing returned panels."
                )
              )
            }
          }
        }
      }

      panels[match(focalRegions, region)][]
    },

    make_ellipse_layout = function(visible_anchors, focalRegions, main_panel, inset_width, inset_height){
      main_panel <- private$normalize_main_panel(main_panel)
      anchors <- as.data.table(copy(visible_anchors))
      anchors <- anchors[region %in% focalRegions]
      if(nrow(anchors) != length(focalRegions)){
        stop("Anchors are missing for one or more focal regions.")
      }
      anchors <- anchors[match(focalRegions, region)]

      cx <- mean(c(main_panel["left"], main_panel["right"]))
      cy <- mean(c(main_panel["bottom"], main_panel["top"]))
      main_w <- main_panel["right"] - main_panel["left"]
      main_h <- main_panel["top"] - main_panel["bottom"]
      main_half_w <- main_w / 2
      main_half_h <- main_h / 2
      panel_half_w <- inset_width / 2
      panel_half_h <- inset_height / 2

      outer_margin <- 0.03
      label_strip <- 0.025
      subPlotMult <- max(inset_width / main_w, inset_height / main_h)
      unitCircStretch <- sqrt(2) + sqrt((subPlotMult / 2)^2 + (subPlotMult / 2)^2)
      distAdjFactor <- 1.0
      n <- nrow(anchors)

      norm_angle <- function(x){
        x %% (2 * pi)
      }
      angle_dist <- function(a, b){
        abs(atan2(sin(a - b), cos(a - b)))
      }

      candidate_base <- seq(0, 2 * pi * (n - 1) / n, length.out = n)
      fallback_thetas <- candidate_base

      anchor_theta <- atan2(
        (anchors$y - cy) / main_half_h,
        (anchors$x - cx) / main_half_w
      )
      anchor_theta <- norm_angle(anchor_theta)
      near_centre <- abs((anchors$x - cx) / main_half_w) < 1e-8 &
        abs((anchors$y - cy) / main_half_h) < 1e-8
      if(any(near_centre)){
        anchor_theta[near_centre] <- fallback_thetas[which(near_centre)]
      }

      anchor_order <- order(anchor_theta, anchors$region)
      anchor_theta_sorted <- anchor_theta[anchor_order]
      anchors_sorted <- anchors[anchor_order]

      offset_candidates <- anchor_theta_sorted - candidate_base
      best_score <- Inf
      best_offset <- offset_candidates[1]
      best_candidate_theta <- NULL
      best_candidate_order <- NULL

      for(offset in offset_candidates){
        candidate_theta <- norm_angle(candidate_base + offset)
        candidate_order <- order(candidate_theta)
        score <- sum(angle_dist(anchor_theta_sorted, candidate_theta[candidate_order]))
        if(score < best_score){
          best_score <- score
          best_offset <- offset
          best_candidate_theta <- candidate_theta
          best_candidate_order <- candidate_order
        }
      }

      assigned_theta <- best_candidate_theta[best_candidate_order]
      px <- cx + main_half_w * unitCircStretch * distAdjFactor * cos(assigned_theta)
      py <- cy + main_half_h * unitCircStretch * distAdjFactor * sin(assigned_theta)

      out <- copy(anchors_sorted)[, `:=`(
        .region_order = match(region, focalRegions),
        .theta = assigned_theta,
        .px = px,
        .py = py
      )]

      out[, `:=`(
        map_left = .px - panel_half_w,
        map_right = .px + panel_half_w,
        map_bottom = .py - panel_half_h,
        map_top = .py + panel_half_h
      )]

      out[, side := fifelse(
        abs(.px - cx) >= abs(.py - cy),
        fifelse(.px >= cx, "right", "left"),
        fifelse(.py >= cy, "top", "bottom")
      )]

      out[, `:=`(
        card_left = map_left,
        card_right = map_right,
        card_bottom = fifelse(side == "bottom", map_bottom - label_strip, map_bottom),
        card_top = fifelse(side == "bottom", map_top, map_top + label_strip)
      )]

      outside_idx <- which(
        out$map_left < 0 | out$map_right > 1 |
          out$map_bottom < 0 | out$map_top > 1 |
          out$card_left < 0 | out$card_right > 1 |
          out$card_bottom < 0 | out$card_top > 1
      )
      if(length(outside_idx) > 0){
        stop(
          "Ellipse layout could not place region `", out$region[outside_idx[1]],
          "` inside the device. Reduce `inset_width`, reduce `inset_height`, adjust `main_panel`, or use `layout = \"slot\"`."
        )
      }

      for(i in seq_len(nrow(out))){
        if(private$rectangles_overlap(
          out$map_left[i], out$map_right[i], out$map_bottom[i], out$map_top[i],
          main_panel["left"], main_panel["right"], main_panel["bottom"], main_panel["top"]
        )){
          stop(
            "Ellipse layout overlaps the central overview for region `", out$region[i],
            "`. Reduce `inset_width`, reduce `inset_height`, adjust `main_panel`, or use `layout = \"slot\"`."
          )
        }
      }

      out <- out[, .(
        region,
        side,
        card_left,
        card_right,
        card_bottom,
        card_top,
        map_left,
        map_right,
        map_bottom,
        map_top
      )]
      labels <- private$compute_panel_label_positions(out)
      out <- labels[out, on = .(region)]
      out[, `:=`(
        left = map_left,
        right = map_right,
        bottom = map_bottom,
        top = map_top
      )]

      if(nrow(out) > 1){
        for(i in seq_len(nrow(out) - 1)){
          for(j in (i + 1):nrow(out)){
            if(private$rectangles_overlap(
              out$card_left[i], out$card_right[i], out$card_bottom[i], out$card_top[i],
              out$card_left[j], out$card_right[j], out$card_bottom[j], out$card_top[j]
            )){
              stop(
                "Ellipse layout produced overlapping panels: ",
                out$region[i], " and ", out$region[j],
                ". Try fewer `focalRegions`, smaller `inset_width`, or `layout = \"slot\"`."
              )
            }
          }
        }
      }

      out[match(focalRegions, region)][]
    },

    validate_composite_layout = function(panels, focalRegions, main_panel){
      if(!(is.data.frame(panels) || data.table::is.data.table(panels))){
        stop("`panels` must be a data.frame or data.table.")
      }

      panels <- as.data.table(copy(panels))
      main_panel <- private$normalize_main_panel(main_panel)
      label_strip <- 0.025
      outer_margin <- 0.03
      gap <- 0.045
      has_map_schema <- all(c("map_left", "map_right", "map_bottom", "map_top") %in% colnames(panels))
      has_old_schema <- all(c("left", "right", "bottom", "top") %in% colnames(panels))
      if(!has_map_schema && !has_old_schema){
        stop("`panels` must contain either `left/right/bottom/top` or `map_left/map_right/map_bottom/map_top`.")
      }
      if(!has_map_schema){
        panels[, `:=`(
          map_left = left,
          map_right = right,
          map_bottom = bottom,
          map_top = top
        )]
      }
      extra_regions <- setdiff(panels$region, focalRegions)
      if(length(extra_regions) > 0){
        stop(
          "`panels` includes region(s) without visible anchors or outside the selected composite set: ",
          paste(extra_regions, collapse = ", ")
        )
      }

      panels <- panels[region %in% focalRegions]
      if(anyDuplicated(panels$region)){
        stop("`panels` must contain at most one row per region.")
      }
      if(!setequal(panels$region, focalRegions)){
        stop("`panels` must contain exactly one row for each focal region.")
      }

      if(!"side" %in% colnames(panels)){
        cx <- mean(c(main_panel["left"], main_panel["right"]))
        cy <- mean(c(main_panel["bottom"], main_panel["top"]))
        panels[, `:=`(
          center_x = (map_left + map_right) / 2,
          center_y = (map_bottom + map_top) / 2
        )]
        panels[, side := fifelse(
          abs(center_x - cx) >= abs(center_y - cy),
          fifelse(center_x >= cx, "right", "left"),
          fifelse(center_y >= cy, "top", "bottom")
        )]
        panels[, c("center_x", "center_y") := NULL]
      }

      if(!all(c("card_left", "card_right", "card_bottom", "card_top") %in% colnames(panels))){
        panels[, `:=`(
          card_left = map_left,
          card_right = map_right,
          card_bottom = fifelse(side == "bottom", map_bottom - label_strip, map_bottom),
          card_top = fifelse(side == "bottom", map_top, map_top + label_strip)
        )]
      }

      need_labels <- !all(c("label_x", "label_y", "label_adj_x", "label_adj_y") %in% colnames(panels))
      if(need_labels){
        labels <- private$compute_panel_label_positions(panels)
        panels <- labels[panels, on = .(region)]
      }

      overlaps_main_idx <- which(
        panels$map_right > main_panel["left"] + 1e-6 &
          panels$map_left < main_panel["right"] - 1e-6 &
          panels$map_top > main_panel["bottom"] + 1e-6 &
          panels$map_bottom < main_panel["top"] - 1e-6
      )
      if(length(overlaps_main_idx) > 0){
        warning(
          paste0(
            "Custom `panels` overlap the central overview: ",
            paste(panels$region[overlaps_main_idx], collapse = ", "),
            ". Adjust panel coordinates or narrow `main_panel`."
          ),
          call. = FALSE
        )
      }

      left_capacity <- unname(main_panel["left"] - gap - outer_margin)
      right_capacity <- unname(1 - outer_margin - main_panel["right"] - gap)
      left_width_excess <- panels[side == "left", .N] > 0 &&
        panels[side == "left", max(map_right - map_left, na.rm = TRUE)] > left_capacity
      right_width_excess <- panels[side == "right", .N] > 0 &&
        panels[side == "right", max(map_right - map_left, na.rm = TRUE)] > right_capacity
      if(isTRUE(left_width_excess) || isTRUE(right_width_excess)){
        warning(
          paste0(
            "Custom `panels` may overlap the central overview: side-panel width exceeds available side capacity ",
            "(left capacity=",
            signif(left_capacity, 4),
            ", right capacity=",
            signif(right_capacity, 4),
            "). Consider reducing panel width or narrowing `main_panel`."
          ),
          call. = FALSE
        )
      }

      panels[, `:=`(
        left = map_left,
        right = map_right,
        bottom = map_bottom,
        top = map_top
      )]

      numeric_cols <- c(
        "map_left", "map_right", "map_bottom", "map_top",
        "card_left", "card_right", "card_bottom", "card_top",
        "label_x", "label_y", "label_adj_x", "label_adj_y"
      )
      for(col_name in numeric_cols){
        if(any(!is.finite(panels[[col_name]]))){
          stop("`panels` coordinates must be finite.")
        }
      }
      if(any(
        panels$map_left < 0 | panels$map_right > 1 |
          panels$map_bottom < 0 | panels$map_top > 1 |
          panels$card_left < 0 | panels$card_right > 1 |
          panels$card_bottom < 0 | panels$card_top > 1
      )){
        stop("`panels` map/card coordinates must fall within [0,1].")
      }
      if(any(panels$map_left >= panels$map_right)){
        stop("Each `panels` row must satisfy `map_left < map_right`.")
      }
      if(any(panels$map_bottom >= panels$map_top)){
        stop("Each `panels` row must satisfy `map_bottom < map_top`.")
      }
      if(any(panels$card_left > panels$map_left | panels$card_right < panels$map_right)){
        stop("Card rectangles must contain map rectangles horizontally.")
      }
      if(any(panels$card_bottom > panels$map_bottom | panels$card_top < panels$map_top)){
        stop("Card rectangles must contain map rectangles vertically.")
      }

      panels[match(focalRegions, region)][]
    },

    make_leader_segments = function(panels, anchors){
      panels <- as.data.table(copy(panels))
      anchors <- as.data.table(copy(anchors))

      anchor_idx <- match(panels$region, anchors$region)
      if(any(is.na(anchor_idx))){
        stop("Anchors are missing for one or more panels.")
      }

      x <- anchors$x[anchor_idx]
      y <- anchors$y[anchor_idx]
      x0 <- fifelse(
        panels$side %in% c("top", "bottom"),
        (panels$map_left + panels$map_right) / 2,
        fifelse(panels$side == "left", panels$map_right, panels$map_left)
      )
      y0 <- fifelse(
        panels$side %in% c("left", "right"),
        (panels$map_bottom + panels$map_top) / 2,
        fifelse(panels$side == "top", panels$map_bottom, panels$map_top)
      )

      data.table(
        region = panels$region,
        x0 = x0,
        y0 = y0,
        x1 = x,
        y1 = y
      )
    },

    render_leader_segments = function(leader_segments){
      leader_segments <- as.data.table(copy(leader_segments))
      if(nrow(leader_segments) == 0){
        return(invisible(NULL))
      }

      private$begin_ndc_overlay()
      with(leader_segments, segments(x0, y0, x1, y1, lty = 2, lwd = 0.8, col = "#00000066"))

      invisible(NULL)
    },

    begin_ndc_overlay = function(){
      par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE, xpd = NA)
      plot.new()
      plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
      invisible(NULL)
    },

    compute_panel_label_positions = function(
      panels,
      panel_label_cex = 0.75
    ){
      panels <- as.data.table(copy(panels))
      if(nrow(panels) == 0){
        return(data.table(region = character(), x = numeric(), y = numeric(), adj_x = numeric(), adj_y = numeric()))
      }

      label_strip <- 0.025
      if(!all(c("map_left", "map_right", "map_bottom", "map_top") %in% colnames(panels))){
        panels[, `:=`(
          map_left = left,
          map_right = right,
          map_bottom = bottom,
          map_top = top
        )]
      }

      device_size <- par("din")
      panels[, label_width := strwidth(region, cex = panel_label_cex, units = "inches") / device_size[1]]
      panels[, label_height := strheight(region, cex = panel_label_cex, units = "inches") / device_size[2]]

      panels[, `:=`(
        label_x = (map_left + map_right) / 2,
        label_y = fifelse(
          side == "bottom",
          map_bottom - label_strip * 0.5,
          map_top + label_strip * 0.5
        ),
        adj_x = 0.5,
        adj_y = 0.5
      )]

      panels[, label_x := pmax(0.02 + label_width / 2, pmin(0.98 - label_width / 2, label_x))]
      panels[, label_y := pmax(0.02 + label_height / 2, pmin(0.98 - label_height / 2, label_y))]

      panels[, .(
        region,
        label_x,
        label_y,
        label_adj_x = adj_x,
        label_adj_y = adj_y
      )]
    },

    render_panel_labels = function(panels, cex = 0.75){
      panels <- as.data.table(copy(panels))
      if(nrow(panels) == 0){
        return(invisible(NULL))
      }

      private$begin_ndc_overlay()
      for(i in seq_len(nrow(panels))){
        text(
          x = panels$label_x[i],
          y = panels$label_y[i],
          labels = panels$region[i],
          cex = cex,
          adj = c(panels$label_adj_x[i], panels$label_adj_y[i]),
          xpd = NA
        )
      }
      invisible(NULL)
    }
  ),

  #### ACTIVE BINDINGS ####
  active = list(

    #' @field distance_matrix
    #' Return a copy of the distance matrix.
    distance_matrix=function(){
      out_dm <- copy(private$m)
      diag(out_dm) <- 0.0
      return(out_dm)
    },


    #' @field region_table
    #' Return a copy of the region table with lat/lon information.
    region_table=function(){
      rt <- copy(private$rt)
      return(rt)
    },

    #' @field run_results
    #' Return the raw results of ReMIXTURE runs (runs must be first performed with `<ReMIXTURE_object>$run()`)
    run_results=function(){
      if(private$runflag==FALSE){
        stop("Analysis has not been run. Perform using `$run()`")
      }
      return(copy(private$results))
    },

    #' @field mds
    #' Return the mds data (created when `<ReMIXTURE_object>$plot_MDS()` is run)
    mds=function(){
      if(is.null(private$mdsPlot)){
        stop("MDS has not been produced. Make it using e.g. `<ReMIXTURE_object>$plot_MDS()`")
      }
      return(copy(private$mdsPlot))
    }
  ),

  #### PUBLIC ####
  public = list(


    #### INTIALISER ####

    #' @description
    #' Create a new ReMIXTURE object.
    #' @param distance_matrix \[NULL\] An all-vs-all, full numeric distance matrix, with rownames and colnames giving the region of origin of the corresponding individual.
    #' @param region_table \[no default\] A data.table or data.frame describing the longitudes/latitudes of each region, with columns named "region" (character), and "lon" and "lat" (numeric or integer). The "region" column must have names corresponding to all the row/column names of the distance matrix.
    #' @param vcf_file \[NULL\] Path to a `.vcf` or `.vcf.gz` file. Provide this instead of `distance_matrix` to build the matrix via `SNPRelate`.
    #' @param sample_table \[NULL\] A data.frame or data.table with character columns named `sample` and `region`, mapping VCF sample IDs to ReMIXTURE regions.
    #' @return A new ReMIXTURE object.
    initialize = function(
      distance_matrix = NULL,
      region_table,
      vcf_file = NULL,
      sample_table = NULL
    ){


      ce("------------------------------------------------")
      ce("Initialising ReMixture object ...")
      ce("------------------------------------------------\n")

      if(!is.null(distance_matrix) && !is.null(vcf_file)){
        stop("Provide either `distance_matrix` or `vcf_file`, not both.")
      }
      if(is.null(distance_matrix) && is.null(vcf_file)){
        stop("Provide either `distance_matrix` or `vcf_file`.")
      }
      if(is.null(distance_matrix)){
        ce("\tConstructing distance matrix from VCF input ...")
        distance_matrix <- private$distance_matrix_from_vcf_inputs(
          vcf_file = vcf_file,
          sample_table = sample_table
        )
      }

      ce("\tModifying distance matrix if necessary ... ")
      bad_diagonal <- is.na(diag(distance_matrix)) | is.infinite(diag(distance_matrix))
      if(any(bad_diagonal)){
        ce("Diagonals of distance matrix contain `NA` or +/-Inf --- these will be replaced by zeroes.")
        diag_values <- diag(distance_matrix)
        diag_values[bad_diagonal] <- 0
        diag(distance_matrix) <- diag_values
      }

      ce("\tValidating input distance matrix ...")
      private$validate_m(distance_matrix)

      ce("\tValidating input region table ...")
      region_table <- private$validate_rt(region_table)

      ce("\tChecking distance matrix and region table compatibility ...")
      private$validate_dm_rt(distance_matrix,region_table)

      if(any(distance_matrix>1)){
        ce("Distance matrix has values > 1, and will now have all entries linearly scaled to fit the range [0,1]. If this is an issue, please provide a pre-scaled distance matrix.")
        distance_matrix <- distance_matrix %>% scale_between(0,1)
        stopifnot(all(diag(distance_matrix)==0))
      }

      ce("\tTrimming region table if necessary ... ")
      region_table <- region_table[region %in% colnames(distance_matrix)]

      ce("\tSaving input to object ...")
      private$m <- distance_matrix
      diag(private$m) <- Inf
      private$rt <- region_table

      ce("\tAdding sample counts to internal region table ...")
      if(!is.null(private$rt$N)){ warning("Column named 'N' discovered in region table. This will be overwritten. To preserve it, please rename it and initialise again.") }
      tmp <- as.data.table(table(colnames(private$m))) %>% setnames(c("V1"),c("region"))
      private$rt <- tmp[private$rt,on=.(region)]

      ce("\n------------------------------------------------")
      ce("Initialisation complete.")
      ce("------------------------------------------------")
    },

    #' @description
    #' Build a pairwise distance matrix directly from a `.vcf` or `.vcf.gz` file using `SNPRelate`.
    #' Direct `.bcf` input is not supported by this helper.
    #'
    #' The VCF is converted to a temporary GDS file with `snpgdsVCF2GDS()`, IBS is computed with `snpgdsIBS()`, and the distance matrix is returned as `1 - ibs`.
    #'
    #' @param vcf_file \[no default\] Path to a `.vcf` or `.vcf.gz` file.
    #' @param sample_table \[no default\] A data.frame or data.table with character columns named `sample` and `region`.
    #' @return A full square numeric distance matrix with region labels in the row and column names.
    distance_matrix_from_vcf = function(vcf_file, sample_table){
      private$distance_matrix_from_vcf_inputs(
        vcf_file = vcf_file,
        sample_table = sample_table
      )
    },

    #### RUN ####
    #' @description
    #' Run the ReMIXTURE algorithm and save the results in the object. Multiple runs can be requested to aid parameter selection (see description for `H`), which is strongly recommended.
    #'
    #' @param iterations \[1000\] The number of times subsamples are drawn, clustered, and counted.
    #' @param subsample_proportions \[0.8\] Size of subsample to draw each iteration. An equal number of samples will be selected from each region, the number being
    #'   \eqn{\text{subsample\_proportions} \times \text{\# of samples in region with the fewest samples}}
    #' (rounded to the nearest integer).
    #' @param H \["auto"\] Controls the clustering cut heights \eqn{H} used in each run. There are three modes possible:
    #' - (Default) Random \eqn{H}, truncated normal distributions. Each iteration, \eqn{H} will be drawn from a truncated normal distribution, truncated to the values given by `H_truncNorm_range`. To use this option, you must provide the means of the distributions you want to try in `H`, and the standard deviations (one corresponding to each mean) in `H_truncNorm_sd`.
    #'   - If `H="auto"` (the default), then a selection of 16 normal distributions will be chosen based on a simple heuristic that typically captures some very good values. Briefly, the values of `H` fully span the interguartile ranges of inter-sample distances of each region, `H_truncNorm_sd` is half the standard deviation of inter-sample distances in the region with the smallest such value (i.e. it is given the same value in all 16 distributions), and `H_truncNorm_range` is just the range of all inter-sample distances in the whole dataset.
    #' - Fixed \eqn{H}. To use this option, provide a list of values to try as \eqn{H}. A run will be conducted with each, using the same value at each iteration.
    #' - Random \eqn{H}, empirical distribution. This samples values of \eqn{H} from the distances in the distance matrix. This should probably never be used, because while it does help produce meaningful plots in cases of extremely structured populations, the results are just not very intuitive. If you use this mode, you should clearly indicate it in any figure or publication, e.g. "ReMIXTURE plot using empirically distributed \eqn{H}". But just don't. This option will be removed.
    #' @param H_truncNorm_sd \[NULL\] See description for `H`.
    #' @param H_truncNorm_range \[upper_tri_ply(private$m,range)\] See description for `H`. By default it is set to the range of values in the distance matrix.
    #' @param diagnosticPlotMDSclusters \[FALSE\] A very useful diagnostic tool, especially for weirdly-structured datasets. If true, at each iteration, will draw the subsampled samples on an MDS plot, and draw hulls around the clusters.
    #'
    #' @return Nothing
    run = function(
      iterations=1000,
      subsample_proportions=c(0.8),
      H="auto",
      H_truncNorm_sd=NULL,
      H_truncNorm_range=upper_tri_ply(private$m,range),
      diagnosticPlotMDSclusters=FALSE
    ){
      ce("------------------------------------------------")
      ce("Running ReMixture analysis ...")
      ce("------------------------------------------------\n")

      if( diagnosticPlotMDSclusters==TRUE & interactive()==FALSE ){
        ce("Diagnostic MDS clustering plots are only available in interactive sessions--disabling this feature.")
        diagnosticPlotMDSclusters <- FALSE
      }

      ce("\tResetting results fields and flags")
      private$results <- list()
      private$runflag <- FALSE
      private$plotflag <- FALSE

      ce("\tSetting up local variables and containers ...")
      #local params
      nits <- iterations
      ind_info <- data.table(
        gp=colnames(private$m),
        idx=1:ncol(private$m)
      )
      gp_vec <- colnames(private$m)
      gp_list <- sort(ind_info[,unique(gp)]) # check this lines up
      gp_info <- data.table(
        gp = sort(gp_list),
        gp_idx = 1:length(gp_list)
      )
      ind_info <- gp_info[ind_info,on="gp"]
      nind <- ncol(private$m)
      ngp <- nrow(gp_info)
      # param_test_insert <- 1:ngp
      results_insert <- 1

      H_type <- if(H[1]=="auto"){
        ce("\t\tCalculating some (hopefully) sensible parameters.")
        tmp <- apply(as.matrix(gp_list),1,function(gp){
          # gp<-"Africa"
          selRegionIdx <- which(colnames(private$m)==gp)
          utIdx <- upper.tri(private$m[selRegionIdx,selRegionIdx])
          matrix(
            c(
              quantile(private$m[selRegionIdx,selRegionIdx][utIdx],0.25),
              quantile(private$m[selRegionIdx,selRegionIdx][utIdx],0.75),
              # median(private$m[selRegionIdx,selRegionIdx][utIdx]),
              sd(private$m[selRegionIdx,selRegionIdx][utIdx])/2
            ),
            ncol=1
          )
        })

        h_min <- tmp[1,which(tmp[1,]==min(tmp[1,]))][1]
        h_max <- tmp[2,which(tmp[2,]==max(tmp[2,]))][1]
        H <- seq(h_min,h_max,length.out=16)
        H_truncNorm_sd <- min(tmp[3,which(tmp[3,]==min(tmp[3,]))][1])/2 %>% rep(16)

        ce("\tRandom H values following a selection of automatically chosen (truncated) normal distributions are requested.")
        ce("\t`H`=c(",paste0(signif(H,4),collapse=","),")")
        ce("\t`H_truncNorm_sd`=c(",paste0(signif(H_truncNorm_sd,4),collapse=","),")")
        "random-truncNormal"
      } else if(H[1]=="random-empirical"){
        "random-empirical"
      } else if (!is.null(H_truncNorm_sd)){
        ce("\tRandom H values following user-supplied (truncated) normal distributions are requested.")
        if(
          !(is.numeric(H) | is.integer(H)) |
          !(is.numeric(H_truncNorm_sd) | is.integer(H_truncNorm_sd)) |
          !(is.numeric(H_truncNorm_range) | is.integer(H_truncNorm_range)) |
          (length(H)!=length(H_truncNorm_sd)) |
          (length(H_truncNorm_range)!=2) |
          !(H_truncNorm_range[1]<H_truncNorm_range[2]) |
          any(H_truncNorm_sd<=0.0)
        ){
          stop("For random (truncated normal) H values, `H` and `H_truncNorm_sd` must be numerical or integers vectors of equal length, and `H_truncNorm_range` must be a length two vector specifying a range of numbers > 0.")
        }
        "random-truncNormal"
      } else {
        ce("\tConstant H requested.")
        if( !(is.numeric(H) | is.integer(H)) ){
          stop("H values must be numeric or integers.")
        }
        "constant"
      }


      # In case we are doing diagnostic MDS plots, get the coords
      if( diagnosticPlotMDSclusters==TRUE ){
        makePlot <- if(!is.null(private$mdsPlot)){
          if (all(private$mdsPlot$axes==c(1,2))){ FALSE } else { TRUE }
        } else { TRUE }

        if(makePlot==TRUE){
          colTable <- data.table(
            region = unique(colnames(private$m))
          )[,col:=colorRampPalette(c("#DD000088","#DDDD0088","#00DD0088","#0000DD88","#DD00DD88"))(.N)]
          mds <- cmdscale( as.dist(private$m) , k=2 )
          private$mdsPlot <- list(axes=c(1,2),mds=data.table(region=rownames(mds),axisA=mds[,1],axisB=mds[,2]),legend=colTable)
        }
      }

      for( pr_samp in subsample_proportions ){
        #dev pr_samp <- 0.8
        nSelectPerGrp <- round(ind_info[,.N,by=.(gp)][,min(N)] * pr_samp)
        if(nSelectPerGrp<2){
          private$results <- list()
          private$runflag <- FALSE
          private$plotflag <- FALSE
          stop(paste0("Some regions do not have enough samples for this value of `subsample_proportions` (",subsample_proportions,"), which results in a subsampling number of ",nSelectPerGrp," samples. Increasing `subsample_proportions` could help but most likely you should think about omitting low-sample regions or combining them into other regions."))
        }
        ce("Each iteration, ",nSelectPerGrp," samples from each region will be used.\n")
        for(i_hcut in seq_along(H)){
          #dev hcut="auto"; pr_samp=0.8
          sayHc <- if(H_type=="random-empirical"){
            H_type
          } else if (H_type=="random-truncNormal"){
            paste0("truncNorm( mean=",H[i_hcut],", sd=",H_truncNorm_sd[i_hcut]," ,range=[",H_truncNorm_range[1],",",H_truncNorm_range[2],"] )")
          } else if (H_type=="constant"){
            round(H[i_hcut],digits = 4)
          } else {
            stop("Something is fundamentally wrong with the universe, email mtrw85@gmail.com and alert Tim.")
          }
          ce( "Begin analysis for H==" , sayHc , " and subsample_proportions==" , pr_samp , " ..." )

          #container for parameter testing output
          # param_test_out <- expand.grid( #innermost loops first
          #   gp_idx=gp_info$gp, #
          #   it=1:nits,
          #   hcut=H,
          #   pr_samp=subsample_proportions,
          #   nclust=integer(1),
          #   run=integer(1)
          # ) %>% setDT()

          #local result containers
          nclust_counts <- rep(0L,length(gp_list))
          nclust_counts_2 <- rep(0L,length(gp_list))
          counts_mat <- matrix(0L,nrow=ngp,ncol=ngp)
          colnames(counts_mat) <- rownames(counts_mat) <- gp_list
          counts_mat_accumulator_empty <- counts_mat
          total_population_diversity <- 0.0
          counts_2_mat <- matrix(0.0,nrow=ngp,ncol=ngp)
          colnames(counts_2_mat) <- rownames(counts_2_mat) <- gp_list
          ind_info[,clust:=NA_integer_]

          skipPlotUntil <- 1L
          ce("\tIterating ...")
          for(it in 1:nits){ #Begin iteration loop
            if(it %% 100 == 0) {
              ce("\t\tBegin iteration: ",it)
            }
            #Subsample
            ss_selector <- ind_info[,.(s=sample(idx,nSelectPerGrp)),by=.(gp_idx)]$s

            # set hc
            hc <- if(H_type=="random-empirical"){
              sample(private$m[ss_selector,ss_selector],1)
            } else if (H_type=="random-truncNormal"){
              truncnorm::rtruncnorm(1,H_truncNorm_range[1],H_truncNorm_range[2],H[i_hcut],H_truncNorm_sd[i_hcut])
            } else if (H_type=="constant"){
              H[i_hcut]
            } else {
              stop("Something is fundamentally wrong with the universe, email mtrw85@gmail.com and alert Tim.")
            }

            #Cluster
              #reset
            ind_info[,clust:=NULL]
            ind_info[ss_selector, clust:=cutree(hclust(as.dist(private$m[ss_selector,ss_selector])),h=hc)]


            # Diagnostic plots if requested
            if( diagnosticPlotMDSclusters==TRUE & skipPlotUntil==it)
            {
              plot(
                private$mdsPlot$mds[ss_selector,]$axisA,
                private$mdsPlot$mds[ss_selector,]$axisB,
                col=private$mdsPlot$legend[data.table(region=private$mdsPlot$mds[ss_selector,]$region),on=.(region)]$col,
                xlim=range(private$mdsPlot$mds$axisA),
                ylim=range(private$mdsPlot$mds$axisB),
                xlab="Axis 1",
                ylab="Axis 2",
                pch=20,
                cex=0.4,
                main=paste0( "Iteration " , it , "; Subsampling " , pr_samp*100,"%; H (",H_type,"): ", round(hc,2) )
              )

              for(cl in unique(ind_info[ss_selector,]$clust)){
                cl_selector <- which(ind_info$clust==cl)
                if(length(cl_selector)==1){
                  points(
                    private$mdsPlot$mds[cl_selector,]$axisA,
                    private$mdsPlot$mds[cl_selector,]$axisB,
                    col="#00000077",
                    pch=20,
                    cex=2
                  )
                } else if(length(cl_selector)==2){
                  lines(
                    private$mdsPlot$mds[cl_selector,]$axisA,
                    private$mdsPlot$mds[cl_selector,]$axisB,
                    col="#00000077",
                    lwd=11
                  )
                } else {
                  hull <- chull(private$mdsPlot$mds[cl_selector,]$axisA, private$mdsPlot$mds[cl_selector,]$axisB)
                  polygon(
                    private$mdsPlot$mds[cl_selector,][hull,]$axisA,
                    private$mdsPlot$mds[cl_selector,][hull,]$axisB,
                    col="#00000033"
                  )
                }
              }
              skipPlotUntil <- skipPlotUntil+1L
              ans <- ask(paste0("[Round ",it,"] Press <return> for the next round, enter an integer N to skip to the Nth round, or 'f' to stop showing clusters and finish all rounds."),YN=FALSE)
              if(!is.na(as.integer(ans)%>%suppressWarnings())){ skipPlotUntil<-(as.integer(ans)%>%suppressWarnings()); ce("You have requested another plot at round ",skipPlotUntil,"; You are currently at round ",it,".")}
              if(ans=="f"){ skipPlotUntil <- -1L; ce("No more plots will be drawn.") }
            }

            counts_mat_accumulator <- counts_mat_accumulator_empty
            nclust_counts <- nclust_counts + (t<-ind_info[ss_selector,.(add=nu(clust)),by=.(gp_idx)][gp_info,on=.(gp_idx)][is.na(add),add:=0L][]$add) # in how many unique clusters does each region occur
            nclust_counts_2 <- nclust_counts_2 + t**2

            # param_test_out[param_test_insert]$run <- results_insert
            # param_test_insert <- param_test_insert + ngp

            ind_info[ss_selector,{ #over clusters
              ugidx <- unique(gp_idx)
              if(length(ugidx) == 1){
                counts_mat_accumulator[gp_idx,gp_idx] <<- counts_mat_accumulator[gp_idx,gp_idx] + 1 #solo cluster--contributes to genetic 'uniqueness'
              } else {
                apply(combn(ugidx,2),2,function(c) { counts_mat_accumulator[c[1],c[2]] <<- counts_mat_accumulator[c[1],c[2]] + 1 } ) #over permutations of members of multigroup cluster
              }
            },by=.(clust)] %>% invisible

            total_population_diversity <- total_population_diversity + ind_info[ss_selector,nu(clust)]

            counts_mat_accumulator <- fold_matrix(counts_mat_accumulator)
            counts_mat <- counts_mat + counts_mat_accumulator
            counts_2_mat <- counts_2_mat + counts_mat_accumulator**2

          } #end iteration loop
          #browser()
          ce("\tSummarising and saving results ...")
          private$results[[results_insert]] <- list(
            subsample_proportion=pr_samp,
            H_type=H_type,
            H=if(H_type %in% c("constant","random-truncNormal") ){H[i_hcut]}else{NULL},
            H_truncNorm_sd=if(H_type %in% c("random-truncNormal") ){H_truncNorm_sd[i_hcut]}else{NULL},
            H_truncNorm_range=if(H_type %in% c("random-truncNormal") ){H_truncNorm_range}else{NULL},
            iterations=nits,
            total_population_diversity = total_population_diversity/nits,
            overlap = counts_mat/nits,
            var_overlap = (counts_2_mat/nits) - (counts_mat/nits)**2,
            diversity = nclust_counts/nits,
            var_diversity = (nclust_counts_2/nits) - (nclust_counts/nits)**2
          )

          results_insert <- results_insert + 1
        } # end hcut loop
      } # end pr_samp loop

      # #Save param test output
      # private$results$parameter_selection_clustercounts <- param_test_out

      private$runflag <- TRUE
      ce("\n------------------------------------------------")
      ce("ReMixture analysis complete ...")
      ce("------------------------------------------------")
    },


    #### PLOTTING ####

    #' @description
    #' Plot the unique and overlapped diversity recorded in a ReMIXTURE run, in something like a heatmap format but using the circle-area conventions as per `plot_maps()`.
    #'
    #' @param runs \[NULL\] For which runs would you like a plot? An integer vector. By default, all runs.
    #' @param ... Additional arguments ultimately passed to `base::plot()`
    #'
    #' @return Nothing
    plot_results_grid = function( runs=NULL ,...){
      if(private$runflag==FALSE){
        stop("Analysis has not been run. Perform using `$run()`")
      }
      if(is.null(runs)){
        ce("Plotting a results grids for all runs.")
        runs <- 1:length(private$results)
      }
      for(i in runs){
        ssp <- private$results[[i]]$subsample_proportion
        hc <- private$results[[i]]$H
        hcsd <- private$results[[i]]$H_truncNorm_sd
        its <- private$results[[i]]$iterations

        ce(paste0("Plotting run ",i,":\n\tSubsample proportion: ",ssp,"\n\tH mean: ",round(hc,digits = 2),"\n\tH sd:",hcsd,"\n\tIterations:",its,"\n------------------------------\n"))

        td <- private$results[[i]]$diversity
        ol_ud <- private$results[[i]]$overlap
        maxDiv <- max(c(td,ol_ud))
        tdScale <- td/maxDiv
        ol_udScale <- ol_ud/maxDiv
        n <- m <- nrow(ol_ud)

        #null_plot(1:(n+1),1:(m+1),xaxt="n",yaxt="n", ylab="This focal region ...",xlab="... is overlapped by this region",main="Results grid\nTotal/unique diversity (diagonal) &\noverlapped diversity (off-diagonal)")
        null_plot(1:(n+1),1:(m+1),xaxt="n",yaxt="n", ylab="This focal region ...",xlab="... is overlapped by this region",main=paste0("Results grid for run ",i,"\nTotal/unique diversity (diagonal) &\noverlapped diversity (off-diagonal)"))#,...)
        axis(2,at=(1:n)+0.5,labels=rownames(ol_ud),las=2,cex.axis=0.6)
        axis(1,at=(1:n)+0.5,labels=colnames(ol_ud),las=3,cex.axis=0.6)
        mtext(paste0("(Maximum circle size=",maxDiv," clusters)"), side = 4,cex=.6)
        for(i in 1:n){
          for(j in 1:m){
            drawUnitSquareTopRight(i,j,col="#00000000",lwd=0.2)
            if(i==j){
              drawUnitCircleTopRight(i,j,scale=tdScale[i]     ,col="#000000")
              drawUnitCircleTopRight(i,j,scale=ol_udScale[i,j],col="#FFFFFF")
            } else {
              drawUnitCircleTopRight(i,j,scale=ol_udScale[i,j],col="#000000")
            }
          }
        }
      }
    },


    #' @description
    #' Plot the cluster counts for each region in each run. This is useful for assessing how different values of \eqn{H} affect the output, in particular checking there aren't regions with very low cluster counts, or that they are not "maxxed out", either of which will probably give dodgy results unreliable.
    #'
    #' @return Nothing
    plot_clustercounts = function(){
      if(private$runflag==FALSE){
        stop("Analysis has not been run. Perform using `$run()`")
      }

      pd <- setDT(ldply(1:length(private$results),function(i){
        data.table(
          run=i,
          region=private$results[[i]]$overlap %>% rownames,
          cluster_count=private$results[[i]]$diversity
        )
      }))

      setkey(pd,run,region)
      pd[,xIdx:=1:.N,by=run]

      #i=1
      null_plot(pd$xIdx,pd$cluster_count,xaxt="n",ylab="Cluster count (run #)")
      axis(1,1:nu(pd$region),pd[run==1,region],las=2,cex.axis=0.6)
      for(i in unique(pd$run)){
        lines(pd[run==i]$xIdx,pd[run==i]$cluster_count,type="b",pch=20)
      }
      pd[region==region[1],text(xIdx,cluster_count,paste0("(",run,")"),cex=0.5,pos=2)]

      invisible(NULL)
    },

    #' @description
    #' These counts are excellent for determining a good value for \eqn{H}, in cases where the default choice is not desirable. A good value typically gives a nice balance of multi- and single-region clusters, typically between where the two values cross over and (as occurs in most datasets) a point where the number of multi-region clusters hits some maximum.
    #'
    #' @return Nothing
    plot_h_optimisation = function(){
      if(private$runflag==FALSE){
        stop("Analysis has not been run. Perform using `$run()`")
      }
      d <- setDT(ldply(private$results,function(r){
        data.table(
          subsample_proportion = r$subsample_proportion,
          H = r$H,
          H_type = r$H_type,
          iterations = r$iterations,
          clustercount_uniq = weighted.mean(diag(r$overlap),w=1/private$rt$N),
          clustercount_shared = weighted.mean(r$diversity-diag(r$overlap),w=1/private$rt$N) #or lower tri
        )
      }))[,run:=1:.N][]

      d[,{
        null_plot(run,c(clustercount_uniq,clustercount_shared),xaxt="n",ylab="Aggregated cluster count",xlab=paste0("Run\n(H) [mean; type=\"",H_type[1],"\"]"),)
        abline(v=run,lty=2,col="#00000044")
        axis(1,run,paste0(run,"\n(",signif(H,3),")."),padj = -0.2,cex.axis=0.6)
        lines(run,clustercount_uniq,type="b",col="#880000",pch=20)
        lines(run,clustercount_shared,type="b",col="#ffcb00",pch=20)
        legend(1,max(c(clustercount_uniq,clustercount_shared)),c("Single-region clusters","Multi-region clusters"),fill=c("#880000","#ffcb00"))
      }]
      invisible(NULL)
    },

    #' @description
    #' Produces the canonical ReMIXTURE plots.
    #'
    #' @param run \[NULL\] If multiple runs were done, choose which run to use. This is best assessed using `<ReMIXTURE Object>$plot_h_optimisation()`, `<ReMIXTURE Object>$plot_results_grid()`, and `<ReMIXTURE Object>$plot_clustercounts()`.
    #' @param focalRegion \[NULL\] A character string naming one of the regions as it occurs in the region table, chosen as the focal region (from which lines showing inter-region overlap will emanate). If left NULL then a map will be plotted with each region as the focus in turn--this is best used after setting up a multi-panel plot with some variant of `par(mfrow=c(<number of rows>,<number of columns>))`.
    #' @param range_lon \[c(-179.0,179.0)\] Limit the map to some range of longitudes. Must be a vector of 2 numbers. For awkward reasons, using 180 or -180 can cause graphical bugs and I am deeply sorry.
    #' @param range_lat \[c(-85.0,85.0)\] Limit the map to some range of latitudes. Must be a vector of 2 numbers.
    #' @param width_max \[10.0\] The maximum width of the circles/lines, in units of lat/lons. This width will correspond to the highest cluster count and the others will be scaled such that zero clusters <=> zero width.
    #' @param alpha_max \[1.0\] As per width_max but controlling the alpha of connecting lines. Setting to NULL will disable alpha and make all lines solid. Can be set above 1.0, with weird results--probably don't do this.
    #' @param diversityCirclesFocalOnly \[FALSE\] Plot the circle representing a region's total/unique/shared diversity at only the focal region. Otherwise, all plots will show the diversities at all regions. This is good when you have multiple plots on the page, see description for `focalRegion`.
    #' @param projection \[EckertIV\] A function the performs projection of the map coordinates. Included in the package are `eckertIV`, `winkelIII`, and `equirectangular`. A custom function can be given--see details.
    #' @param curvature_matrix \[NULL\] A matrix describing how the lines emanating from the focal region should bend. The (i,j)th entry gives the bend angle (in radians) of the line beginning at region i and ending at region j (indexed as per the order in the region table). In practice this is pretty tedious to enter manually. The argument "random" will auto-generate a matrix filled with entries ~ normal(0,0.3), which usually works well.
    #' @param mapData \[mapData110\] The polygon data from which the map is made. For higher resolution, use mapData10. See e.g. ?mapData110.
    #'
    #' @details
    #' A custom projection function should follow these rules: its first argument, named 'dtLL' should take a data.table with columns 'lon' and 'lat'. The second argument, named 'projColNames', should take a length-2 character vector. The output should be a data.table with columns containing the transformed longitude and latitude values. These columns should be named as per the entries in the input 'projColNames'.
    #'
    #'
    #' @return Nothing
    plot_maps = function(
      run=NULL,
      focalRegion=NULL,
      range_lon=c(-179.0,179.0),
      range_lat=c(-85.0,85.0),
      width_max=10.0,
      alpha_max=1.0,
      diversityCirclesFocalOnly=FALSE,
      projection=eckertIV,
      curvature_matrix=NULL,
      mapData=mapData110
    ){
      state <- private$get_map_plot_state(run, width_max, alpha_max, curvature_matrix)
      focalRegions <- if(is.null(focalRegion)){
        state$regions
      } else {
        state$regions[state$regions %in% focalRegion]
      }
      focalRegions <- private$validate_focal_regions(focalRegions, state$regions)

      for(region_name in focalRegions){
        private$render_map_panel(
          state = state,
          focalRegion = region_name,
          range_lon = range_lon,
          range_lat = range_lat,
          projection = projection,
          mapData = mapData,
          overview = FALSE,
          diversityCirclesFocalOnly = diversityCirclesFocalOnly,
          returnAnchors = FALSE
        )
      }
    },

    #' @description
    #' Plot a composite overview map with inset focal-region maps.
    #'
    #' @param run \[NULL\] If multiple runs were done, choose which run to use.
    #' @param focalRegions \[NULL\] Character vector of focal regions to display. If NULL, regions whose overview anchors are visible inside the current central overview crop are used. If provided, the supplied regions are validated first; regions outside the current overview crop may be dropped according to the current validation path.
    #' @param range_lon \[c(-179.0,179.0)\] Limit the map to some range of longitudes.
    #' @param range_lat \[c(-85.0,85.0)\] Limit the map to some range of latitudes.
    #' @param width_max \[10.0\] Maximum width of circles/lines.
    #' @param alpha_max \[1.0\] Maximum alpha scaling for overlap curves.
    #' @param projection \[EckertIV\] Projection function for the map coordinates.
    #' @param curvature_matrix \[NULL\] Optional curvature matrix or `"random"`.
    #' @param mapData \[mapData110\] Polygon data from which the map is made.
    #' @param panels \[NULL\] Optional data.table/data.frame for manual reuse or editing of inset/card layout. Legacy `left/right/bottom/top` input is accepted and interpreted as the map rectangle. The full schema supports `region`, `side`, `card_left`, `card_right`, `card_bottom`, `card_top`, `map_left`, `map_right`, `map_bottom`, `map_top`, `label_x`, `label_y`, `label_adj_x`, `label_adj_y`, plus `left/right/bottom/top` aliases for `map_*`.
    #' @param layout \[`"slot"`\] Layout engine for automatic inset placement. `"slot"` is the default stable layout. `"ellipse"` is an experimental squashed-circle layout that places evenly spaced inset panels around the central map and uses overview-anchor angles to assign regions to those positions.
    #' @param main_panel \[c(left = 0.28, right = 0.72, bottom = 0.30, top = 0.68)\] Normalized device coordinates for the central overview panel.
    #' @param inset_width \[0.18\] Width of each inset panel in normalized device coordinates.
    #' @param inset_height \[NULL\] Height of each inset panel in normalized device coordinates. If NULL, it will be estimated from the current map crop and device aspect ratio. Very small explicit values can visually flatten inset maps on wide devices; publication figures can set this manually or edit the returned `panels`.
    #' @param panel_aspect \[NULL\] Optional target width/height aspect ratio for inset panels. If NULL, a Tim-style bounded aspect based on the crop is used.
    #' @param leader_lines \[TRUE\] If TRUE, draw dashed leader lines between inset panels and the overview anchors.
    #' @param panel_labels \[TRUE\] If TRUE, draw panel labels outside the inset maps.
    #' @param panel_label_cex \[0.75\] Text size for inset panel labels.
    #' @param overview_circle_scale \[0.35\] Scale factor applied to diversity circles in the central overview map. This is a composite-only visual scaling parameter and does not affect analysis results or `plot_maps()`.
    #' @param inset_circle_scale \[0.55\] Scale factor applied to diversity circles in inset maps. This is a composite-only visual scaling parameter and does not affect analysis results or `plot_maps()`.
    #' @param inset_line_scale \[0.75\] Scale factor applied to overlap-curve widths in inset maps. This is a composite-only visual scaling parameter and does not affect analysis results or `plot_maps()`.
    #' @param inset_style \[`"default"`\] Controls visual emphasis for inset panels. `"default"` uses the standard map rendering. `"line_emphasis"` reduces the visual dominance of inset circles and emphasises overlap curves; this is useful when inset panels are small.
    #'
    #' @details
    #' The method draws a slot-based composite figure around a central overview map. Large `inset_width` values require sufficient horizontal space around `main_panel`. For auto-generated layouts, impossible side geometry fails with a clear error instead of drawing an overlapping figure. For user-supplied `panels`, potentially overlapping layouts are warned about rather than automatically modified.
    #' When `focalRegions = NULL`, all visible overview anchors are used. If many regions are visible, the automatic layout may fail clearly rather than drawing overlapping panels; users can pass fewer `focalRegions`, reduce `inset_width`, adjust `main_panel`, or edit the returned `panels`.
    #'
    #' Before using this method for final interpretation, users are encouraged to inspect H behaviour with the existing ReMIXTURE diagnostic tools, such as `plot_h_optimisation()`, `plot_results_grid()`, and `plot_distance_densities()`, then pass a selected run to `plot_maps_composite()`.
    #' A typical workflow is `rm$run(...)`, `rm$plot_h_optimisation()`, `rm$plot_results_grid()`, `rm$plot_distance_densities(HdistFromRun = selected_run)`, and finally `rm$plot_maps_composite(run = selected_run)`.
    #'
    #' The overview map keeps the full circle-based visual encoding. The `"line_emphasis"` inset style is intended for publication-style composite figures where inset maps should highlight overlap curves more than node circles.
    #'
    #' When `panels` is left NULL, `layout = "slot"` uses the default slot-based layout and `layout = "ellipse"` uses an experimental ellipse layout. The returned panels can be edited and passed back manually through `panels`.
    #'
    #' The returned layout uses a full card/map/label schema so it can be edited and passed back through `panels` for publication-style manual refinement.
    #'
    #' @return Invisibly returns a `data.table` describing the composite panel layout. The returned table can be edited and passed back through `panels` to manually fine-tune publication figures. It includes at least `region`, `side`, `card_left`, `card_right`, `card_bottom`, `card_top`, `map_left`, `map_right`, `map_bottom`, `map_top`, `label_x`, `label_y`, `label_adj_x`, `label_adj_y`, plus `left/right/bottom/top` aliases equal to `map_*`.
    plot_maps_composite = function(
      run = NULL,
      focalRegions = NULL,
      range_lon = c(-179.0, 179.0),
      range_lat = c(-85.0, 85.0),
      width_max = 10.0,
      alpha_max = 1.0,
      projection = eckertIV,
      curvature_matrix = NULL,
      mapData = mapData110,
      panels = NULL,
      layout = c("slot", "ellipse"),
      main_panel = c(left = 0.28, right = 0.72, bottom = 0.30, top = 0.68),
      inset_width = 0.18,
      inset_height = NULL,
      panel_aspect = NULL,
      leader_lines = TRUE,
      panel_labels = TRUE,
      panel_label_cex = 0.75,
      overview_circle_scale = 0.35,
      inset_circle_scale = 0.55,
      inset_line_scale = 0.75,
      inset_style = c("default", "line_emphasis")
    ){
      state <- private$get_map_plot_state(run, width_max, alpha_max, curvature_matrix)
      main_panel <- private$normalize_main_panel(main_panel)
      layout <- match.arg(layout)
      inset_style <- match.arg(inset_style)

      op <- par(no.readonly = TRUE)
      on.exit(par(op), add = TRUE)

      plot.new()
      par(fig = unname(main_panel[c("left", "right", "bottom", "top")]), mar = c(0.2, 0.2, 1, 0.2), new = TRUE, xpd = FALSE)
      anchors <- private$render_map_panel(
        state = state,
        focalRegion = focalRegions[1],
        range_lon = range_lon,
        range_lat = range_lat,
        projection = projection,
        mapData = mapData,
        overview = TRUE,
        diversityCirclesFocalOnly = FALSE,
        returnAnchors = TRUE,
        circle_scale = overview_circle_scale,
        clip_to_panel = TRUE,
        border_last = TRUE,
        tight_axes = TRUE,
        draw_grid = FALSE
      )
      visible_anchors <- private$filter_visible_anchors(anchors, main_panel)

      if(!is.null(panels) && is.null(focalRegions)){
        if(!(is.data.frame(panels) || data.table::is.data.table(panels)) || !"region" %in% colnames(panels)){
          stop("When `panels` is supplied and `focalRegions` is NULL, `panels` must contain a `region` column.")
        }
        focalRegions <- as.character(panels$region)
      }

      focalRegions <- if(is.null(focalRegions)){
        visible_anchors$region
      } else {
        focalRegions <- private$validate_focal_regions(focalRegions, state$regions)
        dropped_regions <- setdiff(focalRegions, visible_anchors$region)
        if(length(dropped_regions) > 0){
          warning(
            "Dropping focal regions outside the current overview crop: ",
            paste(dropped_regions, collapse = ", ")
          )
        }
        focalRegions[focalRegions %in% visible_anchors$region]
      }
      if(length(focalRegions) == 0){
        stop("No focal regions are visible within the current `range_lon` / `range_lat`. Expand the crop or explicitly choose visible regions.")
      }

      visible_anchors <- visible_anchors[match(focalRegions, region)]
      inset_height <- private$compute_inset_height(
        inset_width = inset_width,
        inset_height = inset_height,
        range_lon = range_lon,
        range_lat = range_lat,
        panel_aspect = panel_aspect
      )

      panels <- if(is.null(panels)){
        if(layout == "slot"){
          private$make_composite_layout(visible_anchors, focalRegions, main_panel, inset_width, inset_height)
        } else {
          private$make_ellipse_layout(visible_anchors, focalRegions, main_panel, inset_width, inset_height)
        }
      } else {
        private$validate_composite_layout(panels, focalRegions, main_panel)
      }

      leader_segments <- private$make_leader_segments(panels, visible_anchors)
      if(leader_lines){
        private$render_leader_segments(leader_segments)
      }

      panel_regions <- panels$region
      for(i in seq_len(nrow(panels))){
        par(fig = c(panels$map_left[i], panels$map_right[i], panels$map_bottom[i], panels$map_top[i]), mar = c(0.2, 0.2, 1, 0.2), new = TRUE, xpd = FALSE)
        private$render_map_panel(
          state = state,
          focalRegion = panels$region[i],
          range_lon = range_lon,
          range_lat = range_lat,
          projection = projection,
          mapData = mapData,
          overview = FALSE,
          diversityCirclesFocalOnly = TRUE,
          returnAnchors = FALSE,
          circle_scale = inset_circle_scale,
          line_scale = inset_line_scale,
          targetRegions = panel_regions,
          clip_to_panel = TRUE,
          border_last = TRUE,
          tight_axes = TRUE,
          draw_grid = FALSE,
          draw_title = FALSE,
          inset_style = inset_style
        )
      }

      if(panel_labels){
        private$render_panel_labels(panels, cex = panel_label_cex)
      }

      invisible(panels)
    },
    #' @description
    #' Plot inter-sample distances as a density plot over the distance matrix--basically, a convenient way to judge how far apart samples tend to be, and thus how the clustering might behave when various H values (\eqn{H}) are used.
    #'
    #' Each plot tells the story from the perspective of a focal region. Each region is assigned a colour that remains constant across plots. A region's samples' distances to each other are marked with a thick line allowing you to see what colour is for what region.
    #'
    #' @param set_bw \[NULL\] The bandwidth parameter 'bw' passed to `density()`. By default it will be set automatically.
    #' @param set_xlims \[NULL\] Limits on the x-axis. Default is set by `density()`.
    #' @param samePlot \[FALSE\] If TRUE, it will put all the plots on one page.
    #' @param H \[NULL\] Visualise chosen values of (or distribution over) \eqn{H} on the plot. If only this argument is used, a single line will indicate that H value.
    #' @param H_truncNorm_sd \[NULL\] Visualise a given distribution of \eqn{H} on the plot. If this argument is used (in conjunction with `H`), a truncated normal distribution will be shown.
    #' @param H_truncNorm_lims \[upper_tri_ply(private$m,range)\] Visualise your values of \eqn{H} on the plot. If this argument is used, it will set the truncation limits of the normal distribution. By default it will be the range of inter-sample distances in the distance matrix.
    #' @param HdistFromRun \[NULL\] Show the H distribution for a particular run. Obviously, this requires that a run has been done.
    #' @param plotLegend \[FALSE\] If true, the last plot produced will be the colour legend.
    #'
    #' @return Invisibly returns a data.table giving the legend colours.
    plot_distance_densities = function(
      set_bw=0.001,
      set_xlims=upper_tri_ply(private$m,range),
      samePlot=FALSE,
      H=NULL,
      H_truncNorm_sd=NULL,
      H_truncNorm_lims=NULL,
      HdistFromRun=NULL,
      plotLegend=FALSE
    ){

      if(!is.null(HdistFromRun)){
        if(private$runflag==FALSE){
          stop("Analysis has not been run. Perform using `$run()`")
        }
        if(HdistFromRun==TRUE){ HdistFromRun <- 1L }
        ce("Plotting H distribution from run ",HdistFromRun)
        H <- private$results[[HdistFromRun]]$H
        H_truncNorm_sd <- private$results[[HdistFromRun]]$H_truncNorm_sd
        H_truncNorm_range <- private$results[[HdistFromRun]]$H_truncNorm_range
      }

      dt1 <- ldply(unique(colnames(private$m)),function(r){ #dev r = "Africa"
        selr <- rownames(private$m)==r
        ldply(unique(colnames(private$m)),function(c){ #dev r = "Africa"
          selc <- colnames(private$m)==c
          data.table(
            x=(0:5000)/5000,
            y=density(private$m[selr,selc],bw=set_bw,from=0,to=1,n=5001)$y,
            region1=c,
            region2=r
          )
        })
      }) %>% setDT

      colTable <- data.table(
        region2 = unique(colnames(private$m)),
        col     = rgb( t(col2rgb(hsv(seq(0, 0.8, length.out = nu(colnames(private$m))), 1, 1)) / 255) )
      )
      dt1 <- colTable[dt1,on=.(region2)]
      dt1[,col:=rgb(t(col2rgb(col))/255,alpha=fifelse(region1==region2,1,0.4))]
      dt1[,lwd:=fifelse(region1==region2,3,1)]

      if(samePlot==TRUE){
        warning("Altering mfrow parameters for multi-plot graphics. It will be reset to c(1,1).")
        par( mfrow=c(ceiling(sqrt(nrow(private$rt))),ceiling(sqrt(nrow(private$rt)))) )
      }

      rangeY <- range(dt1$y)
      for(r1 in unique(colnames(private$m))){
        sel1 <- rownames(private$m)==r1
        null_plot(set_xlims,rangeY,main=paste0(r1))
        for(r2 in unique(colnames(private$m))){
          sel2 <- rownames(private$m)==r2
          dt1[ region1==r1 & region2==r2 , lines(x,y,col=col,lwd=lwd) ]
        }
        if(!is.null(H)){
          if(!is.null(H_truncNorm_sd)){
            if(is.null(H_truncNorm_lims)){
              ce("Argument `H_truncNorm_lims` not provided. Will be set to +/- Inf.")
              H_truncNorm_lims <- c(-Inf,Inf)
            }
            require(truncnorm)
            polygon(
              x=c(
                seq(set_xlims[1],set_xlims[2],l=1000L),
                set_xlims[2],
                set_xlims[1]
              ),
              y=c(
                ((t<-truncnorm::dtruncnorm(seq(set_xlims[1],set_xlims[2],l=1000L),mean=H,sd=H_truncNorm_sd,a=H_truncNorm_lims[1],b=H_truncNorm_lims[2]))/max(t)) * rangeY[2],
                0,
                0
              ),
              border="#00000000",
              col="#00000022"
            )
          } else {
            abline(v=H,lty=2,col="#00000033")
          }
        }
      }
      if(samePlot==TRUE){
        par( mfrow=c(1,1) )
      }

      if(plotLegend==TRUE){
        null_plot(0:2,0:nrow(colTable),xaxt="n",yaxt="n")
        for(i in nrow(colTable):1){
          rect(0,i-1,1,i,col=colTable[i]$col)
          text(x=1.0,y=i-0.5,label=colTable[i]$region2,pos=4,cex=0.8)
        }
      }
      invisible(return(colTable))
    },
    #' @description
    #' An MDS plot showing relationships among samples
    #'
    #' An MDS ordination plot (in the same flavour as a PCA--in fact in some circumstances they will be the same depending on what your distances represent), useful as a first-pass assessment of the relationships among your samples--what overlaps what, what is over/underrepresented in the sample set, etc.
    #'
    #' These can take a very long time if the number of samples is high. In these cases, I suggest taking a smaller sample of the whole dataset.
    #'
    #' @param axes \[c(1,2)\] A vector of the MDS axes that should be outputted.
    #' @param colPalette \[c("#DD000088","#DDDD0088","#00DD0088","#0000DD88","#DD00DD88")\] A vector collection of colours that will be interpolated and applied to the different
    #' @param showLegend \[TRUE\] Add a legend?
    #' @param doPlot \[TRUE\] Produce a plot? If not, will still return the data for the MDS.
    #' @param pch \[20\] The shape of the points. See `?points`.
    #' @param cex \[0.8\] The size of the points. See `?points`.
    #' @param ... Other arguments passed to `plot()`.
    #'
    #' @return Invisibly returns a 2-list containing the (1) results of the MDS in a data.table and (2) information about the legend (another data.table).
    plot_MDS = function(
      axes=c(1L,2L),
      colPalette=c("#DD000088","#DDDD0088","#00DD0088","#0000DD88","#DD00DD88"),
      showLegend=TRUE,
      doPlot=TRUE,
      pch=20,
      cex=0.8,
      ...
    ){
      if( length(axes)!=2L | !(is.numeric(axes) | is.integer(axes)) ){ stop("`axes` must be a numeric or integer vector of length 2") }
      colTable <- data.table(
        region = unique(colnames(private$m))
      )[,col:=colorRampPalette(colPalette)(.N)]
      dim <- max(axes)

      makePlot <- if(!is.null(private$mdsPlot)){
        if (all(private$mdsPlot$axes==axes)){ FALSE } else { TRUE }
      } else { TRUE }

      if(makePlot==TRUE){
        mds <- cmdscale( as.dist(private$m) , k=dim )  # Also has function of ignoring the Inf diagonals (done for reasons to make the main algorithm work)
        private$mdsPlot <- list(axes=axes,mds=data.table(region=rownames(mds),axisA=mds[,axes[1]],axisB=mds[,axes[2]]),legend=colTable)
      }

      if(doPlot==TRUE){
        plot(
          private$mdsPlot$mds$axisA,
          private$mdsPlot$mds$axisB,
          pch=pch,
          col=colTable[data.table(region=private$mdsPlot$mds$region),on=.(region)]$col,
          cex=cex,
          xlab=paste0("Axis ",axes[1]),
          ylab=paste0("Axis ",axes[2]),
          ...
        )
        if(showLegend==TRUE){
          legend(min(private$mdsPlot$mds$axisA),max(private$mdsPlot$mds$axisB),colTable$region,colTable$col,bg="#FFFFFFAA")
        }
      }

      invisible(return(private$mdsPlot))
    }



  )
)

