.GRcalculate = function(inputData, groupingColumns) {
  log2nn = with(inputData, log2(cell_count/cell_count__time0))
  log2nn_ctrl = with(inputData, log2(cell_count__ctrl/cell_count__time0))
  GR = 2^(log2nn/log2nn_ctrl) - 1
  input_edited = inputData
  input_edited$log10_concentration = log10(input_edited$concentration)
  input_edited$GR = GR
  tmp<-input_edited[,groupingColumns, drop = F]
  experimentNew = (apply(tmp,1, function(x) (paste(x,collapse=" "))))
  if(length(groupingColumns) > 0) {
    input_edited$experiment = as.factor(experimentNew)
  } else {
    input_edited$experiment = as.factor("All Data")
  }
  return(input_edited)
}

.GRlogisticFit = function(inputData, groupingColumns) {
  experiments = levels(inputData$experiment)
  parameters = matrix(data = NA, ncol = 3, nrow = length(experiments))
  parameters = as.data.frame(parameters)
  colnames(parameters) = c('Hill','GRinf','GEC50')
  if(length(groupingColumns) > 0) {
    metadata = matrix(data = NA, ncol = length(groupingColumns), nrow = length(experiments))
    metadata = as.data.frame(metadata)
    colnames(metadata) = groupingColumns
  } else {
    metadata = NULL
  }
  pval = NULL
  GRmax = NULL
  GR_mean = NULL
  AOC = NULL
  R_square = NULL
  for(i in 1:length(experiments)) {
    # print(i)
    data_exp = inputData[inputData$experiment == experiments[i],]
    concs = sort(unique(data_exp$concentration))
    l = length(concs)
    max_concs = data_exp[data_exp$concentration %in% concs[c(l,l-1)],]
    GRmax[i] = min(max_concs$GR)
    #     metadata[i,] = data_exp[1,1:5]
    if(!is.null(metadata)) {
      metadata[i,] = data_exp[1,groupingColumns, drop = F]
    }
    GR_mean[i] = mean(data_exp$GR)
    #===== constrained fit ============
    c = unique(data_exp$concentration)
    priors = c(2, 0.1, median(c))
    lower = c(.1, -1, min(c)*1e-2)
    upper = c(5, 1, max(c)*1e2)
    if(dim(data_exp)[1] > 1) {
      output_model_new = try(drc::drm(GR~concentration, experiment, data=data_exp, fct=LL.3u(names=c('Hill','GRinf','GEC50')), start = priors, lowerl = lower, upperl = upper))
      if(class(output_model_new)!="try-error") {
        parameters[i,] = c(as.numeric(coef(output_model_new)))
        # F-test for the significance of the sigmoidal fit
        Npara = 3 # N of parameters in the growth curve
        Npara_flat = 1 # F-test for the models
        RSS2 = sum(residuals(output_model_new)^2)
        RSS1 = sum((data_exp$GR - mean(data_exp$GR))^2)
        df1 = (Npara - Npara_flat)
        df2 = (length(data_exp$GR) - Npara + 1)
        f_value = ((RSS1-RSS2)/df1)/(RSS2/df2)
        f_pval = pf(f_value, df1, df2, lower.tail = F)
        pval[i] = f_pval
        R_square[i] = 1 - RSS2/RSS1
      }
    }
    #==================================

    # Trapezoid rule for integration of GR_AOC
    GRavg = NULL
    for(j in 1:length(concs)) {
      data_trapz = data_exp[data_exp$concentration == concs[j],]
      GRavg[j] = mean(data_trapz$GR)
    }
    AOC[i] = sum((1 - (GRavg[1:(length(GRavg)-1)]+GRavg[2:length(GRavg)])/2)*diff(log10(concs), lag = 1))/(log10(concs[length(concs)]) - log10(concs[1]))
  }

  # Calculate GR50 from parameters
  parameters$GR50 = with(parameters, GEC50*((1-GRinf)/(0.5-GRinf) - 1)^(1/Hill))
  parameters$GRmax = GRmax
  parameters$GR_AOC = AOC
  parameters$r2 = R_square
  if(is.null(pval)) {pval = NA}
  parameters$pval = pval
  # Re-order rows to match reference_output
  parameters$experiment = experiments
  pcutoff = .05 # Threshold for F-test pval
  parameters$fit = with(parameters, ifelse(pval >= .05 | is.na(GEC50), "flat","sigmoid"))
  parameters$flat_fit = GR_mean
  # Add specified values for flat fits: GEC50 = 0, Hill = 0.01 and GR50 = +/- Inf
  for(i in 1:dim(parameters)[1]) {
    if(parameters$fit[i] == "flat") {
      parameters$GEC50[i] = 0
      parameters$Hill[i] = 0.01
      parameters$GR50[i] = ifelse(parameters$flat_fit[i] > .5, Inf, -Inf)
      parameters$GRinf[i] = parameters$flat_fit[i]
    }
  }
  # Add GR50 = +/-Inf for any curves that don't reach GR = 0.5
  for(i in 1:dim(parameters)[1]) {
    if(is.na(parameters$GR50[i])) {
      parameters$GR50[i] = ifelse(parameters$flat_fit[i] > .5, Inf, -Inf)
    }
  }
  parameters = parameters[,c('GR50','GRmax','GR_AOC','GEC50','GRinf','Hill', 'r2','pval','experiment', 'fit')]
  if(!is.null(metadata)) {
    parameters = cbind(metadata, parameters)
  }
  return(parameters)
}

.GRlogistic_3u = function(c){GRinf + (1 - GRinf)/(1 + (c/EC50)^Hill)}

#' Extract GR parameters from a dataset
#'
#' This function takes in a dataset with the columns 'cell_count__time0', 'cell_count__ctrl',
#' 'cell_count', 'concentration', and additional grouping variables and calculates growth rate
#' inhibition (GR) metrics for each experiment in the dataset.
#'
#' @param inputData a data table consisting of the columns 'cell_count__time0', 'cell_count__ctrl', cell_count', 'concentration', and additional grouping columns.
#' @param groupingColumns columns of inputData that identify the experiment (i.e. which columns to not average over)
#' @return a data table of GR metrics parameters (GR50, GRmax, etc.) as well as goodness of fit measures.
#' @author Nicholas Clark
#' @details
#' By default, only the table of GR metrics is returned. If GRtable equals 'both', then a list with two elements is returned: the table of GR metrics and a table of the initial data with an added column of GR values. If GRtable equals 'only', then only the table with GR values is returned.
#' @seealso \code{drm}
#' @export

GRfit = function(inputData, groupingColumns, GRtable = 'none') {
  gr_table = .GRcalculate(inputData, groupingColumns)
  if(GRtable == 'only') {
    return(gr_table)
  }
  parameter_table = .GRlogisticFit(gr_table, groupingColumns)
  if(GRtable == 'none'){
    return(parameter_table)
  } else if(GRtable == 'both') {
    return(list(gr_table,parameter_table))
  }
}
