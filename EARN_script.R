#---
#  title: "rENA ARM (aka EARN)"
# author: "Zhongtian Huang, Zeyneb Sarioglan & Andres F. Zambrano" (in alphabetical order)
# date: "2025-07-16"
#---


# 1. Install packages and load libraries

#You can download rENA from this link:
#  <https://cran.r-project.org/src/contrib/Archive/rENA/>
  

rm(list = ls())

library(rENA)
library(ggplot2)
library(plotly)
library(dplyr)
library(arules)
library(dplyr)
library(arulesViz)
library(purrr)
library(stringr)
library(randomcoloR)
library(igraph)
library(ggraph)

# 2. Read the functions

# There are two functions:
  
#**ena_arm_func** -- function for original plots of each group
#*/*ena_arm_func_diff** -- function for difference plots between each pair of groups


# function for original plots
ena_arm_func <- function(data, m1, m2, m3, unitCols, codeCols, conversationCols, rotation, rotationCol, MovStanzaSize, threshold, nsize, esize, xlim, ylim, dir, metaCols){
  
  # get all values for rotation column
  rotation_all <- unique(data[[rotationCol]])
  
  # set default values
  if (missing(metaCols)) {
    metaCols <- "metadata"
    i <- 0
    while (metaCols %in% names(data)) {
      i <- i + 1
      metaCols <- paste0("metadata", i)
    }
    data[[metaCols]] <- "none"}
  if (missing(nsize)) {nsize <- 1}
  if (missing(esize)) {esize <- 1}
  if (missing(xlim)) {xlim <- NULL}
  if (missing(ylim)) {ylim <- NULL}
  if (missing(dir)) {dir <- NULL}
  
  if (length(unitCols) == 1) {
    units <- data.frame(data[,unitCols])
    colnames(units)[1] <- unitCols[1]}
  else {units <- data[,unitCols]}
  
  if (length(conversationCols) == 1) {
    conversation <- data.frame(data[,conversationCols])
    colnames(conversation)[1] <- conversationCols[1]}
  else {conversation <- data[,conversationCols]}
  
  if (length(codeCols) == 1) {
    codes <- data.frame(data[,codeCols])
    colnames(codes)[1] <- codeCols[1]}
  else {codes <- data[,codeCols]}
  
  
  # fit ENA model
  accum.ena <- 
    ena.accumulate.data(
      units = units,
      conversation = conversation,
      metadata = data[,metaCols], # optional
      codes = codes,
      window.size.back = MovStanzaSize)
  
  set.ena <- 
    ena.make.set(
      enadata = accum.ena, # the accumulation run above
      rotation.by = ena.rotate.by.mean, # equivalent of mean=TRUE in the ena function
      rotation.params = lapply(rotation, function(cond) {
        eval(parse(text = paste0("accum.ena$meta.data$", rotationCol))) == cond
      })
    )
  
  lineweights_list <- list()
  for (i in 1:length(rotation_all)){lineweights_list[[i]] = as.matrix(eval(parse(text = paste0("set.ena$line.weights$", rotationCol, "$", rotation_all[i]))))}
  
  mean_list <- list()
  for (i in 1:length(rotation_all)){mean_list[[i]] = as.vector(colMeans(lineweights_list[[i]]))}
  
  nodes0 <- set.ena$rotation$nodes
  
  # get the coordinates for nodes
  func <- function(){
    x <- NULL
    y <- NULL
    xend <- NULL
    yend <- NULL
    midx <- NULL
    midy <- NULL
    rulesf <- NULL
    rulesr <- NULL
    
    for (i in 1:length(mean_list[[1]])){
      code1 <- strsplit(colnames(lineweights_list[[1]])[i], " & ")[[1]][1]
      code2 <- strsplit(colnames(lineweights_list[[1]])[i], " & ")[[1]][2]
      rulesf <- c(rulesf, paste0('{', code1, '} => {', code2, '}'))
      rulesr <- c(rulesr, paste0('{', code2, '} => {', code1, '}'))
      x <- c(x, nodes0[nodes0$code == code1, "MR1"]$MR1[1])
      y <- c(y, nodes0[nodes0$code == code1, "SVD2"]$SVD2[1])
      xend <- c(xend, nodes0[nodes0$code == code2, "MR1"]$MR1[1])
      yend <- c(yend, nodes0[nodes0$code == code2, "SVD2"]$SVD2[1])
      midx <- c(midx, (nodes0[nodes0$code == code1, "MR1"]$MR1[1]+nodes0[nodes0$code == code2, "MR1"]$MR1[1])/2)
      midy <- c(midy, (nodes0[nodes0$code == code1, "SVD2"]$SVD2[1]+nodes0[nodes0$code == code2, "SVD2"]$SVD2[1])/2)
    }
    return(list(v1=x, v2=y, v3=xend, v4=yend, v5=midx, v6=midy, v7=rulesf, v8=rulesr))}
  
  
  result2 <- func()
  x <- result2$v1
  y <- result2$v2
  xend <- result2$v3
  yend <- result2$v4
  midx <- result2$v5
  midy <- result2$v6
  rulesf <- result2$v7
  rulesr <- result2$v8
  
  # calculate ARM values
  ARM <- function(rotation){
    nodes <- set.ena$rotation$nodes
    Data <- data[data[rotationCol] == rotation, ]
    rownames(Data) <- NULL
    Data <- Data[codeCols]
    
    # Initialize an empty list to store sequences
    Sequences <- list()
    
    # Find sequences of 1s
    for (i in 1:(nrow(Data) - MovStanzaSize + 1)) {
      # Extract the moving stanza of data
      SubData <- Data[i:(i + MovStanzaSize - 1), ]
      # Get columns which contain at least one 1 in the stanza
      CurrentCodes <- which(colSums(SubData == 1) > 0)
      # Append the column names to the Sequences list
      Sequences[[length(Sequences) + 1]] <- names(Data)[CurrentCodes]
    }
    
    #Contains Support and Cosine
    rules0 <- apriori(Sequences, parameter = list(supp = 0, conf = 0))
    # Calculate additional interest measures
    interestMeasures <- interestMeasure(rules0, measure=c("cosine", "addedValue"), transactions=Transactions)
    
    # Adding calculated measures to the rules
    quality(rules0) <- cbind(quality(rules0), interestMeasures)
    rules0 <- as(rules0, "data.frame")
    write.csv(rules0, file = paste0(dir, rotation, "_all_rules.csv"), row.names = FALSE)
    
    rules <- rules0[rules0[['count']] > 0, ]
    rules <- rules[rules[[m2]] >= threshold, ]
    
    rules_1node <- rules0[grepl("\\{\\}", rules0$rules), ]
    rownames(rules_1node) <- NULL
    remove_chars <- function(x){return(substring(x, 8, nchar(x)-1))}
    rules_1node$code <- sapply(rules_1node$rules, remove_chars)
    rules_1node <- rules_1node[c('code', m1)]
    nodes <- merge(nodes, rules_1node, by = "code")
    nodes$sign <- ifelse(nodes[[m1]] >= 0, "positive", "negative")
    
    linenum1 <- NULL
    linenum2 <- NULL
    linenum3 <- NULL
    for (i in 1:length(rulesf)){
      linenum1 <- c(linenum1, rules0[rules0$rules == rulesf[i], m2])
      linenum2 <- c(linenum2, rules0[rules0$rules == rulesf[i], m3])
      linenum3 <- c(linenum3, rules0[rules0$rules == rulesr[i], m3])
    }
    return(list(v1=linenum1, v2=linenum2, v3=linenum3, nodes=nodes, rules=rules))
  }
  
  # parse ARM values into a list
  node_func <- function(rules_df){
    rules_df <- rules_df[!grepl("\\{\\}", rules_df$rules), ]
    rownames(rules_df) <- NULL
    n <- nrow(rules_df)
    result <- list()
    ivec <- NULL
    j <- 2
    while (TRUE){
      result_temp1 <- list()
      for (k in 1:(j-1)){
        result_temp2 <- list()
        for (i in 1:length(rules_df$rules)){
          str_list <- strsplit(rules_df[i,'rules'], " => ")
          sup <- rules_df[i,m2]
          conf <- rules_df[i,m3]
          comma_count1 <- str_count(str_list[[1]][1], ",")
          comma_count2 <- str_count(str_list[[1]][2], ",")
          if (comma_count1 == j-k-1 & comma_count2 == k-1){
            if (length(result_temp2) == 0){
              for (l in 1:(j+2)){
                if (l <= j-k){result_temp2[[l]] <- strsplit(substr(str_list[[1]][1], 2, nchar(str_list[[1]][1]) - 1), ",")[[1]][l]}
                if (l > j-k & l <= j){result_temp2[[l]] <- strsplit(substr(str_list[[1]][2], 2, nchar(str_list[[1]][2]) - 1), ",")[[1]][l-j+k]}
                if (l == j+1){result_temp2[[l]] <- sup}
                if (l == j+2){result_temp2[[l]] <- conf}
              }}
            else {
              for (l in 1:(j+2)){
                if (l <= j-k){result_temp2[[l]] <- c(result_temp2[[l]], strsplit(substr(str_list[[1]][1], 2, nchar(str_list[[1]][1]) - 1), ",")[[1]][l])}
                if (l > j-k & l <= j){result_temp2[[l]] <- c(result_temp2[[l]], strsplit(substr(str_list[[1]][2], 2, nchar(str_list[[1]][2]) - 1), ",")[[1]][l-j+k])}
                if (l == j+1){result_temp2[[l]] <- c(result_temp2[[l]], sup)}
                if (l == j+2){result_temp2[[l]] <- c(result_temp2[[l]], conf)}
              }}
            ivec <- c(ivec, i)
          }}
        if (length(result_temp2) != 0){result_temp1[[k]] <- result_temp2}
      } 
      if (length(result_temp1) != 0){result[[j]] <- result_temp1}
      if (length(ivec) == n){break}
      else {j <- j + 1}
    } 
    return(result)
  }
  
  # generate plots
  plot_func <- function(N, rotation){
    plot_list <- list()
    if (N == 2){
      textvec <- NULL
      
      df <- data.frame(x=x, y=y, xend=xend, yend=yend, linewidth=linenum1, midx=midx, midy=midy, num1=linenum1, num2=linenum2, num3=linenum3)
      df <- df[df[['linewidth']] >= threshold & !is.na(df[['linewidth']]), ]
      df$sign <- ifelse(df$num1 >= 0, "positive", "negative")
      rownames(df) <- NULL
      
      for (i in 1:length(df$num1)){
        if (df$num2[i] >= df$num3[i]){
          text <- paste0(m2, ': ',as.character(round(df$num1[i],3)),', ', m3, '(forward): ',as.character(round(df$num2[i],3)),', ', m3, '(reverse): ',as.character(round(df$num3[i],3)))
          textvec <- c(textvec, text)}
        else {
          text <- paste0(m2, ': ',as.character(round(df$num1[i],3)),', ', m3, '(forward): ',as.character(round(df$num3[i],3)),', ', m3, '(reverse): ',as.character(round(df$num2[i],3)))
          textvec <- c(textvec, text)}
      }
      
      df$textvec <- textvec
      
      plot <- ggplot() +
        geom_point(data=nodes, aes(x = MR1, y = SVD2, size = .data[[m1]]*nsize, text = code, color = sign)) +
        geom_text(
          data = nodes,
          aes(x = MR1, y = SVD2, label = code),
          vjust = -1.2,
          hjust = 0.5,
          size = 2.5
        ) +
        geom_segment(data=df, aes(x = x, y = y, xend = xend, yend = yend, size=linewidth*esize, color = sign)) +
        scale_color_manual(values = c("positive" = "red", "negative" = "blue")) +
        scale_size_continuous(range = c(0.5, 5)) +
        geom_text(data=df, aes(x = midx, y = midy, label = textvec), color = 'transparent') +
        coord_equal(clip = "off", 
                    xlim = c(min(nodes$MR1) - 0.8, max(nodes$MR1) + 0.8),
                    ylim = c(min(nodes$SVD2) - 0.8, max(nodes$SVD2) + 0.8)) +
        ggtitle(paste0(rotation, " 2-Node Rules")) +
        theme(legend.position = "none")
      
      p <- ggplotly(plot, tooltip = "text")
      arrow_annotations <- lapply(1:nrow(df), function(i) {
        if (df$num2[i] >= df$num3[i]){
          list(
            x = df$midx[i],
            y = df$midy[i],
            ax = df$midx[i]-0.05*(df$xend[i]-df$x[i]),
            ay = df$midy[i]-0.05*(df$yend[i]-df$y[i]),
            text = "",
            xref = "x",     # important: refer to the data coordinate system
            yref = "y",
            axref = "x",
            ayref = "y",
            showarrow = TRUE,
            arrowhead = 5,
            arrowsize = 1,
            arrowcolor = "red"
          )}
        else {
          list(
            x = df$midx[i],
            y = df$midy[i],
            ax = df$midx[i]-0.05*(df$x[i]-df$xend[i]),
            ay = df$midy[i]-0.05*(df$y[i]-df$yend[i]),
            text = "",
            xref = "x",     # important: refer to the data coordinate system
            yref = "y",
            axref = "x",
            ayref = "y",
            showarrow = TRUE,
            arrowhead = 5,
            arrowsize = 1,
            arrowcolor = "red"
          )}
      })
      
      # Add annotations to plotly layout
      p <- p %>%
        layout(annotations = arrow_annotations)
      plot_list[[1]] <- p
    }
    else {
      for (k in 1:length(result[[N]])){
        xlist <- list()
        ylist <- list()
        nlist <- NULL
        label_vec <- NULL
        check_dict <- list()
        check_num <- NULL
        N_vec <- length(result[[N]][[k]][[1]])
        for (i in 1:N_vec){
          check_str <- NULL
          for (l in 1:N){check_str <- c(check_str, result[[N]][[k]][[l]][i])}
          check_str <- sort(check_str)
          if (!any(sapply(check_dict, function(x) identical(x, check_str)))){
            check_num <- c(check_num, result[[N]][[k]][[N+2]][i])
            check_dict[[length(check_dict)+1]] <- check_str
            for (j in 1:N){
              if (i == 1){
                xlist[[j]] <- nodes[nodes$code == result[[N]][[k]][[j]][i], "MR1"]$MR1[1]
                ylist[[j]] <- nodes[nodes$code == result[[N]][[k]][[j]][i], "SVD2"]$SVD2[1]
              }
              else {
                xlist[[j]] <- c(xlist[[j]], nodes[nodes$code == result[[N]][[k]][[j]][i], "MR1"]$MR1[1])
                ylist[[j]] <- c(ylist[[j]], nodes[nodes$code == result[[N]][[k]][[j]][i], "SVD2"]$SVD2[1])
              }}
            label_vec <- c(label_vec, result[[N]][[k]][[N+1]][i])
            for (m in 1:N){
              if (m == 1){
                nlist <- c(nlist, result[[N]][[k]][[m]][i])
              }
              else if (m == N-k+1) {
                nlist[length(nlist)] <- paste(nlist[length(nlist)], result[[N]][[k]][[m]][i], sep=" => ")
              }
              else {
                nlist[length(nlist)] <- paste(nlist[length(nlist)], result[[N]][[k]][[m]][i], sep=", ")
              }
              label_vec[length(nlist)] <- paste0(
                nlist[length(nlist)], "\n",
                m2, " = ", round(result[[N]][[k]][[N+1]][i], 3), ", ",
                m3, " = ", round(result[[N]][[k]][[N+2]][i], 3)
              )
            }}
          else if (result[[N]][[k]][[N+2]][i] > check_num[which(sapply(check_dict, function(x) identical(x, check_str)))]) {
            index <- which(sapply(check_dict, function(x) identical(x, check_str)))
            check_num[index] <- result[[N]][[k]][[N+2]][i]
            check_dict[[index]] <- check_str
            for (j in 1:N){
              xlist[[j]][index] <- nodes[nodes$code == result[[N]][[k]][[j]][i], "MR1"]$MR1[1]
              ylist[[j]][index] <- nodes[nodes$code == result[[N]][[k]][[j]][i], "SVD2"]$SVD2[1]
            }
            label_vec[index] <- result[[N]][[k]][[N+1]][i]
            for (m in 1:N){
              if (m == 1){
                nlist[index] <- result[[N]][[k]][[m]][i]
              }
              else if (m == N-k+1) {
                nlist[index] <- paste(nlist[index], result[[N]][[k]][[m]][i], sep=" => ")
              }
              else {
                nlist[index] <- paste(nlist[index], result[[N]][[k]][[m]][i], sep=", ")
              }
              label_vec[length(nlist)] <- paste0(
                nlist[length(nlist)], "\n",
                m2, " = ", round(result[[N]][[k]][[N+1]][i], 3), ", ",
                m3, " = ", round(result[[N]][[k]][[N+2]][i], 3)
              )
            }
          }
        }
        
        N_vec <- length(label_vec)
        set.seed(123)
        conn3df <- as.data.frame(c(xlist, ylist))
        colnames(conn3df) <- c(paste("xs", 1:(N-k), sep = ""), paste("xe", 1:k, sep = ""), paste("ys", 1:(N-k), sep = ""), paste("ye", 1:k, sep = ""))
        col_names <- colnames(conn3df)
        curvature_values <- runif(N_vec, min = 0.1, max = 1)
        conn3df$curvature <- curvature_values
        colors <- distinctColorPalette(N_vec)
        conn3df$colors <- colors
        conn3df$values <- factor(label_vec)
        
        curve_layers <- list()
        xs <- grep("^xs", col_names, value = TRUE)
        xe <- grep("^xe", col_names, value = TRUE)
        ys <- grep("^ys", col_names, value = TRUE)
        ye <- grep("^ye", col_names, value = TRUE)
        xs_xe <- paste(expand.grid(xs, xe)[,1], expand.grid(xs, xe)[,2], sep = "_")
        ys_ye <- paste(expand.grid(ys, ye)[,1], expand.grid(ys, ye)[,2], sep = "_")
        
        for (l in 1:(k*(N-k))){
          x <- strsplit(xs_xe[l], "_")[[1]][1]
          xend <- strsplit(xs_xe[l], "_")[[1]][2]
          y <- strsplit(ys_ye[l], "_")[[1]][1]
          yend <- strsplit(ys_ye[l], "_")[[1]][2]
          conn3df_temp <- conn3df[c(x, y, xend, yend, "curvature", "values")]
          colnames(conn3df_temp) <- c("x", "y", "xend", "yend", "curvature", "values")
          
          curve_layers[[l]] <- pmap(conn3df_temp, function(x, y, xend, yend, curvature, values) {
            geom_curve(
              aes(x = x, y = y, xend = xend, yend = yend, color = values),
              curvature = curvature,
              arrow = arrow(length = unit(0.15, "inches"), type = "closed"),
              linewidth = 1
            )
          })
        }
        
        plot_list[[k]]  <- ggplot() + 
          geom_point(data=nodes, aes(x = MR1, y = SVD2, size = .data[[m1]]), color = "black") +
          geom_text(
            data = nodes,
            aes(x = MR1, y = SVD2, label = code),
            vjust = -1.2,
            hjust = 0.5,
            size = 2.5
          ) +
          curve_layers +
          scale_color_manual(name = NULL, values = colors) +
          guides(
            color = guide_legend(
              ncol = 1,          # one rule per line
              byrow = TRUE
            ),
            size = guide_legend(
              ncol = 6           # keeps support dots compact
            )
          ) +
          coord_equal(
            clip = "off",
            xlim = c(min(nodes$MR1) - 0.8, max(nodes$MR1) + 0.8),
            ylim = c(min(nodes$SVD2) - 0.8, max(nodes$SVD2) + 0.8)
          ) +
          ggtitle(paste0(rotation, " ", N, "-Node Rules Plot ", k)) +
          theme_minimal() +
          theme(
            legend.box = "vertical",
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 7),
            legend.key.size = unit(0.4, "cm"),
            legend.spacing.y = unit(0.1, "cm"),
            axis.title = element_text(size = 7),
            axis.text = element_text(size = 6)
          )
      }}
    return(plot_list)
  }
  
  # run the functions above
  for (k in 1:length(rotation_all)){
    element <- rotation_all[k]
    ARMresult <- ARM(element)
    linenum1 <- ARMresult$v1
    linenum2 <- ARMresult$v2
    linenum3 <- ARMresult$v3
    nodes <- ARMresult$nodes
    rules <- ARMresult$rules
    cat("\n\n\n")
    result <- tryCatch({node_func(rules)}, error = function(e) {
      cat(paste0("No rule was found for ", element, ". \n\n\n\n"))
      return(NULL)
    })
    if (is.null(result)) {next}
    
    for (i in 2:5){
      tryCatch({
        plot_list <- plot_func(i, element)
        
        if (i == 2) {panel <- "'Viewer' panel"}
        else {panel <- "'Plots' panel"}
        
        for (j in length(plot_list)) {
          plot.new()
          print(plot_list[[j]])  
          cat(paste0("Plot ", j, " for ", element, " ", i, "-node rules shown. Press [Enter] to continue...\n(The plot is shown in the ", panel, ")\n"))
          readline()
        }}, error = function(e) {
          cat(paste0("There are no ", i, "-node rules for ", element, ". Press [Enter] to continue...\n"))
          readline()
        })
    }
  }
}

# function for different plots
ena_arm_func_diff <- function(data, m1, m2, m3, unitCols, codeCols, conversationCols, rotation, rotationCol, MovStanzaSize, threshold, nsize, esize, xlim, ylim, metaCols){
  
  # get all values for rotation column
  rotation_all <- unique(data[[rotationCol]])
  
  # set default values
  if (missing(metaCols)) {
    metaCols <- "metadata"
    i <- 0
    while (metaCols %in% names(data)) {
      i <- i + 1
      metaCols <- paste0("metadata", i)
    }
    data[[metaCols]] <- "none"}
  if (missing(nsize)) {nsize <- 1}
  if (missing(esize)) {esize <- 1}
  if (missing(xlim)) {xlim <- NULL}
  if (missing(ylim)) {ylim <- NULL}
  
  if (length(unitCols) == 1) {
    units <- data.frame(data[,unitCols])
    colnames(units)[1] <- unitCols[1]}
  else {units <- data[,unitCols]}
  
  if (length(conversationCols) == 1) {
    conversation <- data.frame(data[,conversationCols])
    colnames(conversation)[1] <- conversationCols[1]}
  else {conversation <- data[,conversationCols]}
  
  if (length(codeCols) == 1) {
    codes <- data.frame(data[,codeCols])
    colnames(codes)[1] <- codeCols[1]}
  else {codes <- data[,codeCols]}
  
  # fit ENA model
  accum.ena <- 
    ena.accumulate.data(
      units = units,
      conversation = conversation,
      metadata = data[,metaCols], # optional
      codes = codes,
      window.size.back = MovStanzaSize)
  
  set.ena <- 
    ena.make.set(
      enadata = accum.ena, # the accumulation run above
      rotation.by = ena.rotate.by.mean, # equivalent of mean=TRUE in the ena function
      rotation.params = lapply(rotation, function(cond) {
        eval(parse(text = paste0("accum.ena$meta.data$", rotationCol))) == cond
      })
    )
  
  lineweights_list <- list()
  for (i in 1:length(rotation_all)){lineweights_list[[i]] = as.matrix(eval(parse(text = paste0("set.ena$line.weights$", rotationCol, "$", rotation_all[i]))))}
  
  mean_list <- list()
  for (i in 1:length(rotation_all)){mean_list[[i]] = as.vector(colMeans(lineweights_list[[i]]))}
  
  nodes0 <- set.ena$rotation$nodes
  
  # get coordinates for nodes
  func <- function(){
    x <- NULL
    y <- NULL
    xend <- NULL
    yend <- NULL
    midx <- NULL
    midy <- NULL
    rulesf <- NULL
    rulesr <- NULL
    
    for (i in 1:length(mean_list[[1]])){
      code1 <- strsplit(colnames(lineweights_list[[1]])[i], " & ")[[1]][1]
      code2 <- strsplit(colnames(lineweights_list[[1]])[i], " & ")[[1]][2]
      rulesf <- c(rulesf, paste0('{', code1, '} => {', code2, '}'))
      rulesr <- c(rulesr, paste0('{', code2, '} => {', code1, '}'))
      x <- c(x, nodes0[nodes0$code == code1, "MR1"]$MR1[1])
      y <- c(y, nodes0[nodes0$code == code1, "SVD2"]$SVD2[1])
      xend <- c(xend, nodes0[nodes0$code == code2, "MR1"]$MR1[1])
      yend <- c(yend, nodes0[nodes0$code == code2, "SVD2"]$SVD2[1])
      midx <- c(midx, (nodes0[nodes0$code == code1, "MR1"]$MR1[1]+nodes0[nodes0$code == code2, "MR1"]$MR1[1])/2)
      midy <- c(midy, (nodes0[nodes0$code == code1, "SVD2"]$SVD2[1]+nodes0[nodes0$code == code2, "SVD2"]$SVD2[1])/2)
    }
    return(list(v1=x, v2=y, v3=xend, v4=yend, v5=midx, v6=midy, v7=rulesf, v8=rulesr))}
  
  
  result2 <- func()
  x <- result2$v1
  y <- result2$v2
  xend <- result2$v3
  yend <- result2$v4
  midx <- result2$v5
  midy <- result2$v6
  rulesf <- result2$v7
  rulesr <- result2$v8
  
  # calculate ARM values
  ARM <- function(rotation, a, b){
    rules <- list()
    nodes <- set.ena$rotation$nodes
    
    Data <- data[data[rotationCol] == rotation[a], ]
    rownames(Data) <- NULL
    Data <- Data[codeCols]
    
    # Initialize an empty list to store sequences
    Sequences <- list()
    
    # Find sequences of 1s
    for (i in 1:(nrow(Data) - MovStanzaSize + 1)) {
      # Extract the moving stanza of data
      SubData <- Data[i:(i + MovStanzaSize - 1), ]
      # Get columns which contain at least one 1 in the stanza
      CurrentCodes <- which(colSums(SubData == 1) > 0)
      # Append the column names to the Sequences list
      Sequences[[length(Sequences) + 1]] <- names(Data)[CurrentCodes]
    }
    
    #Contains Support and Cosine
    rules[[1]] <- apriori(Sequences, parameter = list(supp = 0, conf = 0))
    # Calculate additional interest measures
    interestMeasures <- interestMeasure(rules[[1]], measure=c("cosine", "addedValue"), transactions=Transactions)
    
    # Adding calculated measures to the rules
    quality(rules[[1]]) <- cbind(quality(rules[[1]]), interestMeasures)
    rules[[1]] <- as(rules[[1]], "data.frame")
    
    Data <- data[data[rotationCol] == rotation[b], ]
    rownames(Data) <- NULL
    Data <- Data[codeCols]
    
    # Initialize an empty list to store sequences
    Sequences <- list()
    
    # Find sequences of 1s
    for (i in 1:(nrow(Data) - MovStanzaSize + 1)) {
      # Extract the moving stanza of data
      SubData <- Data[i:(i + MovStanzaSize - 1), ]
      # Get columns which contain at least one 1 in the stanza
      CurrentCodes <- which(colSums(SubData == 1) > 0)
      # Append the column names to the Sequences list
      Sequences[[length(Sequences) + 1]] <- names(Data)[CurrentCodes]
    }
    
    #Contains Support and Cosine
    rules[[2]] <- apriori(Sequences, parameter = list(supp = 0, conf = 0))
    # Calculate additional interest measures
    interestMeasures <- interestMeasure(rules[[2]], measure=c("cosine", "addedValue"), transactions=Transactions)
    
    # Adding calculated measures to the rules
    quality(rules[[2]]) <- cbind(quality(rules[[2]]), interestMeasures)
    rules[[2]] <- as(rules[[2]], "data.frame")
    
    rules0 <- merge(rules[[1]][c('rules', m2, m3, 'count')], rules[[2]][c('rules', m2, m3, 'count')], by='rules')
    rules0[[m2]] <- rules0[[paste0(m2, '.x')]] - rules0[[paste0(m2, '.y')]]
    rules0[[m3]] <- rules0[[paste0(m3, '.x')]] - rules0[[paste0(m3, '.y')]]
    rules0 <- rules0[c('rules', m2, m3, 'count.x', 'count.y')]
    
    rules00 <- merge(rules[[1]][c('rules', m1)], rules[[2]][c('rules', m1)], by='rules')
    rules00[[m1]] <- rules00[[paste0(m1, '.x')]] - rules00[[paste0(m1, '.y')]]
    
    rules_1node <- rules00[grepl("\\{\\}", rules00$rules), ]
    rownames(rules_1node) <- NULL
    remove_chars <- function(x){return(substring(x, 8, nchar(x)-1))}
    rules_1node$code <- sapply(rules_1node$rules, remove_chars)
    rules_1node <- rules_1node[c('code', m1)]
    nodes <- merge(nodes, rules_1node, by = "code")
    nodes$sign <- ifelse(nodes[[m1]] >= 0, "positive", "negative")
    
    linenum1 <- NULL
    linenum2 <- NULL
    linenum3 <- NULL
    for (i in 1:length(rulesf)){
      linenum1 <- c(linenum1, rules0[rules0$rules == rulesf[i], m2])
      linenum2 <- c(linenum2, rules0[rules0$rules == rulesf[i], m3])
      linenum3 <- c(linenum3, rules0[rules0$rules == rulesr[i], m3])
    }
    return(list(v1=linenum1, v2=linenum2, v3=linenum3, nodes=nodes, rules0=rules0))
  }
  
  # parse ARM values into a list
  node_func <- function(rules_df){
    rules_df <- rules_df[!grepl("\\{\\}", rules_df$rules), ]
    rownames(rules_df) <- NULL
    n <- nrow(rules_df)
    result <- list()
    ivec <- NULL
    j <- 2
    while (TRUE){
      result_temp1 <- list()
      for (k in 1:(j-1)){
        result_temp2 <- list()
        for (i in 1:length(rules_df$rules)){
          str_list <- strsplit(rules_df[i,'rules'], " => ")
          sup <- rules_df[i,m2]
          conf <- rules_df[i,m3]
          comma_count1 <- str_count(str_list[[1]][1], ",")
          comma_count2 <- str_count(str_list[[1]][2], ",")
          if (comma_count1 == j-k-1 & comma_count2 == k-1){
            if (length(result_temp2) == 0){
              for (l in 1:(j+2)){
                if (l <= j-k){result_temp2[[l]] <- strsplit(substr(str_list[[1]][1], 2, nchar(str_list[[1]][1]) - 1), ",")[[1]][l]}
                if (l > j-k & l <= j){result_temp2[[l]] <- strsplit(substr(str_list[[1]][2], 2, nchar(str_list[[1]][2]) - 1), ",")[[1]][l-j+k]}
                if (l == j+1){result_temp2[[l]] <- sup}
                if (l == j+2){result_temp2[[l]] <- conf}
              }}
            else {
              for (l in 1:(j+2)){
                if (l <= j-k){result_temp2[[l]] <- c(result_temp2[[l]], strsplit(substr(str_list[[1]][1], 2, nchar(str_list[[1]][1]) - 1), ",")[[1]][l])}
                if (l > j-k & l <= j){result_temp2[[l]] <- c(result_temp2[[l]], strsplit(substr(str_list[[1]][2], 2, nchar(str_list[[1]][2]) - 1), ",")[[1]][l-j+k])}
                if (l == j+1){result_temp2[[l]] <- c(result_temp2[[l]], sup)}
                if (l == j+2){result_temp2[[l]] <- c(result_temp2[[l]], conf)}
              }}
            ivec <- c(ivec, i)
          }}
        if (length(result_temp2) != 0){result_temp1[[k]] <- result_temp2}
      } 
      if (length(result_temp1) != 0){result[[j]] <- result_temp1}
      if (length(ivec) == n){break}
      else {j <- j + 1}
    } 
    return(result)
  }
  
  # generate plots
  plot_func <- function(N, rotation){
    plot_list <- list()
    if (N == 2){
      textvec <- NULL
      
      df <- data.frame(x=x, y=y, xend=xend, yend=yend, linewidth=abs(linenum1), midx=midx, midy=midy, num1=linenum1, num2=linenum2, num3=linenum3)
      df <- df[df[['linewidth']] >= threshold & !is.na(df[['linewidth']]), ]
      df$sign <- ifelse(df$num1 >= 0, "positive", "negative")
      rownames(df) <- NULL
      
      for (i in 1:length(df$num1)){
        if (abs(df$num2[i]) >= abs(df$num3[i])){
          text <- paste0(m2, ': ',as.character(round(df$num1[i],3)),', ', m3, '(forward): ',as.character(round(df$num2[i],3)),', ', m3, '(reverse): ',as.character(round(df$num3[i],3)))
          textvec <- c(textvec, text)}
        else {
          text <- paste0(m2, ': ',as.character(round(df$num1[i],3)),', ', m3, '(forward): ',as.character(round(df$num3[i],3)),', ', m3, '(reverse): ',as.character(round(df$num2[i],3)))
          textvec <- c(textvec, text)}
      }
      
      df$textvec <- textvec
      
      plot <- ggplot() +
        geom_point(data=nodes, aes(x = MR1, y = SVD2, size = abs(.data[[m1]])*nsize, text = paste0(code, ': ', round(.data[[m1]], 3)), color = sign)) +
        geom_text(data=nodes, aes(x = MR1, y = SVD2+0.15, label = code), vjust = -1.2, hjust = 0.5, size = 3) +
        geom_segment(data=df, aes(x = x, y = y, xend = xend, yend = yend, size=linewidth*esize, color = sign)) +
        scale_color_manual(values = c("positive" = "red", "negative" = "blue")) +
        scale_size_continuous(range = c(0.5, 5)) +
        geom_text(data=df, aes(x = midx, y = midy, label = textvec), color = 'transparent') +
        coord_equal(clip = "off", 
                    xlim = c(min(nodes$MR1) - 0.8, max(nodes$MR1) + 0.8),
                    ylim = c(min(nodes$SVD2) - 0.8, max(nodes$SVD2) + 0.8)) +
        ggtitle(paste0(rotation, " 2-Node Rules")) +
        theme(legend.position = "none")
      
      p <- ggplotly(plot, tooltip = "text")
      arrow_annotations <- lapply(1:nrow(df), function(i) {
        if (abs(df$num2[i]) >= abs(df$num3[i]) & df$num1[i] >= 0){
          list(
            x = df$midx[i],
            y = df$midy[i],
            ax = df$midx[i]-0.05*(df$xend[i]-df$x[i]),
            ay = df$midy[i]-0.05*(df$yend[i]-df$y[i]),
            text = "",
            xref = "x",     # important: refer to the data coordinate system
            yref = "y",
            axref = "x",
            ayref = "y",
            showarrow = TRUE,
            arrowhead = 5,
            arrowsize = 1,
            arrowcolor = "red"
          )}
        else if (abs(df$num2[i]) >= abs(df$num3[i]) & df$num1[i] < 0){
          list(
            x = df$midx[i],
            y = df$midy[i],
            ax = df$midx[i]-0.05*(df$xend[i]-df$x[i]),
            ay = df$midy[i]-0.05*(df$yend[i]-df$y[i]),
            text = "",
            xref = "x",     # important: refer to the data coordinate system
            yref = "y",
            axref = "x",
            ayref = "y",
            showarrow = TRUE,
            arrowhead = 5,
            arrowsize = 1,
            arrowcolor = "blue"
          )}
        else if (abs(df$num2[i]) < abs(df$num3[i]) & df$num1[i] >= 0){
          list(
            x = df$midx[i],
            y = df$midy[i],
            ax = df$midx[i]-0.05*(df$x[i]-df$xend[i]),
            ay = df$midy[i]-0.05*(df$y[i]-df$yend[i]),
            text = "",
            xref = "x",     # important: refer to the data coordinate system
            yref = "y",
            axref = "x",
            ayref = "y",
            showarrow = TRUE,
            arrowhead = 5,
            arrowsize = 1,
            arrowcolor = "red"
          )}
        else {
          list(
            x = df$midx[i],
            y = df$midy[i],
            ax = df$midx[i]-0.05*(df$x[i]-df$xend[i]),
            ay = df$midy[i]-0.05*(df$y[i]-df$yend[i]),
            text = "",
            xref = "x",     # important: refer to the data coordinate system
            yref = "y",
            axref = "x",
            ayref = "y",
            showarrow = TRUE,
            arrowhead = 5,
            arrowsize = 1,
            arrowcolor = "blue"
          )}
      })
      
      # Add annotations to plotly layout
      p <- p %>%
        layout(annotations = arrow_annotations)
      plot_list[[1]] <- p
    }
    else {
      for (k in 1:length(result[[N]])){
        xlist <- list()
        ylist <- list()
        nlist <- NULL
        label_vec <- NULL
        check_dict <- list()
        check_num <- NULL
        N_vec <- length(result[[N]][[k]][[1]])
        for (i in 1:N_vec){
          check_str <- NULL
          for (l in 1:N){check_str <- c(check_str, result[[N]][[k]][[l]][i])}
          check_str <- sort(check_str)
          if (!any(sapply(check_dict, function(x) identical(x, check_str)))){
            check_num <- c(check_num, result[[N]][[k]][[N+2]][i])
            check_dict[[length(check_dict)+1]] <- check_str
            for (j in 1:N){
              if (i == 1){
                xlist[[j]] <- nodes[nodes$code == result[[N]][[k]][[j]][i], "MR1"]$MR1[1]
                ylist[[j]] <- nodes[nodes$code == result[[N]][[k]][[j]][i], "SVD2"]$SVD2[1]
              }
              else {
                xlist[[j]] <- c(xlist[[j]], nodes[nodes$code == result[[N]][[k]][[j]][i], "MR1"]$MR1[1])
                ylist[[j]] <- c(ylist[[j]], nodes[nodes$code == result[[N]][[k]][[j]][i], "SVD2"]$SVD2[1])
              }}
            label_vec <- c(label_vec, result[[N]][[k]][[N+1]][i])
            for (m in 1:N){
              if (m == 1){
                nlist <- c(nlist, result[[N]][[k]][[m]][i])
              }
              else if (m == N-k+1) {
                nlist[length(nlist)] <- paste(nlist[length(nlist)], result[[N]][[k]][[m]][i], sep=" => ")
              }
              else {
                nlist[length(nlist)] <- paste(nlist[length(nlist)], result[[N]][[k]][[m]][i], sep=", ")
              }
              label_vec[length(nlist)] <- paste0(
                nlist[length(nlist)], "\n",
                m2, " = ", round(result[[N]][[k]][[N+1]][i], 3), ", ",
                m3, " = ", round(result[[N]][[k]][[N+2]][i], 3)
              )
            }}
          else if (result[[N]][[k]][[N+2]][i] > check_num[which(sapply(check_dict, function(x) identical(x, check_str)))]) {
            index <- which(sapply(check_dict, function(x) identical(x, check_str)))
            check_num[index] <- result[[N]][[k]][[N+2]][i]
            check_dict[[index]] <- check_str
            for (j in 1:N){
              xlist[[j]][index] <- nodes[nodes$code == result[[N]][[k]][[j]][i], "MR1"]$MR1[1]
              ylist[[j]][index] <- nodes[nodes$code == result[[N]][[k]][[j]][i], "SVD2"]$SVD2[1]
            }
            label_vec[index] <- result[[N]][[k]][[N+1]][i]
            for (m in 1:N){
              if (m == 1){
                nlist[index] <- result[[N]][[k]][[m]][i]
              }
              else if (m == N-k+1) {
                nlist[index] <- paste(nlist[index], result[[N]][[k]][[m]][i], sep=" => ")
              }
              else {
                nlist[index] <- paste(nlist[index], result[[N]][[k]][[m]][i], sep=", ")
              }
              label_vec[length(nlist)] <- paste0(
                nlist[length(nlist)], "\n",
                m2, " = ", round(result[[N]][[k]][[N+1]][i], 3), ", ",
                m3, " = ", round(result[[N]][[k]][[N+2]][i], 3)
              )
            }
          }
        }
        
        N_vec <- length(label_vec)
        set.seed(123)
        conn3df <- as.data.frame(c(xlist, ylist))
        colnames(conn3df) <- c(paste("xs", 1:(N-k), sep = ""), paste("xe", 1:k, sep = ""), paste("ys", 1:(N-k), sep = ""), paste("ye", 1:k, sep = ""))
        col_names <- colnames(conn3df)
        curvature_values <- runif(N_vec, min = 0.1, max = 1)
        conn3df$curvature <- curvature_values
        colors <- distinctColorPalette(N_vec)
        conn3df$colors <- colors
        conn3df$values <- factor(label_vec)
        
        curve_layers <- list()
        xs <- grep("^xs", col_names, value = TRUE)
        xe <- grep("^xe", col_names, value = TRUE)
        ys <- grep("^ys", col_names, value = TRUE)
        ye <- grep("^ye", col_names, value = TRUE)
        xs_xe <- paste(expand.grid(xs, xe)[,1], expand.grid(xs, xe)[,2], sep = "_")
        ys_ye <- paste(expand.grid(ys, ye)[,1], expand.grid(ys, ye)[,2], sep = "_")
        
        for (l in 1:(k*(N-k))){
          x <- strsplit(xs_xe[l], "_")[[1]][1]
          xend <- strsplit(xs_xe[l], "_")[[1]][2]
          y <- strsplit(ys_ye[l], "_")[[1]][1]
          yend <- strsplit(ys_ye[l], "_")[[1]][2]
          conn3df_temp <- conn3df[c(x, y, xend, yend, "curvature", "values")]
          colnames(conn3df_temp) <- c("x", "y", "xend", "yend", "curvature", "values")
          
          curve_layers[[l]] <- pmap(conn3df_temp, function(x, y, xend, yend, curvature, values) {
            geom_curve(
              aes(x = x, y = y, xend = xend, yend = yend, color = values),
              curvature = curvature,
              arrow = arrow(length = unit(0.15, "inches"), type = "closed"),
              linewidth = 1
            )
          })
        }
        
        plot_list[[k]]  <- ggplot() + 
          geom_point(data=nodes, aes(x = MR1, y = SVD2, size = abs(.data[[m1]])), color = "black") +
          geom_text(data=nodes, aes(x = MR1, y = SVD2, label = paste0(code, ': ', round(.data[[m1]], 3)), vjust = -1, hjust = 0.5, size = 3)) +
          curve_layers +
          #scale_fill_manual(values = c("positive" = "black", "negative" = "darkblue")) +
          scale_color_manual(name = NULL, values = colors) +
          guides(
            color = guide_legend(
              ncol = 1,          # one rule per line
              byrow = TRUE
            ),
            size = guide_legend(
              ncol = 6           # keeps support dots compact
            )
          ) +
          coord_equal(
            clip = "off",
            xlim = c(min(nodes$MR1) - 0.8, max(nodes$MR1) + 0.8),
            ylim = c(min(nodes$SVD2) - 0.8, max(nodes$SVD2) + 0.8)
          ) +
          ggtitle(paste0(rotation, " ", N, "-Node Rules Plot ", k)) +
          theme_minimal() +
          theme(
            legend.box = "vertical",
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 7),
            legend.key.size = unit(0.4, "cm"),
            legend.spacing.y = unit(0.1, "cm"),
            axis.title = element_text(size = 7),
            axis.text = element_text(size = 6)
          )
      }}
    return(plot_list)
  }
  
  # run the functions above
  for (k in 1:length(rotation_all)){
    for (m in 1:length(rotation_all)){
      if (k<m){
        element <- paste0(rotation_all[k], '-', rotation_all[m])
        ARMresult <- ARM(rotation_all, k, m)
        linenum1 <- ARMresult$v1
        linenum2 <- ARMresult$v2
        linenum3 <- ARMresult$v3
        nodes <- ARMresult$nodes
        rules0 <- ARMresult$rules0
        rules <- rules0[rules0[['count.x']] > 0 & rules0[['count.y']] > 0, ]
        rules <- rules[abs(rules[m2]) >= threshold, ]
        cat("\n\n\n")
        result <- tryCatch({node_func(rules)}, error = function(e) {
          cat(paste0("No rule was found for ", element, ". \n\n\n\n"))
          return(NULL)
        })
        if (is.null(result)) {next}
        
        for (i in 2:5){
          tryCatch({
            plot_list <- plot_func(i, element)
            
            if (i == 2) {panel <- "'Viewer' panel"}
            else {panel <- "'Plots' panel"}
            
            for (j in length(plot_list)) {
              plot.new()
              print(plot_list[[j]])  
              cat(paste0("Plot ", j, " for ", element, " ", i, "-node rules shown. Press [Enter] to continue...\n(The plot is shown in the ", panel, ")\n"))
              readline()
            }}, error = function(e) {
              cat(paste0("There are no ", i, "-node rules for ", element, ". Press [Enter] to continue...\n"))
              readline()
            })
        }}}}
}


# 3. Set your parameters

#**data** — dataset for analysis\
#**m1** - ARM metric for the nodes, can be any metric except for addedValue (support, confidence, coverage, lift, cosine)\
#**m2** - ARM metric for the edges, can be any symmetric metric (support,lift, cosine)\
#**m3** - ARM metric for the edges, can be any asymmetric metric (confidence, coverage, addedValue)\
#**unitCols**, **codeCols**, **conversationCols** —parameters for rENA model\
#**rotation** — the two groups whose differences you want to directly represent on the x-axis, must be a vector with two values, and these values must belong to the group column\
# 'rotation' continued: Select the pair of units/groups to apply means rotation
#**rotationCol** — Select the column with the unit/group you want to visualize\
#**MovStanzaSize** — moving window size for ARM\
#**threshold** — the cut-off to decide whether the rules should be shown (based on the symmetric metric)\
#**colors** — must be a vector whose length is the number of groups (default: all red)\
#**nsize** — coefficient to adjust the size of nodes (default: 1)\
#**esize** — coefficient to adjust the size of edges (default: 1)\
#**xlim** — a vector controlling the upper and lower limit of x-axis (default: NULL) **ylim** — a vector controlling the upper and lower
#limit of y-axis (default: NULL)\
#**dir** — directory for saving the csv files for all the rules (default:NULL)

# Following is an example with a dataset in rENA package:
  
## Dataset #1: data from the rENA package
# data = rENA::RS.data  
# m1 = 'support'
# m2 = 'support'
# m3 = 'confidence'
# unitCols = c("Condition", "UserName")
# codeCols = c('Data', 'Technical.Constraints', 'Performance.Parameters', 'Client.and.Consultant.Requests', 'Design.Reasoning', 'Collaboration')
# conversationCols = c("Condition", "GroupName", "ActivityNumber")
# rotation <- c("FirstGame", "SecondGame")
# rotationCol <- "Condition"
# MovStanzaSize <- 8
# threshold <- 0.3
# esize <- 0.3
# threshold_diff <- 0.1
# esize_diff <- 0.1
# xlim <- c(-2.8, 2.5)
# ylim <- c(-2.5, 3.2)
# dir <- "~/Desktop/Research/PCLA/EARN/"


## Dataset #2: hamletRomeoJuliet data 
setwd("~/Desktop/Research/PCLA/ENA/")
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
 data <- read.csv("hamletRomeoJuliet.csv") 
 m1 = 'cosine'
 m2 = 'cosine'
 m3 = 'addedValue'
 unitCols = c("Play", "Act", "Speaker")
 codeCols = c('Love', 'Beauty', 'Death', 'Fear', 'Friendship', 'Hate', "Honor", "Men", "Women", "Pride")
 conversationCols = c("Play", "Act", "Scene")
 data$Play <- make.names(data$Play)
 data$Act <- paste0("Act_", data$Act) # to make sure the levels are converted to text labels
 rotation <- c("Act_1", "Act_3", "Act_2")
 rotationCol <- "Act" # can also try the Play column with "Hamlet" and "Romeo.and.Juliet", just make sure to run this code which is a couple lines above: data$Play <- make.names(data$Play)
 MovStanzaSize <- 4
 threshold <- 0.5
 esize <- 0.4
 threshold_diff <- 0.01
 esize_diff <- 0.01
 xlim <- c(-2.8, 2.5)
 ylim <- c(-2.5, 3.2)
 dir <- "~/Desktop/Research/PCLA/EARN/"


# 4. Run the functions
ena_arm_func(data, m1, m2, m3, unitCols, codeCols, conversationCols, rotation=rotation, rotationCol=rotationCol, MovStanzaSize=MovStanzaSize, threshold=threshold, esize=esize)
ena_arm_func_diff(data, m1, m2, m3, unitCols, codeCols, conversationCols, rotation, rotationCol, MovStanzaSize, threshold_diff, esize=esize_diff)



