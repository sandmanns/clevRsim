server <- function(input, output, session) {
  ##supporting functions copied from clevRvis (invisible in package)
  .getRelatedColors <- function(seaObject) {
    parents <- NULL
    
    lenBranch <- max(seaObject@nestLevels) + 1
    colorMatrix <- .getColorMatrix(lenBranch = lenBranch)
    nclones <- length(seaObject@cloneLabels)
    cloneCols <- rep(NA, nclones)
    parentClone <- -1
    parentX <- c(1, rep(NA, nclones))
    parentY <- c(0, rep(NA, nclones))
    cloneBranches <- .getAllBranches(seaObject@parents)
    
    ## first stablish the color for the initial clones (parent = 0)
    indepClones <- which(seaObject@parents == 0)
    if (length(indepClones) > 1) {
      step <-
        25 %/% length(indepClones) #as far away as possible in the matrix
      x <- 0
      for (clone in indepClones) {
        x <- x + step
        y <- 1
        if (x > 25) {
          x <- x %/% 25
        }
        cloneCols[clone] <- colorMatrix[y, x]
        parentX[clone + 1] <- x
        parentY[clone + 1] <- y
      }
      ##get which clones where not initial
      branchedClones <-
        seq_len(nclones)[-which(seq_len(nclones) %in% indepClones)]
    } else{
      branchedClones <- seq_len(nclones)
    }
    
    ##Deal with the non-initial clones
    for (clone in branchedClones) {
      parent <- seaObject@parents[clone]
      ##color change depends on the number of branches
      if (sum(seaObject@parents == parent) > 4) {
        branchStep <- 0
      } else {
        branchStep <- 2
      }
      ## If it's a branch (parent already in the vector)
      if (parent %in% parentClone) {
        counts <- sum(parentClone == parent) + branchStep
        modX <- FALSE
        ## alternate between the two sides of the matrix to variate the
        ## color, darken the color
        if (counts %% 2 == 1) {
          x <- parentX[parent + 1] + ((counts + 1) / 2)
          
          if (x %in% parentX) {
            if (sum((which(parentX == x) - 1) %in% cloneBranches[[clone]]) == 0) {
              parentsInBranch <-
                which(names(table(parents)) %in% cloneBranches[[clone]])
              counts <- counts +
                sum(table(parents)[parentsInBranch] > 1)
              if (counts %% 2 == 1) {
                x <- parentX[parent + 1] +
                  ((counts + 1) / 2) + branchStep
              } else{
                x <- parentX[parent + 1] - (counts / 2) -
                  branchStep
              }
              modX <- TRUE
            }
          }
        } else if (counts %% 2 == 0 & !modX) {
          x <- parentX[parent + 1] - (counts / 2)
          if (x %in% parentX) {
            if (sum((which(parentX == x) - 1) %in% cloneBranches[[clone]]) == 0) {
              parentsInBranch <-
                which(names(table(parents)) %in% cloneBranches[[clone]])
              counts <-
                counts + sum(table(parents)[parentsInBranch] > 1)
              x <- parentX[parent + 1] + (counts / 2)
              modX <- TRUE
            }
          }
        }
        ##For linear evolution (parent not in the vector) --> darken the
        ##color
      } else {
        x <- parentX[parent + 1]
      }
      y <- parentY[parent + 1] + 1
      
      ##Fix x and y values to be inside the matrix
      if (x > 25) {
        x <- x %/% 25
      }
      if (y > lenBranch) {
        y <- y %/% lenBranch
      }
      if (x == 0) {
        x <- 25
      }
      if (x < 0) {
        x <- 25 + x
      }
      
      
      ##check another clone doesn't have the same color
      repColor <- colorMatrix[y, x] %in% cloneCols
      if (repColor) {
        while (repColor) {
          if (x < 12) {
            x <- x + 1 + branchStep
          } else {
            x <- x - 1 - branchStep
          }
          
          repColor <- colorMatrix[y, x] %in% cloneCols
        }
      }
      
      cloneCols[clone] <- colorMatrix[y, x] ##Add color to the vector
      parentClone <-
        c(parentClone, seaObject@parents[clone]) #update parents vec
      parentX[clone + 1] <-
        x #update parents vectors with the clone x and y positions
      parentY[clone + 1] <- y
    }
    return(cloneCols)
  }
  
  .getColorMatrix <- function(lenBranch) {
    ordr <- c(seq(18, 25), seq(1, 17))
    rbw <- rainbow_hcl(25, l = 85)[ordr]
    allColors <- c()
    for (col in rbw) {
      darkCol <- darken(col, 0.7)
      darkColors <- colorRampPalette(c(col, darkCol))(lenBranch)
      allColors <- c(allColors, darkColors)
    }
    colMatrix <- matrix(allColors, ncol = length(rbw), byrow = FALSE)
    
    return(colMatrix)
  }
  
  .getAllBranches <- function(parents) {
    branch <- list()
    for (i in seq_len(length(parents))) {
      branch[[i]] <- c(.getBranch(parents, i), i)
    }
    return(branch)
  }
  
  .getBranch <- function(parents, x) {
    #sanity checks
    if (x > length(parents)) {
      stop(
        paste(
          "cannot have a parent that does not exist in list. parent =",
          x,
          ", length(parents) =",
          length(parents)
        )
      )
    }
    if (x < 0) {
      stop("cannot have a value in parents of less than zero")
    }
    
    if (parents[x] == 0) {
      return(c(0)) #return initial clone (normal cell)
    } else {
      ##return clone's parent and compute again
      return(c(.getBranch(parents, parents[x]), parents[x]))
    }
  }
  
    # Initialization #####
    # Hide all outputs
    hide_outputs <- function(){
        output$fishPlot <- renderPlot({plot.new()})
        hide("clonetable")
        hide("mutationtable")
        hide("updatetext1")
        hide("updatetext2")
        hide("updatebtn1")
        hide("updatebtn2")
        hide("uiaddclone")
        hide("uiaddsample")
        hide("uideletesamplevalue")
        hide("uideleteclonevalue")
        hide("uideleteclone")
        hide("uideletesample")
        hide("addcnvbtn")
        hide("deletecnvvalue")
        hide("deletecnvbtn")
        hide("cnvtable")
        hide("puritytable")
        output$button <- renderUI({h5("Choose Basic Settings first and simulate")})
    }
    #Hide outputs at beginning
    #hide_outputs()
    
    #Empty mutation table that is diplayed before simulation
    output$mutationtable <- renderDataTable((data.frame(ID = 0, Chr = 0, Start = 0, End = 0, Ref = 0, Var = 0, VAF = 0)))
    #show("mutationtable")
    

    output$independentpng <- renderImage({
        
        # Return a list containing the filename
        return(list(
            src = "images/independent.png",
            contentType = 'image/png',
            height = 120,
            width = 200,
            alt = "Logo of IMI"))
    }, deleteFile = FALSE)
    
    output$linearpng <- renderImage({
        
        # Return a list containing the filename
        return(list(
            src = "images/linear.png",
            contentType = 'image/png',
            height = 120,
            width = 200,
            alt = "Logo of IMI"))
    }, deleteFile = FALSE)
    
    output$parallelpng <- renderImage({
        
        # Return a list containing the filename
        return(list(
            src = "images/parallel.png",
            contentType = 'image/png',
            height = 120,
            width = 200,
            alt = "Logo of IMI"))
    }, deleteFile = FALSE)
    
    
    # Update mutationcount if clones is modified so that Variants >= Clones
    observe({
        req(input$clonecount)
        clon <- as.numeric(input$clonecount)
        val <- as.numeric(input$mutationcount)
        if(val < clon){
            val = ceiling(clon*2)
        }
        updateNumericInput(session, "mutationcount", label = "#Variants (SNVs+CNVs)", min= clon, max=200, value = val)
    })
    
    # Update CNV slider so that max CNV == max Variants
    observe({
        req(input$mutationcount)
        mut <- as.numeric(input$mutationcount)
        val <- as.numeric(input$cnvcount)
        if(val > mut){
            val = 0
        }
        updateSliderInput(session, "cnvcount", label = "#CNVs (included in #Variants)", min = 0, max = mut, value = val)
    })
    
    # Unselect Checkboxes if CNV = 0 
    observe({
        req(input$cnvcount)
        cnv <- input$cnvcount
        if(cnv == 0){
            updateCheckboxGroupInput(session, "cnvcheckbox", "(opt.) Specification", choices = c("Duplication","Deletion","LOH"), 
                                     selected = NULL,inline = T)
        }
    })
    
    inputvalues <- reactiveValues(data = NULL)
    
    #Control variables which influence outputs (bad for fishplot)
    rsamples <- eventReactive(input$simulatebtn, {input$samplevalue})
    rclones <- eventReactive(input$simulatebtn, {input$clonecount})
    rclonetable <- eventReactive(input$updatebtn1, {input$clonetable})
    rclonetable <- eventReactive(input$updatebtn2, {input$clonetable})
    
    # Simulation ####
    observeEvent(input$simulatebtn, {
        if(input$evolutionvalue == "branched dependent" && input$clonecount <= 2){
            #hide_outputs()
            showNotification("Please select at least 3 clones for branched dependant evolution.", duration = 5, closeButton = TRUE, type = "error")
        } else {
            req(input$mutationcount)
            req(input$clonecount)
            id <- showNotification("Simulation is running. This can take up to 2 minutes.", duration = NULL, closeButton = FALSE, type = "message")
            # Try to create a phylogeny with all chosen basic settings
            phylogeny <<- create_phylogeny(input$mutationcount, input$clonecount, input$samplevalue, input$evolutionvalue, input$detectionvalue, 0, input$distancevalue)
            
            ##check for clevRvis rules
            if(!is.null(phylogeny)){
              clones <- phylogeny[-1,-1]
              samplecount = ncol(clones)-2
              clonecount = nrow(clones)
              fracTable = as.matrix(clones[,-(1:2)]) #P arent+Variants cutted
              timepoints <- 1:samplecount
              parents <- as.vector(clones[,1])
              cloneLabels<-paste0("clone ",1:clonecount)
              count_clevRvis<-0
              
              abc=data.frame(x=1)
              data<-lapply(abc,function(x){
                tryCatch({
                  seaObject <- createSeaObject(fracTable, parents, timepoints = timepoints, 
                                               cloneLabels=cloneLabels, timepointInterpolation = T,originTimepoint = 0)
                },error=function(e) NULL)
              })
              while(is.null(data[[1]])&&count_clevRvis<5){
                phylogeny <<- create_phylogeny(input$mutationcount, input$clonecount, input$samplevalue, input$evolutionvalue, input$detectionvalue, 0, input$distancevalue)
                clones <- phylogeny[-1,-1]
                samplecount = ncol(clones)-2
                clonecount = nrow(clones)
                fracTable = as.matrix(clones[,-(1:2)]) #P arent+Variants cutted
                timepoints <- 1:samplecount
                parents <- as.vector(clones[,1])
                cloneLabels<-paste0("clone ",1:clonecount)
                count_clevRvis<-count_clevRvis+1
                
                abc=data.frame(x=1)
                data<-lapply(abc,function(x){
                  tryCatch({
                    seaObject <- createSeaObject(fracTable, parents, timepoints = timepoints, 
                                                 cloneLabels=cloneLabels, timepointInterpolation = T,originTimepoint = 0)
                  },error=function(e) NULL)
                })
              }
            }
            
            if(is.null(phylogeny)){# no phylogeny found
                #hide_outputs()
                removeNotification(id)
                showNotification("No resolution found. Try fewer clones or more samples.", duration = 5, closeButton = TRUE, type = "error")
                clones <<- NULL
                mutations <<-NULL
                cnvs <<- NULL
                purities <<-NULL
                showcnvs <<- NULL
            }else{# phylogeny successfully created
                showcnvs <<- NULL
                clones <<- phylogeny[-1,-1]# delete simulation artifacts (Number (column), Ancestor (row))

                #Create CNVs
                if (input$cnvcount > 0 ){
                    cnvs <<- create_cnvs(input$cnvcheckbox, input$cnvcount, clones)
                    showcnvs <<- create_showcnvs(cnvs)
                } else{# no CNVs chosen
                    cnvs <<- NULL
                    showcnvs <<- NULL
                }
                #Initial Puritytable
                purities <<- data.frame(matrix(data = as.numeric(input$purityvalue), nrow = 1, ncol = as.numeric(rsamples())))
                colnames(purities) <<- paste0("t", 1:rsamples())

                # Create SNVs
                mutations <<- create_snvs(clones, input$coveragevalue, purities, 
                                          input$distancevalue, cnvs)[[1]]
                removeNotification(id)
                # Update all tables and plots
                
                update_outputs()
            }
        }
    })
    
    
    # Modularized Updatefunction (Called by Update- & Delete Buttons) ####
    update  <- function(pclones, pcnvs = NULL){
        #hide("updatetext1") # Delete old messages.
        #hide("updatetext2")
        u_clones = pclones # needed parameter
        
        if(any(is.na(u_clones))){
            output$updatetext1 <- renderText("Please fill all cells (Fine-tuning: Clones).")
            output$updatetext2 <- renderText("Please fill all cells (Fine-tuning: Clones).")
            #show("updatetext1")
            #show("updatetext2")
            return()
        }
        if(!is.null(input$puritytable)){
            u_purities <- round(hot_to_r(input$puritytable,0))
            if(any(is.na(u_purities))){
                output$updatetext1 <- renderText("Please fill all cells (Fine-Tuning: Purity).")
                output$updatetext2 <- renderText("Please fill all cells (Fine-Tuning: Purity).")
                #show("updatetext1")
                #show("updatetext2")
                return()
            }
        }else{
            u_purities <- NULL
        }
        warning <-NULL
        if(!is.null(input$cnvtable)){
            if(is.null(pcnvs)){
                u_cnvs <- hot_to_r(input$cnvtable)
            } else {
                u_cnvs = pcnvs #If a CNV was deleted, the cnv table is given as parameter for checks before deletion.
                if(nrow(u_cnvs) < 1){#Last CNV was deleted
                    u_cnvs <- NULL
                    u_showcnvs <- NULL
                }
            }
            if(nrow(u_cnvs) > 0){
                if(any(is.na(u_cnvs[,-8]))){
                    output$updatetext1 <- renderText("Please fill all cells (Fine-Tuning: CNVs).")
                    output$updatetext2 <- renderText("Please fill all cells (Fine-Tuning: CNVs).")
                    #show("updatetext1")
                    #show("updatetext2")
                    return()
                }
                # Check CNVs for validity
              
                cnvcheck <- check_cnvs(u_cnvs, u_clones[,"#Variants"], (sum(u_clones[,"#Variants"])-nrow(u_cnvs)), u_clones[,-(1:2)])
                if(!cnvcheck[[1]]){# Something went wrong. Errormessage in cnvcheck[[2]]
                    output$updatetext1 <- renderText(cnvcheck[[2]])
                    output$updatetext2 <- renderText(cnvcheck[[2]])
                    #show("updatetext1")
                    #show("updatetext2")
                    return() # no update of intern variables happened
                } else {
                    # CNVs are checked successfully
                    warning <- cnvcheck[[2]]
                    u_cnvs <- cnvcheck[[3]]
                    # Create CNV input for mutation table
                    u_showcnvs <- create_showcnvs(u_cnvs)
                }
            } else{
                u_cnvs <- NULL
                u_showcnvs <- NULL
            }
        }else{
            u_cnvs <- NULL
            u_showcnvs <- NULL
        }
        # CNV and purity are checked. Now Phylogeny table.
        clonecheck <- check_phylogeny(u_clones, input$detectionvalue, input$mutationcount, input$distancevalue, u_cnvs)
        if(clonecheck[[1]]){# CLONES VALID
            clones <<- as.data.frame(clonecheck[[3]])
            if(!is.null(u_purities)){
                purities <<- u_purities# update intern variables
            }
            if(!is.null(u_cnvs)){
                cnvs <<- u_cnvs# update intern variables
            }
            if(!is.null(u_showcnvs)){
                showcnvs <<- u_showcnvs# update intern variables
            }
            id <- showNotification("Updating (creating new variant data may take some seconds).", duration = NULL, closeButton = FALSE, type = "message")
            # Create new SNV data with updated clones, puritites, cnvs
            mutationresult <- create_snvs(clones, input$coveragevalue, purities, 
                                          input$distancevalue, cnvs)
            removeNotification(id)
            mutations <<- mutationresult[[1]] # update intern variables
            warning <- c(warning, mutationresult[[2]])
            
        }
        update_outputs()
        
        
        output$updatetext1 <- renderText(paste(warning,clonecheck[[2]])) # Show Warning or success message
        output$updatetext2 <- renderText(paste(warning,clonecheck[[2]]))
        #show("updatetext1")
        #show("updatetext2")
        #print("UPDATE COMPLETE!") # intern print
    }
    
    # Update-Button ####
    observeEvent({
        input$updatebtn1 
        input$updatebtn2
    }, {
        if(!is.null(input$updatebtn1) || !is.null(input$updatebtn2)){
            if(!is.null(input$clonetable)){
                u_clones = round(hot_to_r(input$clonetable),2)
                colnames(u_clones)[2]<-"#Variants"
                update(u_clones) # modularized update function
                ######################################Fehler abfangen?! bei update-funktion
                
            } #else { 
               # print("No Clonetable found") # intern print
            #}
        }
    }, ignoreNULL = FALSE)
    
    
    # Add Clone ####
    observeEvent(input$addclonebtn,{
        phylodata <- round(hot_to_r(input$clonetable),2)
        new_row <-  data.frame(matrix(c(0,0,1, rep(0, times = ncol(phylodata)-3)), nrow = 1))
        colnames(new_row) <- colnames(phylodata)
        rownames(new_row) <- nrow(phylodata)+1
        DF <- rbind(phylodata, new_row)
        update_outputs(DF)
        output$updatetext1 <- renderText("Empty Clone added. Please fill in values.")
        #show("updatetext")
    })
    
    
    # Add Sample ####
    observeEvent(input$addsamplebtn,{
        phylodata <- round(hot_to_r(input$clonetable),2)
        new_column <- data.frame(matrix(rep(0, times = nrow(phylodata)), ncol = 1))
        colnames(new_column) <- paste0("t",ncol(phylodata)-1)
        rownames(new_column) <- rownames(phylodata)
        clones <<- cbind(phylodata, new_column) 
        purities <<- cbind(purities, as.numeric(input$purityvalue))
        colnames(purities) <<- paste0("t",1:ncol(purities))
        update_outputs()
        output$updatetext1 <- renderText("Empty Sample added. Please fill in values.")
        #show("updatetext1")
    })
    
    
    # Add CNV ####
    observeEvent(input$addcnvbtn,{
        cnvdata <- hot_to_r(input$cnvtable)
        clonedata <- hot_to_r(input$clonetable)
        samplecount <- ncol(clonedata)-2
        samples <- matrix(rep(0, times = samplecount), nrow = 1)
        names <- paste0("t", 1:samplecount)
        colnames(samples) <- names
        new_row <-  cbind(data.frame(ID = paste("CNV", nrow(cnvdata)+1, sep="_"), Chr = 1, Start = 1, End = 1, Type = NA,  Overlap = FALSE, "SNVs" = 0, Scenario = NA, in_Clone = 0), samples)
        rownames(new_row) <- nrow(cnvdata)+1
        DF <- rbind(cnvdata, new_row)
        update_outputs(pcnvs = DF)
        output$updatetext2 <- renderText("Empty CNV added. Please fill in values.")
        #show("updatetext2")
    })
    
    
    # Delete Sample ####
    observeEvent(input$deletesamplebtn,{
        samplenumber <- as.numeric(substr(input$deletesamplevalue,2,nchar(input$deletesamplevalue)))
        maxsamples <- (ncol(hot_to_r(input$clonetable))-2)
        if(maxsamples < 2){
            output$updatetext1 <- renderText(paste("Last sample can not be deleted."))
            return()
        }
        if(samplenumber < 1 || samplenumber > maxsamples){
            range <- paste0(1,"-",maxsamples)
            output$updatetext1 <- renderText(paste("Invalid Sample Number. Range is", range))
        }else{
            DF <- round(hot_to_r(input$clonetable),2)
            DF <- DF[,-(samplenumber+2)]
            # Check if phylogeny would still be valid after deleting clone
            clonecheck <-  check_phylogeny(DF, input$detectionvalue, input$mutationcount, input$distancevalue, u_cnvs)
            if(clonecheck[[1]]){# Clone can be deleted
                clones <<- as.data.frame(clonecheck[[3]]) # update intern variables
                colnames(clones) <<- c("Parent", "#Variants", paste0("t",1:(ncol(clones)-2)))
                purities <<- data.frame(purities[,-(samplenumber)]) # update intern variables
                if(length(purities)>1){
                    colnames(purities) <<- paste0("t",1:ncol(purities))
                }else{
                    pur <- cbind(purities, 0)
                    colnames(pur) <- c("t1","t2")
                    purities <<- pur[,1]
                }
                cnvs <<- cnvs[,-(samplenumber+9)] # update intern variables
                colnames(cnvs) <<- c("ID", "Chr", "Start", "End", "Type", "Overlap", "SNVs",
                                     "Scenario", "in_Clone", paste0("t",1:(ncol(clones)-2)))
                showcnvs <<- create_showcnvs(cnvs) # update intern variables
                mutations <<- create_snvs(clones, input$coveragevalue, purities, 
                                          input$distancevalue, cnvs)[[1]] # update intern variables
                update_outputs() # update all tables and plots
                output$updatetext1 <- renderText("Sample deleted successfully.")
            } else { # Clone can not be deleted. Reason in clonecheck[[2]]
                output$updatetext1 <- renderText(clonecheck[[2]])
            }
            #show("updatetext1")
        }
    })
    
    
    # Delete Clone ####
    observeEvent(input$deleteclonebtn,{
        clonenumber <- as.numeric(input$deleteclonevalue)
        phylodata <- round(hot_to_r(input$clonetable),2)
        maxclones <- nrow(phylodata)
        if(maxclones < 2){
            output$updatetext1 <- renderText(paste("Last Clone can not be deleted."))
            return()
        }
        if(clonenumber < 1 || clonenumber > maxclones){
            range <- paste0(1,"-",maxclones)
            output$updatetext1 <- renderText(paste("Invalid clone number. Range is", range))
        }else{
            DF <- delete_clone(phylodata, clonenumber) # Helper Function in CheckTable.R
            update(DF) # updated with DF corrected by helper function
        }
    })
    
    # Delete CNV ####
    observeEvent(input$deletecnvbtn,{
        cnvnumber <- as.numeric(input$deletecnvvalue)
        maxcnv <- nrow(hot_to_r(input$cnvtable))
        if(maxcnv == 0){
            output$updatetext2 <- renderText(paste("No CNV to be deleted."))
        } else if(cnvnumber < 1 || cnvnumber > maxcnv){
            range <- paste0(1,"-",maxcnv)
            output$updatetext2 <- renderText(paste("Invalid CNV Number. Range is", range))
        }else{
            DF <- hot_to_r(input$cnvtable)
            u_cnvs <- DF[-cnvnumber,]
            if(nrow(u_cnvs)>0){#New IDs for ongoing numeration
                rownames(u_cnvs) <- 1:nrow(u_cnvs)
                u_cnvs$ID <- paste("CNV", 1:nrow(u_cnvs), sep="_")
            } 
            DF <- hot_to_r(input$clonetable)
            update(DF, u_cnvs)
        }
    })
    
    # Round values ####
    observeEvent(rclonetable(), { # 2 digits clonefreqs        
      clones <<- round(hot_to_r(rclonetable()),2)
        draw_colortable2(clones)
    })
    
    #observeEvent(rclonetable(), { # 2 digits clonefreqs        
    #  clones <<- round(hot_to_r(input$clonetable),2)
    #  draw_colortable2(clones)
    #})
    
    
    
    observeEvent(input$puritytable, { # 0 digits purity 
        purities <<- round(hot_to_r(input$puritytable),0)
        output$puritytable <- renderRHandsontable({
            rhandsontable(data.frame(purities),rowHeaderWidth = 75, rowHeaders = NULL,
                          colHeaders = colnames(purities)) %>%
                hot_table(highlightCol = TRUE, highlightRow = TRUE,  contextMenu = FALSE ) %>%
                hot_cols(format = "0", halign = "htCenter") %>%
                hot_validate_numeric(cols = 1:ncol(purities), min = 1, max = 100) 
        })
    })
    
    # Load session ####
    observeEvent(input$upload, { 
        infile <- input$upload
        if(is.null(infile)){
            return(NULL)
        } else {
            #hide_outputs()
            numfiles <- nrow(infile)
        }
        if(numfiles != 6){
            output$updatetext1 <- renderText("Please choose exact 6 RDS files (inputs, clones, variants, cnvs, showcnvs, purities).")
            output$updatetext2 <- renderText("Please choose exact 6 RDS files (inputs, clones, variants, cnvs, showcnvs, purities).")
            #show("updatetext1")
            #show("updatetext2")
            return()
        }
        # All 6 files loaded
        inputs <- NULL
        for(i in 1:numfiles){
            switch (infile$name[i],
                    "CCF-matrix.RDS" =  clones <<- readRDS(infile$datapath[i]),
                    "SNVs.RDS" = mutations <<- readRDS(infile$datapath[i]),
                    "CNVs.RDS" =  cnvs <<- readRDS(infile$datapath[i]),
                    "CNVs_VariantCalls.RDS" = showcnvs <<- readRDS(infile$datapath[i]),
                    "Purity.RDS" =  purities <<- readRDS(infile$datapath[i]),
                    "Input_Shiny.RDS" = inputs <- i,
                    # 6 files chosen, but not every file is correct
                    {  output$updatetext1 <- renderText(paste("Unknown file:", infile$name[i],". Please choose the following 6 RDS files: inputs, clones, variants, cnvs, showcnvs, purities."))
                    output$updatetext2 <- renderText(paste("Unknown file:", infile$name[i],". Please choose the following 6 RDS files: inputs, clones, variants, cnvs, showcnvs, purities."))
                    #show("updatetext1")
                    #show("updatetext2")
                    return()
                    }
            )
        }
        savedInputs <- readRDS(infile$datapath[inputs])
        inputvalues <- unlist(savedInputs)
        inputIDs <- names(inputvalues)
        # restore all setting values
        for (i in 1:length(inputvalues)) { 
            session$sendInputMessage(inputIDs[i],  list(value=inputvalues[[i]]) )
        }
        update_outputs() #update all table and plot outputs
    })
    
    
    # Update Outputs (Visualization)####
    # All intern variables have to be uptaded before. This function only visualizes all outputs.
    update_outputs <- function(pclones = clones, ppurities = purities, pcnvs = cnvs, pmutations = mutations, pshowcnvs = showcnvs){ 
        output$updatetext1 <- renderText("")
        output$updatetext2 <- renderText("")

        # Clonetable
        draw_colortable2(pclones)
        #show("clonetable")
        
        # Mutationtable
        if(is.null(pshowcnvs)){ # no CNVs selected
            allmutations <<- rbind(pmutations, pshowcnvs) 
        }else{ # CNVs present
            allmutations <<- rbind(pmutations,setNames(pshowcnvs, names(pmutations)))
        }
        rownames(allmutations) <- 1:nrow(allmutations)
        output$mutationtable <-renderDataTable(allmutations, options = list(
            pageLength = 10,
            autoWidth = TRUE,
            columnDefs = list(list(width = '10px', targets = c(1:5)))
        ))
        #show("mutationtable")
        
        # Puritytable
        output$uipuritytitle <- renderUI({
            h4("Fine-tuning: Purity")
        })
        if(length(ppurities)>1){
            if(ncol(ppurities)>ncol(pclones)-2){#Puritytable exploited, shorten
                ppurities <<- ppurities[,1:(ncol(pclones)-2)]
                purities <<- ppurities[,1:(ncol(pclones)-2)]
            } 
            cols <- 1:ncol(purities)
        }else{
            cols <- 1
        }
        t1 <- ppurities
        output$puritytable <- renderRHandsontable({
            rhandsontable(data.frame(t1), rowHeaderWidth = 75, rowHeaders = NULL) %>%
                hot_table(highlightCol = TRUE, highlightRow = TRUE,  contextMenu = FALSE ) %>%
                hot_cols(format = "0", halign = "htCenter") %>%
                hot_validate_numeric(cols = cols, min = 1, max = 100) %>%
                htmlwidgets::onRender("
                   function(el, x) {
                    var hot = this.hot
                    $('a[data-value=\"Output\"').click(function() {
                      setTimeout(function() {hot.render();}, 0);
                     })
          }")
        })
        
        #show("puritytable")
        outputOptions(output, "puritytable", suspendWhenHidden = FALSE) # Needed to create table without changing into Data tab
        
        # CNVtable
        if(is.null(pcnvs)){ # no CNVs present, provide an empty table
            pcnvs <- data.frame(ID = NA, Chr = NA, Start = NA, End = NA, Type = NA,  Overlap = FALSE, "SNVs" = 0, Scenario = NA, in_Clone = 0)
            pcnvs <- pcnvs[-1,]
        }
        
        output$cnvtable <- renderRHandsontable({
            rhandsontable(data.frame(pcnvs),height=(length(pcnvs[,1])*30+100)) %>%
                hot_table(highlightCol = TRUE, highlightRow = TRUE,  contextMenu = FALSE ) %>%
                hot_cols(format = "0", halign = "htCenter", allowInvalid = FALSE, readOnly = TRUE) %>%
                hot_col(col = "Type", type = "dropdown", source = c("Duplication", "Deletion", "LOH"), strict = TRUE,readOnly = FALSE)%>%
                hot_col(col = "Scenario", type = "dropdown", source = c("CNV first","SNV first (affected)","SNV first (un-affected)", "Parallel"), strict = TRUE, readOnly = FALSE)%>%
                hot_col(col = c("in_Clone","Overlap","SNVs"), readOnly = FALSE) %>%
                hot_col(col = c("Chr","End","Start"), readOnly = FALSE) %>%
                hot_validate_numeric(cols = 2, min = 1, max = 23) %>%
                hot_validate_numeric(cols = 3, min = 1) %>%
                hot_validate_numeric(cols = 4, min = 1) %>%
                hot_col(col = "ID", readOnly = TRUE) %>%
                htmlwidgets::onRender("
                   function(el, x) {
                    var hot = this.hot
                    $('a[data-value=\"Output\"').click(function() {
                      setTimeout(function() {hot.render();}, 0);
                     })
          }")
        })
        
        #show("addcnvbtn")
        #show("deletecnvvalue")
        #show("deletecnvbtn")
        #show("uicnvtitle")
        #show("cnvtable")
        
        outputOptions(output, "cnvtable", suspendWhenHidden = FALSE) # Needed to create table without changing into Data tab
        
        #Fishplot
        output$fishPlot <- renderPlot({
            draw_fishplot2(pclones) # helper function
        })
        output$fishPlot2 <- renderPlot({
          draw_fishplot2b(pclones) # helper function
        })
        #show("fishPlot")
        #show("fishPlot2")
        
        # Add all buttons to the screen (Not seen before first simulation)
        output$uiheader <- renderUI({
          h4("Fine-tuning: Clones")
        })
        output$uiaddclone <- renderUI({
            actionButton("addclonebtn","Add clone")
        })
        output$uiaddsample <- renderUI({
            actionButton("addsamplebtn","Add time point")
        })
        output$uideletesamplevalue <- renderUI({
            #numericInput("deletesamplevalue", "Choose time point to delete", value = 1, min = 1, width = 500)
          selectInput("deletesamplevalue", label=NULL,choices = paste0("t",1:(length(pclones[1,])-2)),selected = "t1", width = "80%")
          
        })
        output$uiaddsamplevalue <- renderUI({
          selectInput("addsamplevalue", label=NULL,choices = paste0("t",length(pclones[1,])-1),selected = paste0("t",ncol(purities)+1), width = "80%")
        })
        
        output$uideletesample <- renderUI({
            actionButton("deletesamplebtn", "Delete time point")
        })

        #output$uideleteclonevalue <- renderUI({
        #    selectInput("deleteclonevalue", label=NULL, selected = 1, choices = 1:nrow(round(hot_to_r(input$clonetable),2)), width = "80%")
        #})
        #output$uiaddclonevalue <- renderUI({
        #  selectInput("addclonevalue", label=NULL, selected = (nrow(round(hot_to_r(input$clonetable),2))+1), 
        #               choices = (nrow(round(hot_to_r(input$clonetable),2))+1), width = "80%")
        #})
        
        output$uideleteclonevalue <- renderUI({
          selectInput("deleteclonevalue", label=NULL, selected = 1, choices = 1:max(allmutations$Clone), width = "80%")
        })
        output$uiaddclonevalue <- renderUI({
          selectInput("addclonevalue", label=NULL, selected = (max(allmutations$Clone)+1), 
                      choices = (max(allmutations$Clone)+1), width = "80%")
        })
        
        output$uideleteclone <- renderUI({
            actionButton("deleteclonebtn", "Delete clone")
        })
        output$uicnvtitle <- renderUI({
            h4("Fine-tuning: CNVs")
        })
        output$uicnvtitle2 <- renderUI({
           h5("To simulate overlapping CNVs:")
        })
        output$uicnvtitle2a <- renderUI({
          h5(HTML(rep("&nbsp",4)),"* check 'Overlap'")
        })
        output$uicnvtitle2b <- renderUI({
          h5(HTML(rep("&nbsp",4)),"* define number of 'SNVs' to overlap (>0)")
        })
        output$uicnvtitle2c <- renderUI({
          h5(HTML(rep("&nbsp",4)),"* select a 'Scenario' of overlap (CNV first, SNV first (affected), SNV first (un-affected), Parallel)")
        })
        output$uicnvtitle2d <- renderUI({
          h5(HTML(rep("&nbsp",4)),"* make sure 'in_Clone' does not contain less SNVs than defined in 'SNVs' (otherwise: change 'in_Clone')")
        })
        output$uicnvtitle2e <- renderUI({
          h5(HTML(rep("&nbsp",4)),"* (opt.) change 'Chr', 'Start', 'End', 'Type'")
        })
        
        output$uicnvtitle_scenarios1 <- renderUI({
          h5("The influence of 'Type' and 'Scenario' on the SNP genotype (Normal: AB (AA->AB)):")
        })

        output$uicnvtitle_scenarios3a <- renderUI({
          h5("CNV first")
        })
        output$uicnvtitle_scenarios3b <- renderUI({
          h5(HTML(rep("&nbsp",2)),"Deletion: B (AA->A->B)")
        })
        output$uicnvtitle_scenarios3c <- renderUI({
          h5(HTML(rep("&nbsp",2)),"Duplication: AAB (AA->AAA->AAB)")
        })
        output$uicnvtitle_scenarios3d <- renderUI({
          h5(HTML(rep("&nbsp",2)),"LOH: AB (AA->AA(LOH)->AB)")
        })
        
        output$uicnvtitle_scenarios4a <- renderUI({
          h5("SNV first (affected)")
        })
        output$uicnvtitle_scenarios4b <- renderUI({
          h5(HTML(rep("&nbsp",2)),"Deletion: AB, A (AA->AB->A)")
        })
        output$uicnvtitle_scenarios4c <- renderUI({
          h5(HTML(rep("&nbsp",2)),"Duplication: AB, ABB (AA->AB->ABB)")
        })
        output$uicnvtitle_scenarios4d <- renderUI({
          h5(HTML(rep("&nbsp",2)),"LOH: AB, BB (AA->AB->BB)")
        })
        
        output$uicnvtitle_scenarios5a <- renderUI({
          h5("SNV first (un-affected)")
        })
        output$uicnvtitle_scenarios5b <- renderUI({
          h5(HTML(rep("&nbsp",2)),"Deletion: AB, B (AA->AB->B)")
        })
        output$uicnvtitle_scenarios5c <- renderUI({
          h5(HTML(rep("&nbsp",2)),"Duplication: AB, AAB (AA->AB->AAB)")
        })
        output$uicnvtitle_scenarios5d <- renderUI({
          h5(HTML(rep("&nbsp",2)),"LOH: AB, AA (AA->AB->AA)")
        })
        
        output$uicnvtitle_scenarios6a <- renderUI({
          h5("Parallel")
        })
        output$uicnvtitle_scenarios6b <- renderUI({
          h5(HTML(rep("&nbsp",2)),"Deletion: AB (AA->AB; AA->A)")
        })
        output$uicnvtitle_scenarios6c <- renderUI({
          h5(HTML(rep("&nbsp",2)),"Duplication: AB (AA->AB; AA->AAA)")
        })
        output$uicnvtitle_scenarios6d <- renderUI({
          h5(HTML(rep("&nbsp",2)),"LOH: AB (AA->AB; AA->AA(LOH)")
        })
        
        
        
        
        
         output$uiaddcnvvalue <- renderUI({
          selectInput("addcnvvalue", label=NULL,choices = (length(pcnvs[,1])+1),selected = (length(pcnvs[,1])+1), width = "100%")
        })
        output$uideletecnvvalue <- renderUI({
          selectInput("deletecnvvalue", label=NULL,choices = c(1:length(pcnvs[,1])),selected = 1, width = "100%")
        })
        

        output$uiexport2 <- renderUI({
            downloadButton("exportbtn1","Export data")
        })
        
        #show("updatebtn1")
        #show("updatebtn2")
        #show("uiaddclone")
        #show("uiaddsample")
        #show("uideletesamplevalue")
        #show("uideleteclonevalue")
        #show("uideleteclone")
        #show("uideletesample")
        
    }
    
    setColors<-function(color){
      rgb_col<-col2rgb(color)
      hsv_col<-rgb2hsv(rgb_col)[,1]
      if(hsv_col[3]>0.5){
        return("#000000")
      }else{
        return("#ffffff")
      }
    }
    
    # Draw Colortable ####
    draw_colortable <-function(table){
        colors <- grDevices::rainbow(nrow(table),s = 0.7)
        hot <- rhandsontable(table, ColorCode = colors, rowHeaderWidth = 75) %>% 
            hot_cols(renderer = color_renderer, format = "00.00", max = 100, min = 0, halign = "htCenter", allowInvalid = FALSE) %>%
            hot_table(highlightCol = TRUE, highlightRow = TRUE,  contextMenu = FALSE ) %>%
            hot_col("Parent", format = "0", halign = "htCenter") 
        output$clonetable <- renderRHandsontable({hot})
        #show("clonetable")
    }
    
    draw_colortable2 <-function(table){
      samplecount = ncol(table)-2
      clonecount = nrow(table)
      fracTable = as.matrix(table[,-(1:2)]) #P arent+Variants cutted
      timepoints <- 1:samplecount
      
      parents <- as.vector(table[,1])
      cloneLabels<-paste0("clone ",1:clonecount)
      
      
      abc=data.frame(x=1)
      data<-lapply(abc,function(x){
        tryCatch({
          seaObject <- createSeaObject(fracTable, parents, timepoints = timepoints, 
                                       cloneLabels=cloneLabels, timepointInterpolation = T,originTimepoint = 0)
        },error=function(e) NULL)
      })
      if(!is.null(data[[1]])){
        seaObject <- createSeaObject(fracTable, parents, timepoints, 
                                     cloneLabels=cloneLabels,timepointInterpolation = T,originTimepoint = 0)
        colors<-.getRelatedColors(seaObject)
        
        colors2 <- c()
        for(i in 1:length(colors)){
          colors2[i]<-setColors(colors[i])
        }
        
        hot <- rhandsontable(table, ColorCode = colors, rowHeaderWidth = 75, ColorCode2 = colors2) %>% 
          hot_cols(renderer = color_renderer, format = "00.00", max = 100, min = 0, halign = "htCenter", allowInvalid = FALSE) %>%
          hot_table(highlightCol = TRUE, highlightRow = TRUE,  contextMenu = FALSE ) %>%
          hot_col("Parent", format = "0", halign = "htCenter") 
        output$clonetable <- renderRHandsontable({hot})
        #show("clonetable")
      }else{
        output$updatetext1 <- renderText("Please check your CCF table (basic assumptions violated: 
        1) children-clones cannot exceed their parents,
        2) clones developing from normal cells cannot add up to more than 100%, 
        3) a clone cannot reappear, 
        4) a parent-clone being thoroughly replaced by its children-clone(s) cannot reappear.).")
      }
    }
    
    draw_fishplot2 <- function(pclones){
      samplecount = ncol(pclones)-2
      clonecount = nrow(pclones)
      fracTable = as.matrix(pclones[,-(1:2)]) #P arent+Variants cutted
      timepoints <- 1:samplecount
      vlab = paste0("t", 1:samplecount)
      vlines = timepoints
      
      parents <- as.vector(pclones[,1])
      cloneLabels<-paste0("clone ",1:clonecount)

      abc=data.frame(x=1)
      data<-lapply(abc,function(x){
        tryCatch({
          seaObject <- createSeaObject(fracTable, parents, timepoints = timepoints, 
                                       cloneLabels=cloneLabels, timepointInterpolation = T,originTimepoint = 0)
        },error=function(e) NULL)
      })
      if(!is.null(data[[1]])){
        seaObject <- createSeaObject(fracTable, parents, timepoints, 
                                     cloneLabels=cloneLabels,timepointInterpolation = T,originTimepoint = 0)
        dolphinPlot(seaObject, shape = 'spline', borderCol = 'black', vlines = timepoints,
                    vlab = vlab, vlabSize = 7, main = "", xlab = 'Timepoints',
                    ylab = '', pad.left = 0.02, showLegend = F,mainSize = 8,separateIndependentClones = T)
      
      }else{
        output$updatetext1 <- renderText("Please check your CCF table (basic assumptions violated: 
        1) children-clones cannot exceed their parents,
        2) clones developing from normal cells cannot add up to more than 100%, 
        3) a clone cannot reappear, 
        4) a parent-clone being thoroughly replaced by its children-clone(s) cannot reappear.).")
      }
    }
    
    draw_fishplot2b <- function(pclones){
      samplecount = ncol(pclones)-2
      clonecount = nrow(pclones)
      fracTable = as.matrix(pclones[,-(1:2)]) #P arent+Variants cutted
      timepoints <- 1:samplecount
      vlab = paste0("t", 1:samplecount)
      vlines = timepoints
      
      parents <- as.vector(pclones[,1])
      cloneLabels<-paste0("clone ",1:clonecount)
      
      abc=data.frame(x=1)
      data<-lapply(abc,function(x){
        tryCatch({
          seaObject <- createSeaObject(fracTable, parents, timepoints = timepoints, 
                                       cloneLabels=cloneLabels, timepointInterpolation = T,originTimepoint = 0)
        },error=function(e) NULL)
      })
      if(!is.null(data[[1]])){
        seaObject <- createSeaObject(fracTable, parents, timepoints = timepoints, 
                                     cloneLabels=cloneLabels, timepointInterpolation = T,originTimepoint = 0)
        extSharkPlot(seaObject, main = "", showLegend = T, timepoints = timepoints,
                     width = 16,interactivePlot = F) 
      }else{
        output$updatetext1 <- renderText("Please check your CCF table (basic assumptions violated: 
        1) children-clones cannot exceed their parents,
        2) clones developing from normal cells cannot add up to more than 100%, 
        3) a clone cannot reappear, 
        4) a parent-clone being thoroughly replaced by its children-clone(s) cannot reappear.).")
      }
    }
    
    
    # Data Export ####
    output$exportbtn1  <- downloadHandler(
        filename = function() {paste("Simulation-", Sys.time(), ".zip", sep="")}, # For unique downloadfiles
        content = function(file) {
            
            # EXCELFILE with all data
            wb <- createWorkbook()
            addWorksheet(wb, "Input settings")
            addWorksheet(wb, "Visualization")
            addWorksheet(wb, "CCF-matrix")
            addWorksheet(wb, "Variant calls")
            addWorksheet(wb, "Variant calls (SNVs)")
            addWorksheet(wb, "Variant calls (CNVs)")
            #addWorksheet(wb, "Purity")
            
            # Settings 
            # All selected values are compared to their actual value. Both are printed.
            if(is.null(input$cnvcheckbox)){ # handling checkbox input
                selvalues = c(input$evolutionvalue, input$samplevalue, input$clonecount, input$mutationcount, input$cnvcount, "-",
                              input$coveragevalue, input$purityvalue, input$detectionvalue, input$distancevalue)
                
            } else {
                cnvvalues <- input$cnvcheckbox
                cnvvalues <- paste(cnvvalues, collapse = ",")
                selvalues = c(input$evolutionvalue, input$samplevalue, input$clonecount, input$mutationcount, input$cnvcount, cnvvalues,
                              input$coveragevalue, input$purityvalue, input$detectionvalue, input$distancevalue)
            }
            if(is.null(cnvs)){
                actvalues = c("", ncol(clones)-2, nrow(clones), nrow(allmutations), 0, "", "", 
                              "", "", "")
            } else {
                actvalues = c("", ncol(clones)-2, nrow(clones), nrow(allmutations), nrow(cnvs), "", "", 
                              "", "", "")
            }
            settings <- data.frame(
                Parameter = c("Model of evolution", "#Time points", "#Clones", "#Variants (SNV+CNV)", "#CNVs", "Specification",
                            "Mean coverage", "Purity", "Detection threshold", "Minimum clonal distance"),
                Selected = selvalues,
                Actual = actvalues
            )
            names(settings)<-c("Parameter","Initially selected","Fine-tuned")
            writeData(wb, "Input settings", settings)
            #writeData(wb, "Input settings", paste("Warning: These Settings could have been overwritten by your manual changes on tables (see differences in Actual)."), startCol = 1, startRow = 15)
            #writeData(wb, "Input settings", paste("- means that parameter is not trackable."), startCol = 1, startRow = 16)
            setColWidths(wb, "Input settings", c(1,2,3), widths = "auto")
            style_head<-createStyle(border = "bottom",textDecoration = "bold",fontSize = 12,borderStyle = "medium")
            style_line<-createStyle(border="bottom")
            addStyle(wb, "Input settings",style=style_head,rows=c(1),cols=c(1:3))
            addStyle(wb, "Input settings",style=style_line,rows=c(7),cols=c(1:3))
            
            #addStyle(wb, "Settings",style = createStyle(halign = "center"), rows = 1:11, cols = 1:3, gridExpand = TRUE)
            
            # Fishplot
            fishplot <- draw_fishplot2(clones)
            print(fishplot)
            insertPlot(wb, "Visualization", width = 12.51969, height = 5,startCol = 1)
            fishplot2 <- draw_fishplot2b(clones)
            print(fishplot2)
            insertPlot(wb, "Visualization", width = 8, height = 5,startCol = 16)
            dev.off()
            
            # Phylogeny
            writeData(wb, "CCF-matrix", clones)
            addStyle(wb, "CCF-matrix",style=style_head,rows=c(1),cols=c(1:length(clones[1,])))
            
            writeData(wb, "Variant calls", allmutations)
            addStyle(wb, "Variant calls",style=style_head,rows=c(1),cols=c(1:length(allmutations[1,])))
            
            mutations_export<-mutations
            mutations_export$Overlap[is.na(mutations_export$Overlap)]<-"NA"
            writeData(wb, "Variant calls (SNVs)", mutations_export)
            addStyle(wb, "Variant calls (SNVs)",style=style_head,rows=c(1),cols=c(1:length(mutations[1,])))
            
            cnvs_export<-cnvs
            cnvs_export$Overlap<-as.character(cnvs_export$Overlap)
            cnvs_export$Scenario[is.na(cnvs_export$Scenario)]<-"NA"
            writeData(wb, "Variant calls (CNVs)", cnvs_export)
            addStyle(wb, "Variant calls (CNVs)",style=style_head,rows=c(1),cols=c(1:length(cnvs[1,])))
            #writeData(wb, "Purity", purities)
            
            # Files are saved in zip and in local directory, therefore I created directory Output to seperate from code files
            if (!file.exists("Output")){
              dir.create(file.path("Output"))
            }
            saveWorkbook(wb, file = "Output/Data.xlsx", overwrite = TRUE)
            
            # Single tables for easy further processing
            write.table(clones, file = "Output/CCF-matrix.tsv", row.names = FALSE, sep ="\t")
            write.table(cnvs, file = "Output/CNVs.tsv", row.names = FALSE, sep ="\t")
            write.table(mutations, file = "Output/SNVs.tsv", row.names = FALSE, sep ="\t")
            #write.table(purities, file = "Output/purity.tsv", row.names= FALSE, sep="\t")
            
            if (!file.exists("Output/RDS")){
                dir.create(file.path("Output/RDS"))
            }
            
            # Save current session in RDS directory for restoring possibility
            saveRDS( reactiveValuesToList(input) , file = 'Output/RDS/Input_Shiny.RDS')
            saveRDS( clones, "Output/RDS/CCF-matrix.RDS")
            saveRDS( mutations, "Output/RDS/SNVs.RDS")
            saveRDS( cnvs, "Output/RDS/CNVs.RDS")
            saveRDS( purities, "Output/RDS/Purity.RDS")
            saveRDS( showcnvs, "Output/RDS/CNVs_VariantCalls.RDS")
            
            fs = c("Data.xlsx", "Output/Data.xlsx", "Output/CCF-matrix.tsv", "Output/CNVs.tsv", "Output/SNVs.tsv","Output/RDS/Input_Shiny.RDS", 
                   "Output/RDS/CCF-matrix.RDS", "Output/RDS/SNVs.RDS",
                   "Output/RDS/CNVs.RDS", "Output/RDS/Purity.RDS", "Output/RDS/CNVs_VariantCalls.RDS")
            
            zip(zipfile = file, files = fs)
        },
        contentType = "appclication/zip"
    )
    
}