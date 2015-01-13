#######################################################
### R code for classification of Response Dynamics ####
########## Kristoffer Vitting-Seerup feb 2014 #########
#######################################################

library('plyr')

### The input is simply a data.frame where colum names are time in minuts, rownames are feature (promoter, enhancer etc. - they must be unique) and the values are the log2FC compared to t0. (0 min should be included and a pseudocount of 1 should be added)
### An example could be:
# myExampleData <- data.frame( 
#                               '0'  =c(0,0,0) ,
#                               '30' =c(1,2,0) ,
#                               '60' =c(0,2,0) ,
#                               '120'=c(0,2,0) ,
#                               '240'=c(0,2,0) ,
#                               '360'=c(0,2,0)
#                             )
# colnames(myExampleData) <- c('0','30','60','120','240','360') # note that this toy dataset have fewer timepoints than the data analyzed
# rownames(myExampleData) <- paste('feature_',1:3,sep='')


### The function that does the actual classification each entry
classifyDynamics <- function(myDf) {
    
    # extract timepoints
    timePoints <- as.numeric(colnames(myDf))
    
    # Devide timepoints (their indexes) into sections/groups that can be used assess specific parts of the timecourse
    rapidTimepointIndex     <- which( timePoints  != 0   & timePoints <= 30  )
    earlyTimepointIndex     <- which( timePoints  >  30  & timePoints <= 60  )
    midTimePointIndex       <- which( timePoints  >  60  & timePoints <= 120 )
    earlyMidtTimePointIndex <- which( timePoints  >  30  & timePoints <= 120 )
    lateTimePointIndex      <- which( timePoints  >= 100 & timePoints <= 150 )
    endTimePointIndex       <- which( timePoints  >= 240                     )
    
    ### Split the dataframe into each feature (since a feature migth be calssified into multiple calsses)
    mySplit <- split(myDf, f=rownames(myDf))
    mySplit <- lapply(mySplit, function(x) as.vector(unlist(x)) )
    
    ### Functions to analyze each of the 8 core responses identified
    analyzeRapidShortResponse       <- function(myVec) {
        # Since there are two posibilities (respons at 15 min and response at 30 min ) and the directionality depends on it I have an if/else if statements to deside this:
        earlyResp <- FALSE  # Rappid respons variable
        
        # Response possiblility 1: response at 15 min:
        if( 
            head(myVec[rapidTimepointIndex],1)  >=  1     |
            head(myVec[rapidTimepointIndex],1)  <= -1                               # Must have expression above FC 1 at 15 min
        ) {
            earlyResp <- TRUE                                                       # set rappid response variable
            directionOfDynamics <- sign( head(myVec[rapidTimepointIndex],1) )       # set direction variable
        
        # Response possiblility 2: response at 30 min: 
        } else {
            if( sign( tail(myVec[rapidTimepointIndex],1) ) == 1 ) {               
                if(
                    head(myVec[rapidTimepointIndex],1)  >= -0.25 &                  # 15 min have to be larger than -0.25 &
                    tail(myVec[rapidTimepointIndex],1)  >=  1                       # 30 min must be larger than 1
                ) {
                    earlyResp <- TRUE                                               # set rappid response variable
                    directionOfDynamics <- 1                                        # set direction variable
                }
            } else {
                if(
                    head(myVec[rapidTimepointIndex],1)  <=  0.25 &                  # 15 min have to be smaler than 0.25 &
                    tail(myVec[rapidTimepointIndex],1)  <= -1                       # 30 min must be smaler than 1
                ) {
                    earlyResp <- TRUE                                       # set rappid response variable
                    directionOfDynamics <- -1                               # set direction variable
                }
            }
        }
        # Use the result of the rappid reponse inqury above
        if(
            earlyResp       
        ) {
            # Upregulated
            if( directionOfDynamics == 1 ) {
                if( 
                    mean( myVec[ lateTimePointIndex ] ) <= max(myVec[rapidTimepointIndex]) *0.50        # Does the average FC of the late dynamics go below 1/2 of the max initial FC
                ) {
                    return('1 - Rapid Short Response - upregulated')
                }
                # Downregulated 
            } else {
                if( 
                    mean( myVec[ lateTimePointIndex ] ) >= min(myVec[rapidTimepointIndex]) *0.50            # Does the average FC of the late dynamics go below 1/2 of the max initial FC
                ) { 
                    return('1 - Rapid Short Response - downregulated')
                }               
            }
        }
    }
    analyzeRapidLongResponse        <- function(myVec) {
        ### Test the early large response (same as for RapidShortResponse above)
        # Since there are two posibilities (respons at 15 min and response at 30 min ) and the directionality depends on it I have an if/else if statements to deside this
        earlyResp <- FALSE  # Rappid respons variable
        # Possiblility 1: response at 15 min:
        if( 
            head(myVec[rapidTimepointIndex],1)  >=  1     |
            head(myVec[rapidTimepointIndex],1)  <= -1                               # Must have expression above 1 at 15 min
        ) {
            earlyResp <- TRUE                                                       # set rappid response variable
            directionOfDynamics <- sign( head(myVec[rapidTimepointIndex],1) )       # set direction variable
            # Possiblility 2: response at 30 min: 
        } else {
            if( sign( tail(myVec[rapidTimepointIndex],1) ) == 1 ) {               
                if(
                    head(myVec[rapidTimepointIndex],1)  >= -0.25 &                  # 15 min have to be larger than -0.25 &
                    tail(myVec[rapidTimepointIndex],1)  >=  1                       # 30 min must be larger than 1
                ) {
                    earlyResp <- TRUE                                               # set rappid response variable
                    directionOfDynamics <- 1                                        # set direction variable
                }
            } else {
                if(
                    head(myVec[rapidTimepointIndex],1)  <=  0.25 &                  # 15 min have to be smaler than 0.25 &
                    tail(myVec[rapidTimepointIndex],1)  <= -1                       # 30 min must be smaler than 1
                ) {
                    earlyResp <- TRUE                                               # set rappid response variable
                    directionOfDynamics <- -1                                       # set direction variable
                }
            }
        }
        
        ### Use the result of the rappid reponse inqury above
        if(
            earlyResp       
        ) {
            # Upregulated
            if( directionOfDynamics == 1 ) {
                ### Rapid Respons 2) the long response
                if (
                    all( myVec[ earlyMidtTimePointIndex  ] >= 0.25 )                # Does the midrage FC stay large ?
                ) {
                    return( '2 - Rapid Long Response - upregulated' )
                }
                # Downregulated 
            } else {
                ### Rapid Respons 2) the long response
                if (
                    all( myVec[ earlyMidtTimePointIndex  ] <= -0.25 )               # Does the midrage FC stay large ?
                ) {
                    return( '2 - Rapid Long Response - downregulated' )
                }
            }
        }
    }
    analyzeEarlyStandardResponse    <- function(myVec) {
        if (
            # Is there a peak in the early response region
            ! all( c(1,-1) %in% sign(myVec[earlyTimepointIndex]) )                                              # are all early timepoints regulated in the same way
        ) {
            
            # Upregulated
            if( 1 %in% sign(myVec[ earlyTimepointIndex ]) ) {
                if( 
                    # Small start
                    mean( myVec[rapidTimepointIndex    ])   <= 1                  &                             # Do they start small
                    # large response
                    all(  myVec[earlyTimepointIndex    ]    >= 0.25 )             &                             # Is all of them is above 0.25 ?
                    # Down agian
                    all(  myVec[midTimePointIndex]          <= max(myVec[earlyTimepointIndex]) *0.50 )          # Does the expression go down again?
                ) {
                    return('3 - Early Standard Response - upregulated')
                }
            # Downregulated
            } else {
                if( 
                    # Small start
                    mean( myVec[rapidTimepointIndex    ])   >= -1                  &                            # Do they start small
                    # large response
                    all(  myVec[earlyTimepointIndex    ]    <= -0.25           )   &                            # Is all of them is above 0.25 ?
                    # Down agian
                    all(  myVec[midTimePointIndex]          >= min(myVec[earlyTimepointIndex]) *0.50 )          # Does the expression go down again?
                ) {
                    return('3 - Early Standard Response - downregulated')
                }
            }
        }
    }
    analyzeLateStandardResponse     <- function(myVec) {
        if (
            ! all( c(1,-1) %in% sign(myVec[midTimePointIndex]) )                                            # are all midpoints regulated in the same way values?
        ) {
            # Upregulated
            if( 1 %in% sign(myVec[ midTimePointIndex ]) ) {
                if( 
                    # Small start
                    mean( myVec[ which( timePoints  >=   0 & timePoints < 60  ) ])  <=  1           &       # Do they start small
                    # Large response
                    all(  myVec[midTimePointIndex  ]    >= 0.25            )                        &       # Is all of the midpoints is above 0.25 ?
                    # Down again
                    all(  myVec[endTimePointIndex  ]    <=  max(myVec[midTimePointIndex]) *0.50  )          # Does the expression go down again?
                ) {
                    return( '4 - Late Standard Response - upregulated' )
                }
            # Downregulated
            } else {
                if( 
                    # Small start
                    mean( myVec[ which( timePoints  >=   0 & timePoints < 60  ) ])  >=  -1          &       # Do they start small
                    # Large response
                    all(  myVec[midTimePointIndex  ]    <= -0.25           )                        &       # Is all of the midpoints is larger than 0.25 ?
                    # Down again
                    all(  myVec[endTimePointIndex  ]    >=  min(myVec[midTimePointIndex]) *0.50  )          # Does the expression go down again?
                ) {
                    return( '4 - Late Standard Response - downregulated')
                }
            }
        }
    }
    analyzeLongResponse             <- function(myVec) {
        if( 
            ! all( c(1,-1) %in% sign(myVec[ earlyMidtTimePointIndex ]) )            # are all early and midt timepoints regulated in the same way
        ) {
            if( 1 %in% sign(myVec[ earlyMidtTimePointIndex]) ) {
                if( 
                    all(  myVec[ earlyMidtTimePointIndex ]    >=  0.25  )           # Are all of them is above 0.25 ?
                ) {
                    return( '5 - Long Response - upregulated' ) 
                }
            } else {
                if(
                    all(  myVec[ earlyMidtTimePointIndex ]    <= -0.25  )           # Are all of them is above 0.25 ?
                ) {
                    return( '5 - Long Response - downregulated' ) 
                }
                
            }
        }
    }
    analyzeLateResponse             <- function(myVec) {
        # Upregulated
        if( mean( myVec[endTimePointIndex] ) > 0 ) {
            if( 
                # make sure there are large change in late timepoints
                mean( myVec[endTimePointIndex   ] )  >= max(myVec) *0.50        &       # Does the end have large changes ?
                all(  myVec[endTimePointIndex   ]    >= 0.25              )             # Does the end have any small FCs?
            ) {
                return( '6 - Late Response - upregulated' )
            }
            # Downregulated
        } else {
            if ( 
                # make sure there are large change in late timepoints
                mean( myVec[endTimePointIndex   ] )  <= min(myVec) *0.50         &       # Does the end have large changes ?
                all(  myVec[endTimePointIndex   ]    <= -0.25              )           # Does the end have any small FCs?
            ) {
                return( '6 - Late Response - downregulated' )
            }
        }
    }
    analyzeLateFlatResponse         <- function(myVec) {
        if( 
            mean( myVec[which(timePoints >= 120) ] ) <=   0.75    &       # Is the larger values small 
            mean( myVec[which(timePoints >= 120) ] ) >=  -0.75            # Is the larger values small 
        ) {
            return( '7 - Late Flat Response - NA')
        }
    }
    analyzeEarlyFlatResponse        <- function(myVec) {
        if( 
            mean( myVec[ which(timePoints <= 60) ] ) <=   0.75    &       # Is the larger values small 
            mean( myVec[ which(timePoints <= 60) ] ) >=  -0.75            # Is the larger values small 
        ) {
            return( '8 - Early Flat Response - NA' )
        }
    }
    
    ### Wraper function to use the functions above to analyze all of the 8 core responses
    heracicalClassification <- function(vec) {
        ### Use the functions to analyze the different responses
        resultOfRapidShortResponse      <- analyzeRapidShortResponse(vec)
        resultOfRapidLongResponse       <- analyzeRapidLongResponse(vec)
        resultOfEarlyStandardResponse   <- analyzeEarlyStandardResponse(vec)
        resultOfLateStandardResponse    <- analyzeLateStandardResponse(vec)
        resultOfLongResponse            <- analyzeLongResponse(vec)
        resultOfLateResponse            <- analyzeLateResponse(vec)
        resultOfLateFlatResponse        <- analyzeLateFlatResponse(vec)
        resultOfEarlyFlatResponse       <- analyzeEarlyFlatResponse(vec)
        
        # Vector to store clasification(s)
        classification <- NULL
        
        ### Check for each of the allowed responses
        # Raid short responses
        if( !is.null(resultOfRapidShortResponse) ) {
            # + late flat
            if( !is.null(resultOfLateFlatResponse) ) {
                classification <- c(classification, paste(resultOfRapidShortResponse, resultOfLateFlatResponse, sep=',') )
            }
            # + long
            if( !is.null(resultOfLongResponse) ) {
                classification <- c(classification, paste(resultOfRapidShortResponse, resultOfLongResponse,     sep=',') )
            }
            # + late
            if( !is.null(resultOfLateResponse) ) {
                classification <- c(classification, paste(resultOfRapidShortResponse, resultOfLateResponse,     sep=',') )
            }
        }
        # rapid long response
        if( !is.null(resultOfRapidLongResponse) ) {
            classification <- c(classification, resultOfRapidLongResponse)
        }
        # Early Standard response
        if( !is.null(resultOfEarlyStandardResponse) ) {
            # + late flat
            if( !is.null(resultOfLateFlatResponse) ) {
                classification <- c(classification, paste(resultOfEarlyStandardResponse, resultOfLateFlatResponse, sep=',') )
            }
            # + late
            if( !is.null(resultOfLateResponse) ) {
                classification <- c(classification, paste(resultOfEarlyStandardResponse, resultOfLateResponse,     sep=',') )
            }
        }
        # Late standard Response
        if( !is.null(resultOfLateStandardResponse) ) {
            classification <- c(classification, resultOfLateStandardResponse)
        }
        # Long response
        if( !is.null(resultOfLongResponse) ) {
            classification <- c(classification, resultOfLongResponse)
        }
        # Late response
        if( !is.null(resultOfLateResponse) & !is.null(resultOfEarlyFlatResponse) ) { # make sure when to classifiy the late response alone its a flat start
            classification <- c(classification, resultOfLateResponse)
        }
        # If the response dynamic have not been classified by now return unclassified
        if(length(classification) == 0) {
            classification <- '9 - Unclassified - NA'
        }
        
        myResult <- data.frame(t(vec), classification=classification, stringsAsFactors=F)
        return(myResult)
    } # end of Heracical Classification function
    
    ### Use the Heracical classification function on the data.frame
    myResults <- ldply(mySplit, heracicalClassification)
    colnames(myResults)[-ncol(myResults)] <- c('feature', colnames(myDf))
    
    # Return result
    return(myResults)
}

### To use the classification function on the example data simply run:
# myClassificationResult <- classifyDynamics(myExampleData)
