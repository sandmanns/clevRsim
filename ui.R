ui <- fluidPage(
  theme = shinytheme("flatly"),
  titlePanel(div(img(height = 100, width = 100, src = "Logo.png",class="pull-left"),
                 "clevRsim")
             
            # span(div(
            #   column(1, img(src = "Logo.png", width = '80px')),
            #   column(10,"clevR-sim",style = "font-size:40px;"),
            #   column(1,img(src = "IMI.png", width = "100px"))
             #))
  ),
  tabsetPanel(
                 #Settings_UI ####
                 tabPanel("Input", 
                          fluidRow(
                              column(7,
                                     plotOutput("fishPlot")
                                     ),
                              column(4,
                                     plotOutput("fishPlot2")
                                     ),
                              
                              column(2,
                                     h3("Basic settings"),
                                     selectInput('evolutionvalue', 'Model of evolution', 
                                                 c("linear","branched dependent","branched independent"),selected = "linear", width = "90%"),
                                     numericInput('samplevalue', '#Time points', 
                                                 min=1, max=50, value=3, width = "90%"),
                                     numericInput('clonecount', '#Clones', min=1, max=50, value = 5, width = "90%"),
                                     numericInput('mutationcount', '#Variants (SNVs+CNVs)', 
                                                 min=1, max=200, value =20, width = "90%"),
                                     sliderInput("cnvcount", "#CNVs (included in #Variants)", min=0, max=30, value=0, 
                                                 step=1, ticks = FALSE, width = "90%"),
                                     checkboxGroupInput("cnvcheckbox", "(opt.) Specification", choices = c("Deletion", "Duplication", "LOH"), 
                                                        selected = NULL,inline = T, width = "90%")#,
                                     #h6("When no specific type of CNV is defined, random CNVs are simulated.", width = "90%")
                                     #checkboxInput('fp', 'Include false positives', width = "100%")
                              ),
                              column(2, 
                                     br(),
                                     h4("Advanced settings"),
                                     numericInput('coveragevalue', 'Mean coverage', min = 1, max = 1000,
                                                  value = 300, width = "90%"),
                                     sliderInput("purityvalue", "Purity [%]", min = 1, max = 100, value = 100,
                                                 ticks = FALSE, width = "90%"),
                                     numericInput('detectionvalue', 'Detection threshold [%CCF]', min = 0, max = 50,
                                                  value = 2, width = "90%"),
                                     numericInput('distancevalue', 'Minimum clonal distance [%CCF]', min = 0, max = 50,
                                                  value = 4, width = "90%")
                              ),
                              column(2, 
                                     br(),br(),br(),br(),
                                     actionButton("simulatebtn","Run simulation",width = "90%",class = "btn-primary"),
                                     br(),br(),br(),
                                     fileInput("upload", label=("Load session (6 RDS-files)"), accept = ".RDS", multiple = TRUE, width = "90%"),
                                     # actionButton("loadbtn", "Load Session"),
                              ),
                              column(5,
                                     br(),
                                     uiOutput("uiheader"),
                                     br(),
                                     rHandsontableOutput("clonetable"), 
                                     br(),
                                     actionButton("updatebtn1","Update simulation",class = "btn-primary"),
                                     br(), br(),
                                     uiOutput("updatetext1"),
                                     br(), br(),
                                     fluidRow(
                                       column(3,
                                              uiOutput("uiaddsample"),
                                              br(),
                                              uiOutput("uideletesample"),
                                              br(),br(),br()),
                                       column(2,
                                              uiOutput("uiaddsamplevalue"),
                                              br(),
                                              uiOutput("uideletesamplevalue")),
                                       column(1),
                                       column(3,
                                              uiOutput("uiaddclone"),
                                              br(),
                                              uiOutput("uideleteclone")),
                                       column(2,
                                              uiOutput("uiaddclonevalue"),
                                              br(),
                                              uiOutput("uideleteclonevalue"))
                                         #column(3, br(), 
                                        #        uiOutput("uiaddsample"),
                                        #        br(),
                                        #        uiOutput("uiaddclone"),
                                        #        br() 
                                        # ),
                                         #column(3, #offset = 1,
                                        #        uiOutput("uideletesamplevalue"),
                                        #        uiOutput("uideletesample")
                                        #),
                                         #column(3,
                                        #        uiOutput("uideleteclonevalue"),
                                        #        uiOutput("uideleteclone")
                                        # )
                                     )
                              )
                          )
                 ),
                 # Data_UI ####
                 tabPanel("Output",  
                          #h3("Table of all Variants"),
                          dataTableOutput("mutationtable"), 
                          hr(),
                          #fluidRow(
                          #    column(11),
                          #    column(1,
                          #           uiOutput("uiexport2")
                          #           )
                          #),
                          fluidRow(
                            column(12,
                                     fluidRow(
                                       column(5,
                                              uiOutput("uicnvtitle"), 
                                              uiOutput("uicnvtitle2"), 
                                              uiOutput("uicnvtitle2a"), 
                                              uiOutput("uicnvtitle2b"),
                                              uiOutput("uicnvtitle2c"),
                                              uiOutput("uicnvtitle2d"),
                                              uiOutput("uicnvtitle2e")
                                              ),
                                       column(7,
                                              br(),
                                              fluidRow(
                                                column(10,
                                                       uiOutput("uicnvtitle_scenarios1")),
                                                column(2,
                                                       uiOutput("uiexport2")
                                                       )
                                              ),br(),
                                              fluidRow(
                                                column(3,
                                                       uiOutput("uicnvtitle_scenarios3a"),
                                                       uiOutput("uicnvtitle_scenarios3b"),
                                                       uiOutput("uicnvtitle_scenarios3c"),
                                                       uiOutput("uicnvtitle_scenarios3d")
                                                       ),
                                                column(3,
                                                       uiOutput("uicnvtitle_scenarios4a"),
                                                       uiOutput("uicnvtitle_scenarios4b"),
                                                       uiOutput("uicnvtitle_scenarios4c"),
                                                       uiOutput("uicnvtitle_scenarios4d")
                                                       ),
                                                column(3,
                                                       uiOutput("uicnvtitle_scenarios5a"),
                                                       uiOutput("uicnvtitle_scenarios5b"),
                                                       uiOutput("uicnvtitle_scenarios5c"),
                                                       uiOutput("uicnvtitle_scenarios5d")
                                                ),
                                                column(3,
                                                       uiOutput("uicnvtitle_scenarios6a"),
                                                       uiOutput("uicnvtitle_scenarios6b"),
                                                       uiOutput("uicnvtitle_scenarios6c"),
                                                       uiOutput("uicnvtitle_scenarios6d")
                                                )
                                              )
                                       )
                                     ),
                                     
                                     br(), 
                                     rHandsontableOutput("cnvtable"),
                                     br(),
                                     fluidRow(
                                         column(1,
                                                actionButton("addcnvbtn", "Add CNV"),
                                                br(),br(),
                                                actionButton("deletecnvbtn", "Delete CNV")
                                         ),
                                         column(1,
                                                uiOutput("uiaddcnvvalue"),
                                                #selectInput("addcnvvalue", NULL, selected = 1, choices = 1, width = "80%"),
                                                br(),
                                                uiOutput("uideletecnvvalue")
                                                #selectInput("deletecnvvalue", NULL, selected = 1, choices = 1:10, width = "80%")
                                         ),
                                         column(1),
                                         column(4, offset = 1,
                                                actionButton("updatebtn2", "Update simulation",class = "btn-primary"),
                                                br(), br(),
                                                uiOutput("updatetext2")
                                         )
                                     ),
                                     br()
                              )#,
                              #column(1,  br(), uiOutput("uiexport2"), br(), br()
                              #)
                          )#,
                          #fluidRow(
                          #    column(4, br(),
                          #           uiOutput("uipuritytitle"), br(), 
                          #           rHandsontableOutput("puritytable"),
                          #           br()
                          #    )
                          #)
                 ),
                 tabPanel("Manual",
                          fluidRow(
                            column(7,
                                   h1("clevRsim: Simulating clonal evolution in R"),
                                   br(),
                                   h5("clevRsim provides an R shiny interface for simulation of clonal evolution. The main focus of the tool lies on the
                                      realistic simulation of variant calls (point mutations as well as copy number variants) characterizing each clone. The output
                                      generated by clevRsim can be directly taken as input to train and test algorithms for variant clustering and clonal evolution
                                      tree reconstruction."),
                                   br(),br(),
                                   h2(HTML("<u>Input</u>")),br(),
                                   h4(HTML("<b>Basic settings</b>")),
                                   h5(HTML(rep("&nbsp",4),"<b>Model of evolution:</b>")),
                                   h5(HTML(rep("&nbsp",8)),"Linear: Linear development of clonal evolution."),
                                   h5(HTML(rep("&nbsp",8)),img(height = 150, width = 600,src="Linear.png")),
                                   h5(HTML(rep("&nbsp",8)),"Branched dependent: Branching development of clonal evolution. All clones
                                      develop from one joined parent clone (=developing from normal cells). A minimum of #Clones=3 ",HTML("<br>"),HTML(rep("&nbsp",8)),
                                      "is required. The number of branches is randomly simulated, but may be customized subsequently."),
                                   h5(HTML(rep("&nbsp",8)),img(height = 150, width = 600,src="Dependent.png")),
                                   h5(HTML(rep("&nbsp",8)),"Branched independent: Branching development of clonal evolution. At least two parent clones (=developing
                                      from normal cells) exist. A minimum of #Clones=2 is required. ",HTML("<br>"),HTML(rep("&nbsp",8)),
                                      "The number of branches is randomly simulated, but may be customized subsequently."),
                                   h5(HTML(rep("&nbsp",8)),img(height = 150, width = 600,src="Independent.png")),
                                   
                                   h5(HTML(rep("&nbsp",4),"<b>#Time points:</b>")," Number of time points to be simulated."),
                                   h5(HTML(rep("&nbsp",4),"<b>#Clones:</b>")," Number of clones to be simulated."),
                                   h5(HTML(rep("&nbsp",4),"<b>#Variants (SNV+CNV):</b>")," Number of variants to be simulated (summing up SNVs and CNVs)"),
                                   h5(HTML(rep("&nbsp",8)),"#CNVs: Number of CNVs to be simulated (minimum=0, maximum=#Variants)"),
                                   h5(HTML(rep("&nbsp",8)),"(opt.) Specification: Specification of the CNVs to be simulated. Possible values are: Deletion, 
                                      Duplication, LOH. More than 1 specification can be selcted. If no specification is",HTML("<br>"),HTML(rep("&nbsp",8)),
                                      "selected, types of CNVs are simulated randomly."),
                                   br(),
                                   
                                   h4(HTML("<b>Advanced settings</b>")),
                                   h5(HTML(rep("&nbsp",4),"<b>Mean coverage:</b>")," Mean coverage to be used for simulating SNVs (#Reference reads, #Variant
                                   reads). Coverage is simulated assuming a log-normal distriution."),
                                   h5(HTML(rep("&nbsp",4),"<b>Purity:</b>")," Purity of all samples to be simulated."),
                                   h5(HTML(rep("&nbsp",4),"<b>Detection threshold:</b>")," Thresholds for the detection of variants. It is assumed that only 
                                      variants with CCF>Detection threshold can be detected. Thus, no variants with CCF below",HTML("<br>"),HTML(rep("&nbsp",4)),
                                      "this threshold are simualted."),
                                   h5(HTML(rep("&nbsp",4),"<b>Minimum clonal distance:</b>")," Minimum difference in CCF between two clones. 
                                      If variants characterizing two clones differ by less than Minimum clonal distance, it is assumed that",HTML("<br>"),HTML(rep("&nbsp",4)),
                                      "they can no longer beclassified as belonging to two different clones. The threshold is applied to 
                                      all clones at every time point. To allow presence of two clones with the same ",HTML("<br>"),HTML(rep("&nbsp",4)),
                                      "CCFs, Minimum clonal distance may be set to zero."),
                                   br(),
                                   
                                   h4(HTML("<b>Fine-tuning: Clones</b>")),
                                   h5("The following options are available after running the simulation for the first time."),
                                   h5(HTML(rep("&nbsp",4),"<b>Modify CCF-matrix:</b>")),
                                   h5(HTML(rep("&nbsp",8)),"Change CCFs: Optionally change the randomly simulated CCFs at each time point for each clone. When
                                      clicking Update simulation a validity check is performed. If manually",HTML("<br>"),HTML(rep("&nbsp",8)),
                                      "defined CCFs do not allow for simulating clonal evolution under the given scenario, an error-message reportes the 
                                      problematic clone(s) and time point(s). The initially",HTML("<br>"),HTML(rep("&nbsp",8)),
                                      "simulated clonal evolution is not updated."),
                                   h5(HTML(rep("&nbsp",8)),"Change Parent:"," Optionally change the simulated parent for each clone. When
                                      clicking Update simulation a validity check is performed. If manually defined parents do not",HTML("<br>"),HTML(rep("&nbsp",8)),
                                      "allow for simulating clonal evolution under the given scenario, an error-message reportes the 
                                      problematic clone(s) and time point(s). The initially simulated clonal evolution is",HTML("<br>"),HTML(rep("&nbsp",8)),
                                      "not updated."),
                                   h5(HTML(rep("&nbsp",8)),"Change #Variants:"," Optionally change the simulated number of variants for each clone. When
                                      clicking Update simulation a validity check is performed. If the sum of the",HTML("<br>"),HTML(rep("&nbsp",8)),
                                      "manually defined number of variants for each clone does not match #Variants, an error-message is reported. 
                                      The initially simulated clonal evolution is not updated."),
                                   h5(HTML(rep("&nbsp",4),"<b>Add time point:</b>"),"Add a time point to the initially simulated clonal evolution. 
                                   The time point is always inserted after the last time point (at that time), with
                                      CCF=0 for all clones."),
                                   h5(HTML(rep("&nbsp",4),"<b>Add clone:</b>"),"Add a clone to the initially simulated clonal evolution. The clone is always inserted 
                                   after the last clone (at that time), deriving from normal cells (Parent=0),",HTML("<br>"),HTML(rep("&nbsp",4)),
                                      "characterized by no variants (#Variants=0), CCF=1 at time point 1 (t1) and CCF=0 at all other time points. 
                                      These default settings can be changed by modifying the CCF-matrix."),
                                   h5(HTML(rep("&nbsp",4),"<b>Delete time point:</b>"),"Delete a time point of the simulated clonal evolution. 
                                   Any time point can be chosen."),
                                   h5(HTML(rep("&nbsp",4),"<b>Delete clone:</b>"),"Delete a clone of the simulated clonal evolution. Any clone can be
                                      chosen."),
                                   br(),
                                   h4(HTML("<b>Visualization</b>")),
                                   h5("Two plots visualizing simulated clonal evolution are available after running the simulation for the first time. 
                                      The plots are updated every time a new simulation is run, or the current simulation is updated according to parameters
                                      available in Fine-tuing: Clones. "),
                                   h5("clevR-vis is used for simulation. A dolphin plot (phylogeny + CCFs + time course) is visualized along with an extended
                                      shark plot (CCFs + time course)."),
                                   br(),
                                   
                                   h4(HTML("<b>Import</b>")),
                                   h5(HTML(rep("&nbsp",4),"<b>Load session (6 RDS-files):</b>"),"A simulation can be retrieved by uploading the 6 RDS-files generated 
                                      using Output > Export data."),
                                   
                                   br(),br(),
                                   h2(HTML("<u>Output</u>")),br(),
                                   h4(HTML("<b>Variant calls</b>")),
                                   h5(HTML(rep("&nbsp",4),"<b>Basic characteristics:</b>")),
                                   h5(HTML(rep("&nbsp",8)),"ID: Unique variant identifier (for SNVs and CNVs)."),
                                   h5(HTML(rep("&nbsp",8)),"Chr: Simulated chromosomal location (for SNVs and CNVs)."),
                                   h5(HTML(rep("&nbsp",8)),"Start: Start position of the simulated variant (for SNVs and CNVs)."),
                                   h5(HTML(rep("&nbsp",8)),"End: End postiion of the simulated variant (for SNVs and CNVs)."),
                                   h5(HTML(rep("&nbsp",4),"<b>Advanced characteristics:</b>")),
                                   h5(HTML(rep("&nbsp",8)),"Genotype: Genotype of the simulated variant (only for SNVs). In case of no CNV overlapping an SNV, the
                                   genotype is reported as 'AB (AA->AB)'. For scenario 'CNV first' the",HTML("<br>"),HTML(rep("&nbsp",8)),
                                      "genotype is reported as 'B (AA->A->B)' for deletions as 'AAB (AA->AAA->AAB)' for duplications and 'AB (AA->AA(LOH)->AB)' 
                                      for LOH. For scenario 'SNV first (affected)'",HTML("<br>"),HTML(rep("&nbsp",8)),"
                                     the genotype is reported as 'AB, A (AA->AB->A)' for deletions, as 'AB, ABB (AA->AB->ABB)' for duplications and as 'AB, BB (AA->AB->BB)' 
                                      for LOH. For scenario 'SNV first",HTML("<br>"),HTML(rep("&nbsp",8)),
                                      "(un-affected)' the genotype is reported as 'AB, B (AA->AB->B)' for deletions, as
                                      'AB, AAB (AA->AB->AAB)' for duplications and as 'AB, AA (AA->AB->AA)' for LOH. For",HTML("<br>"),HTML(rep("&nbsp",8)),
                                      "scenario 'Parallel' the genotype is
                                      reported as 'AB (AA->AB; AA->A)' for deletions, as 'AB (AA->AB; AA->AAA)' for duplications and as 'AB (AA->AB; AA->AA(LOH)'
                                      for LOH."),
                                   h5(HTML(rep("&nbsp",8)),"Overlap: ID of the overlapping CNV (only for SNVs)."),
                                   h5(HTML(rep("&nbsp",8)),"Clone: Number of the clone to which a variant belongs (for SNVs and CNVs). If CNVs are defined to overlap
                                       SNVs, the 'Clone' is automatically updated to fit the selected",HTML("<br>"),HTML(rep("&nbsp",8)),
                                      "scenario."),
                                   h5(HTML(rep("&nbsp",4),"<b>For each time point:</b>")),
                                   h5(HTML(rep("&nbsp",8)),"Ref: Number of simulated reference reads, detected at the position of an SNV (only for SNVs). The influence of
                                      overlapping CNV, affecting the number of reference reads, is",HTML("<br>"),HTML(rep("&nbsp",8)),
                                      "automatically considered."),
                                   h5(HTML(rep("&nbsp",8)),"Var: Number of simulated variant reads characterizing an SNV (only for SNVs). The influence of
                                      overlapping CNV, affecting the number of variant reads, is automatically",HTML("<br>"),HTML(rep("&nbsp",8)),
                                      "considered."),
                                   h5(HTML(rep("&nbsp",8)),"VAF: Observed variant allele frequency for an SNV, calculated on the basis of 'Ref' and 'Var' (only for SNVs)."),
                                   h5(HTML(rep("&nbsp",8)),"CCF: Simulated cancer cell fraction of the clone that harbours the variant (for SNVs and CNVs). 
                                      Discrepancies between 2*VAF and CCF are intentionally. They represent",HTML("<br>"),HTML(rep("&nbsp",8)),
                                      "natural variation in the number of simulated reference 
                                      and variant reads. Furthermore, overlapping CNVs can affect the common relation between VAF and CCF."),
                                   br(),
                                   
                                   h4(HTML("<b>Fine-tuning: CNVs</b>")),
                                   h5(HTML(rep("&nbsp",4),"<b>Change characteristics:</b>")),
                                   h5(HTML(rep("&nbsp",8)),"Chr: Simulated chromosomal location (only for CNVs). Different from SNVs, this value is fixed when 
                                      updating the simulation."),
                                   h5(HTML(rep("&nbsp",8)),"Start: Start position of the simulated variant (only for CNVs). Different from SNVs, this value is fixed when 
                                      updating the simulation."),
                                   h5(HTML(rep("&nbsp",8)),"End: End position of the simulated variant (only for CNVs). Different from SNVs, this value is fixed when 
                                      updating the simulation."),
                                   h5(HTML(rep("&nbsp",8)),"Type: Type of the CNV; can be one of Deletion, Duplication or LOH."),
                                   
                                   h5(HTML(rep("&nbsp",4),"<b>Define overlap:</b>")),
                                   h5(HTML(rep("&nbsp",8)),"Overlap: Checkbox defining whether a CNV shall be simulated to overlap at least one SNV."),
                                   h5(HTML(rep("&nbsp",8)),"SNVs: Number of SNVs to overlap a CNV. The coordinates of the overlapping SNPs are automatically
                                      adjusted according to the CNV's location."),
                                   h5(HTML(rep("&nbsp",8)),"Scenario: Szenario for the overlap; can be one of 'CNV first', 'SNV first (affected)', 'SNV first (un-affected)',
                                      'Parallel'."),
                                   h5(HTML(rep("&nbsp",8)),"in_Clone: Number of the clone to which a CNV belongs. Can be changed optionally. If a change in 'in_Clone' 
                                      leads to a clone characterized by no variants at all, 'Clone' of the",
                                      HTML("<br>"),HTML(rep("&nbsp",8)),"SNVs are automatically updated. "),
                                   
                                   h5(HTML(rep("&nbsp",4),"<b>Add CNV:</b>"),"Add an empty CNV to the initially simulated clonal evolution. The CNV is always inserted 
                                   after the last CNV (at that time). All parameters (Chr, Start, End, Type,",HTML("<br>"),HTML(rep("&nbsp",4)),
                                      "Overlap, SNVs, Scenario, in_clone, CCFs for each time
                                      point) have to be defined manually."),
                                   h5(HTML(rep("&nbsp",4),"<b>Delete CNV:</b>"),"Delete a CNV of the simulated clonal evolution. Any CNV can be
                                      chosen. The total number of clones and variants is automatically adjusted to match initially defined",
                                      HTML("<br>"),HTML(rep("&nbsp",4)),"values of '#Clones' and '#Variants'."),
                                   br(),
                                   
                                   h4(HTML("<b>Export</b>")),
                                   h5(HTML(rep("&nbsp",4),"<b>Export data:</b> Simulated data can be exported. Automatically, an 'Output'-folder is generated.
                                           Files 'CCF-matrix.tsv', 'CNVs.tsv', 'SNVs.tsv' and 'Data.xlsx' are saved and can be used ",HTML("<br>"),HTML(rep("&nbsp",4)),
                                           "as direct input 
                                           for algorithms performing variant clustering and/or clonal evolution tree reconstruction. Additionally, a subfolder
                                           'RDS' is generated. Files 'CCF-matrix.RDS',",HTML("<br>"),HTML(rep("&nbsp",4)),
                                           "'CNVs.RDS', 'CNVs_VariantCalls.RDS', 'Input_Shiny.RDS', 'Purity.RDS' 
                                           and 'SNVs.RDS' are saved.")),
                                   h5(HTML(rep("&nbsp",8)),"CCF-matrix.tsv: tsv-file containing information as present in Input > Fine-tuning: Clones > CCF-matrix 
                                      (columns: Parent, #Variants, CCFs for",HTML("<br>"),HTML(rep("&nbsp",8))," each simulated time point)."),
                                   h5(HTML(rep("&nbsp",8)),"CNVs.tsv: tsv-file containing information as present in Output > Fine-tuning: CNVs 
                                      (columns: All parameters (Chr, Start, End, Type, Overlap, SNVs, Scenario, in_clone, CCFs 
                                      for",HTML("<br>"),HTML(rep("&nbsp",8)),"each time point)."),
                                   h5(HTML(rep("&nbsp",8)),"SNVs.tsv: tsv-file containing information as present in Output > Variant calls, excluding information on 
                                   simulated CNVs (columns: ID, Chr, Start, End, Genotype, Overlap,",HTML("<br>"),HTML(rep("&nbsp",8)),"Clone,
                                   Ref, Var, VAF, CCF for each time point)."),
                                   h5(HTML(rep("&nbsp",8)),"Data.xlsx: xlsx-file containing information detailed information on the simulation.  
                                      Sheet 'Input settings' contains information on the initiall defined basic and advanced",HTML("<br>"),HTML(rep("&nbsp",8)),
                                      "settings, as well as on the
                                      fine-tuned settings. Sheet 'Visualization' contains the dolphin plot and the extended shark plot visualizing the 
                                      simulated clonal evolution. Sheet",HTML("<br>"),HTML(rep("&nbsp",8)),"'CCF-matrix' contains the CCF-matrix as present in Input > Fine-tuing: Clones > CCF-matrix. 
                                      Sheet 'Variant calls' contains all variant calls as present in Output > Variant calls,",HTML("<br>"),HTML(rep("&nbsp",8)),
                                      "including information on simulated SNVs and CNVs. 
                                      Sheet 'Variant calls (SNVs)' contains all variant calls as present in Output > Variant calls, excluding information on",HTML("<br>"),HTML(rep("&nbsp",8)),
                                      "simulated CNVs. Sheet 'Variant calls (CNVs)' contains all CNVs as present in Output > Fine-tuning: CNVs."),
                                   br(),
                                   h5(HTML(rep("&nbsp",8)),"CCF-matrix.RDS: RDS-file containing information as present in Input > Fine-tuning: Clones > CCF-matrix 
                                      (R data.frame object with columns: Parent, #Variants, CCFs for each",HTML("<br>"),HTML(rep("&nbsp",8)),"simulated time point)."),
                                   h5(HTML(rep("&nbsp",8)),"CNVs.RDS: RDS-file containing information as present in Output > Fine-tuning: CNVs 
                                      (R data.frame object with columns: Chr, Start, End, Type, Overlap, SNVs, Scenario,",HTML("<br>"),HTML(rep("&nbsp",8)),
                                      "in_clone, CCFs for each time point)."),
                                   h5(HTML(rep("&nbsp",8)),"CNVs_VariantCalls.RDS: RDS-file containing information as present in Output > Variant Calls, excluding 
                                   information on simulated SNVs (R data.frame object with columns:",HTML("<br>"),HTML(rep("&nbsp",8)),
                                      "ID, Chr, Start, End, Genotype, Overlap, Clone, Ref, Var, VAF, CCF for each time point; of note Genotype, Overlap, Ref Var and VAF are left empty for CNVs)."),
                                   h5(HTML(rep("&nbsp",8)),"Input_Shiny.RDS: RDS-file containing reactive R shiny input value."),
                                   h5(HTML(rep("&nbsp",8)),"Purity.RDS: RDS-file containing information as present in Input > Advanced settings > Purity."),
                                   h5(HTML(rep("&nbsp",8)),"SNVs.RDS: RDS-file containing information as present in Output > Variant calls, excluding information on 
                                   simulated CNVs (R data.frame object with columns: ID, Chr, Start,",HTML("<br>"),HTML(rep("&nbsp",8)),"End, Genotype, Overlap, Clone,
                                   Ref, Var, VAF, CCF for each time point)."),
                                   
                                   br()
                            )
                          )
                 ),
                 tabPanel("Contact",
                          fluidRow(
                            column(3,
                                   h3("Contact"),
                                   br(),
                                   h5("PD Dr. Sarah Sandmann"),
                                   h5("sarah.sandmann@uni-muenster.de"),
                                   br(),
                                   h5("University of Münster"),
                                   br(),
                            )
                          )
                          )
                 ),
                 #for showing simulation message in better position
                 tags$head(
                     tags$style(
                         HTML(".shiny-notification {
                       position:fixed;
                       top: calc(40%);
                       left: calc(30%);
                       width: 30em;
                       }
                       "
                         )
                     )
                 ),
                 #for show and hide()
                 useShinyjs()
)








