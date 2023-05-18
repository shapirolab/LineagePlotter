library("ape")
library("dplyr") 
library("stats")
library("ggplot2")
library("ggtree")
library("ggforce")
library("yaml")
library(ggnewscale)
library(stringr)
library("ape")
library("dplyr") 
library("ggplot2")
library("ggtree")
library("ggforce")
library("yaml")
library(stringr)
library("stats")
library("stats")
library(ggnewscale)


func.print.lineage.tree <- function(conf_yaml_path,
                                    width=170,height= 60,laderize_flag=FALSE,simulate.p.value=TRUE,
                                    man_adjust_elipse=0.005, man_multiply_elipse= 1.3,
                                    man_adj_second_legend= 0,man_space_second_legend = -0.02,
                                   add_date_to_text_flag=TRUE,
                                    man_adjust_image_of_second_legend=0, man_adj_heat_loc=0,man_boot_x_offset=0,
                                    man_adj_heat_loc2=0, man_adj_heat_loc3=0,
                                    id_tip_trim_flag= TRUE, 
                                    id_tip_trim_start= 3,
                                   id_tip_trim_end=7,
                                    id_tip_prefix='',
                                    debug_mode=FALSE, 
                                   debug_print_data_tree= TRUE,
                                   man_multiply_second_legend= 1,
                                   man_multiply_second_legend_text=1.7,
                                   man_multiply_first_legend_text= 1,
                                   man_multiply_first_legend_title_size=1,
                                   man_space_second_legend_multiplier=1,
                                   man_offset_for_highlight_legend_x=0,
                                   a_4_output=FALSE) {
  
    yaml_file<- func.read_yaml(conf_yaml_path)
    
    #get configuration from yaml file
    
    

    
    out_list<-c()
    out_index <-1
    
    if (debug_mode== TRUE) {
        print("##DEBUG##")
        print("yaml is")
        print(yaml_file)
        
    }
    units_out <- 'cm'
    if (a_4_output== TRUE) {
        print("output in A4 format")
        width <- 297
        height <- 210
        units_out <-  "mm"
        man_multiply_second_legend_text <- 0.34
        man_multiply_first_legend_title_size <- 0.25
        
    }

    individual<- yaml_file[['Individual general definitions']]$Individual
    
    if (is.na(individual)) {
        stop("Missing individual name")
    }

    csv_path<- yaml_file[['Individual general definitions']]$'mapping csv file'
    
    if(is.na(csv_path)) {
        stop("Missing csv file")
    }
    
    #get csv file for mapping subgroups
    print(paste0("Get mapping csv from: ",csv_path))
    readfile <- read.csv(csv_path)
    


    tree_path_list<- yaml_file[['Individual general definitions']]$'tree path'

    title.id <- yaml_file[['Mapping exl renaming titles']]$'ID column'
    
    if (debug_mode== TRUE) {
        print("title.id")
        print(title.id)
        
    }
    
    
    trees_number <- length(tree_path_list)
    if (trees_number==0) {
        stop("No tree path was given")
    }

    flag_short_tips <- yaml_file[['visual definitions']]$'trim tips'$'display'
    tips_length <- as.numeric(yaml_file[['visual definitions']]$'trim tips'$'length')
    path_base <- yaml_file[['Individual general definitions']]$out_file$base_path
    
    
    edge_width_multiplier<-yaml_file[['visual definitions']]$'edge_width_multiplier'$'size'
    if (is.na(edge_width_multiplier)) {
        print("Missing edge_width_multiplier in yaml file")
        print("setting to default value")
        edge_width_multiplier<-1
    }
    size_tip_text<-yaml_file[['visual definitions']]$'font_size'$'tips'
    if (is.na(size_tip_text)) {
        print("Missing size_tip_text in yaml file")
        print("setting to default value")
        size_tip_text <- 10
    }
    
    size_font_legend_title <- yaml_file[['visual definitions']]$'font_size'$'legend_title'
    if(is.na(size_font_legend_title)) {
        print("Missig font size for legend title")
        print("setting to default value")
        size_font_legend_title <- 50
    }
    size_font_legend_text <- yaml_file[['visual definitions']]$'font_size'$'legend_text' 
    if (is.na(size_font_legend_text)) {
        print("Missing font size for legend text")
        print("setting to default value")
        size_font_legend_text <- 35
    }
    size_font_legend_box <- yaml_file[['visual definitions']]$'font_size'$'legend_box' 
    if (is.na(size_font_legend_box)) {
        print("Missing font size for lengend box")
        print("setting to ddefault value")
        size_font_legend_box <- 30
    }
    
    size_font_heat_map_text <- yaml_file[['visual definitions']]$'font_size'$'heat_map_title' 
    if (is.na(size_font_legend_box)) {
        size_font_heat_map_text <- 25
    }

    
    for (tree_index in 1:trees_number) {
        print(paste0('Tree number ',tree_index))
        tree_path <- tree_path_list[tree_index]
    
        ###get newick file and create tree
        print(paste0("Get tree from: ",tree_path))
        tree440 <- read.tree(tree_path)
        
        print(paste0("Individual name is ",individual))

        if (debug_mode==TRUE) {
            print("individual is")
            print(individual)
        }
        if ("individual" %in% names(readfile)) {
            idv_title <- "individual"
        } else {
            idv_title <- "Individual"
        }

        #make subframe for the data of specific individual
        if (idv_title %in% names(readfile)) {   
            readfile440 <- readfile[readfile[[idv_title]]==individual,]
        } else {
            readfile440 <- readfile
        }
        

    
        if (debug_mode== TRUE) {
            print("##DEBUG##")
            print("readfile440 is")
            print(readfile440)
        }
        
        ids_list<-readfile440[[title.id]]
        
        if (debug_mode== TRUE) {
            print(readfile440)
            print("title.id is")
            print(title.id)
            print("ids_list from csv is")
            print(ids_list)
        }
        
        readfile440 <- fix.readfile440.with.missing.leaves(readfile440,title.id,tree440,ids_list,"na",
                                                           id_tip_trim_flag,id_tip_trim_start,id_tip_trim_end)
        
        if (debug_mode==TRUE){
            print("readfile440 is is")
            print(readfile440)
        }
        

        
        keys_visual <- names(yaml_file[['visual definitions']])
        disp_opt_num <- length(yaml_file[['visual definitions']]$'classification')
        if (is.na(disp_opt_num)){
            stop("Missing classification in visual definitions")
        }
        if (disp_opt_num==0) {
            stop("Missing classification in visual definitions")
        }
        
        
        rot1<-yaml_file[['visual definitions']]$'rotation1'$'display'
        rot2<-yaml_file[['visual definitions']]$'rotation2'$'display'
    
        if (debug_mode== TRUE) {
            print("##DEBUG##")
            print("rot1 is")
            print(rot1)
            print("rot2 is")
            print(rot2)
        
        }
        
  
        
  
        list_rotate <- func.calc.rotation.definitions(rot1,rot2,yaml_file,title.id,ids_list,tree440,readfile440,debug_mode)
        

        
        rotate_flag_for_title <-list_rotate$rotate_flag_for_title
        rotate_str <- list_rotate$rotate_str
        rotation_params1 <- list_rotate$rotation_params1
        rotation_params2 <- list_rotate$rotation_params2
        rotate_flag <- list_rotate$rotate_flag
        rotate1_types_list_dx.rx <- list_rotate$rotate1_types_list_dx.rx
        list_weight_frac1 <- list_rotate$list_weight_frac1
        
        

        
        #highlight used to be here, moved inside
        
        show_boot <- yaml_file[['visual definitions']]$'Bootstrap'$'display'
        if (func.check.bin.val.from.conf(show_boot)== TRUE) {
            show_boot_flag <- TRUE
        } else {
            show_boot_flag <- FALSE
        }
        
        if (! 'classification' %in% keys_visual) {
            stop("Missing classification in yaml")
        }

        
        for (disp_index in 1:disp_opt_num) {
            print("disp_index is")
            print(disp_index)
            disp_index_str <- as.character(disp_index)
            temp <- yaml_file[['visual definitions']]$'classification'[[disp_index]]
            
            att<- names(yaml_file[['visual definitions']]$'classification'[[disp_index]])


            ttt <- yaml_file[['visual definitions']]$'classification'[[disp_index]]
            ind1<- as.character(disp_index)
            att1 <- names(ttt[[ind1]])
            heat_map_title <- ""
            
            disp_indx_ch <- as.character(disp_index)
            
       
            
            FDR_perc <- yaml_file[['visual definitions']]$'classification'[[disp_index]][[disp_indx_ch]]$'FDR_perc'

        
            no_name <-  yaml_file[['visual definitions']]$'classification'[[disp_index]][[disp_indx_ch]]$'non_cluster_title'
        
            no_name_color <- yaml_file[['visual definitions']]$'classification'[[disp_index]][[disp_indx_ch]]$'non_cluster_color'
        
        
            labels_not_in_legend <- yaml_file[['visual definitions']]$'classification'[[disp_index]][[disp_indx_ch]]$'not_in_legend'
        
            title_replace <- yaml_file[['visual definitions']]$'classification'[[disp_index]][[disp_indx_ch]]$'title'
            
          
            
            dxdf440_for_heat <- NA
            heat_flag <- FALSE
            heat_map_title_list <-c()
           
   
            heat_display_vec=c()
            heat_display_params_list <- c()
            
            if ('heatmap_display' %in% att1) {
                print("IN heat")
                
                heat_definitions <- yaml_file[['visual definitions']]$'classification'[[disp_index]][[disp_indx_ch]]$heatmap_display
                heat_list_len <- length(heat_definitions)
                dxdf440_for_heat <- c()
                
                 print("heat_list_len is")
                print(heat_list_len)

                for (inx in 1:heat_list_len) {
      
                    
                    ind <- as.character(inx)
                    heat_map_i_def <- heat_definitions[[inx]][[ind]]

                    att2 <- heat_map_i_def$display
                    print("att2 is")
                    print(att2)
                    
                    param<- c()
                    if (func.check.bin.val.from.conf(att2)== TRUE) {
                        print("In  check")
                        
                        if ('is_discrete' %in% names(heat_map_i_def)) {
                             param['is_discrete'] <- func.check.bin.val.from.conf(heat_map_i_def$is_discrete)
                        } else {
                            param['is_discrete'] <- FALSE
                        } 
                        
                        if (param['is_discrete']== TRUE) {
                            if ('man_define_colors' %in% names(heat_map_i_def)) {
                                    param['man_define_colors'] <- func.check.bin.val.from.conf(heat_map_i_def$man_define_colors)
                            } else {
                                    param['man_define_colors'] <- FALSE
                            }
                            
                            
                            if ('color_scale_option' %in% names(heat_map_i_def)) {
                                param['color_scale_option'] <- 'A'
                                param[['color_scale_option1']] <- heat_map_i_def['color_scale_option']
                            } else {
                                param['color_scale_option'] <- 'B'
                            }
                            if ('color_scale_range_start' %in% names(heat_map_i_def)) {
                                param['man'] <- TRUE
                                param['color_scale_range_start'] <- as.numeric(heat_map_i_def$color_scale_range_start)
                            } else {
                                 param['man'] <- FALSE
                            }
                            if ('color_scale_range_end' %in% names(heat_map_i_def)) { 
                                param['color_scale_range_end'] <- as.numeric(heat_map_i_def$color_scale_range_end)
                            } else {
                                param['color_scale_range_end'] <- 300
                            }
                            
                        } else {
                            param['low'] <- "beige"
                            param['mid'] <- "seashell2"
                            param['high'] <- "firebrick4"
                            param['midpoint']<- .02
                            if ('low' %in% names(heat_map_i_def)) {
                                param['low'] <-heat_map_i_def$low
                            }
                            if ('mid' %in% names(heat_map_i_def)) {
                                param['mid'] <-heat_map_i_def$mid
                            }
                            if ('high' %in% names(heat_map_i_def)) {
                                param['high'] <-heat_map_i_def$high
                            }
                            if ('midpoint' %in% names(heat_map_i_def)) {
                                param['midpoint'] <-as.numeric(heat_map_i_def$midpoint)
                            }
                        }
  

         
                        heat_display_params_list[[inx]] <- param
                        
                        
                        
                        heat_display_vec <- c(heat_display_vec, TRUE)
                        heat_flag <- TRUE
                        acc_heat_list <- heat_map_i_def$according
                        heat_map_title <- heat_map_i_def$title
                        heat_map_title_list <- c(heat_map_title_list,heat_map_title)
                        l_titles_for_heat <- c()
                        ind <-1
                        for (j in acc_heat_list) {

                            j1 <- names(j)

                            ind<- ind+1
                            j2<- j[[j1]]

                            l_titles_for_heat <- c(l_titles_for_heat,j2)
                        }
                        
                        
              
                        df_heat_temp <-readfile440[, c(title.id,l_titles_for_heat)]
                  
                        g_check <- ggtree(tree440)$data
                        g_check_tip <- subset(g_check,isTip==TRUE)
                        if (id_tip_trim_flag== TRUE) {
                             tip_list <- substr(g_check_tip$label,id_tip_trim_start, id_tip_trim_end)
                        } else {
                            tip_list <- g_check_tip$label
                        }
                        
                        
                        df_heat_temp <- df_heat_temp[match(tip_list, df_heat_temp[[title.id]]),]
               
                        
                         if (id_tip_trim_flag== TRUE) {
                              rownames(df_heat_temp) <- paste0(id_tip_prefix,df_heat_temp[[title.id]])
                             } else {
                             rownames(df_heat_temp) <- paste0("",df_heat_temp[[title.id]])
                         }
                       
           
                        if (length(l_titles_for_heat)>1) {
                        
                            df_heat_temp <- df_heat_temp[, c(l_titles_for_heat)]
                        } else {
                       
                            df_heat_temp <- df_heat_temp[, c(l_titles_for_heat,l_titles_for_heat)]

                            df_heat_temp <- df_heat_temp[ -c(2) ]
                           
                        }
                
                        
                        dxdf440_for_heat[[inx]] <- df_heat_temp  
                       
                    } else {
                        heat_display_vec <- c(heat_display_vec, FALSE)
                    }
                    temp <- dxdf440_for_heat[[inx]]
                    
                    
                    
                    if ("with" %in% names(heat_map_i_def)) {
                        
                        with_title <- heat_map_i_def$with$title
                        with_value <- heat_map_i_def$with$value
                        col_with <- readfile440[, with_title]
                        print("with_title is")
                        print(with_title)
                        print("with_value is")
                        print(with_value)
                        print("col is")
                        print(col_with)
                        
                        print("colnames(temp) is")
                        print(colnames(temp))
                  
                        for (col_index in colnames(temp)) {
                            
                            old_col <- temp[,col_index]
                            print("old_col is")
                            print(old_col)
                            print("length(old_col) is")
                            print(length(old_col))
                            print("col_with is")
                            print(col_with)
                            print("length(col_with) is")
                            print(length(col_with))
                            if (length(col_with) < length(old_col) ) {
                                gap <- length(old_col) - length(col_with)
                                col_with <- c(col_with, rep(NA,gap))
                            }
                        
                            new_col <-c()
                            for (row_index in 1:length(old_col)) {
                                val <- col_with[[row_index]]
                                if (func.check.if.id.in.sub.class(val, with_title, with_value) == TRUE) {
                                    new_col <- c(new_col, old_col[[row_index]])
                                } else {
                                    new_col <- c(new_col,NA)
                                }
                            }
                        temp[,col_index] <- new_col
                        
                        }
                    dxdf440_for_heat[[inx]] <- temp

                    }

                }
            
            }
            
            #print("dxdf440_for_heat is")
            #print(dxdf440_for_heat)
            
            how_many_hi <- 0
        
            
            FLAG_BULK_DISPLAY <- FALSE
    
            if ('highlight' %in% att1) {
                hi_def <- yaml_file[['visual definitions']]$'classification'[[disp_index]][[disp_indx_ch]]$highlight
                hi<- hi_def$display
 
                if (func.check.bin.val.from.conf(hi)== FALSE) {
                    #no highlight
                } else {
                    FLAG_BULK_DISPLAY <- TRUE
    
                    highlight.params.NEW <- func.make.highlight.params.NEW(yaml_file,title.id,ids_list,tree440,readfile440,hi_def,debug_mode)

                    how_many_hi <- highlight.params.NEW$how_many_hi
                    high_label_list<- highlight.params.NEW$high_label_list
                    high_color_list<- highlight.params.NEW$high_color_list
                    high_title_list<- highlight.params.NEW$high_title_list
                    lists_list_hi<- highlight.params.NEW$lists_list_hi

                    
                }

            }
            
          

            ######
            according_list <- yaml_file[['visual definitions']]$'classification'[[disp_index]][[disp_indx_ch]]$'according'
            cls_list <- c()
            cls_renaming_list <- c()
            colors_scale1 <- c()

            for (j in 1:length(according_list)) {
                acc <- according_list[[j]]
  
                acc1 <- acc[[as.character(j)]]

                j_ch <- as.character(j)
                cls_list <- c(cls_list, j_ch)
                
                nams <- names(acc[[j_ch]])
                for (cri in c("color","display_name","title1","value1")) {
                    if (cri %in% nams) {
                        
                    } else {
                        print(paste0(paste0(paste0("Missing field in according ",acc),": "),cri))
                    }
                }
                          

                color <- acc[[j_ch]]$color

                
                display_name <- acc[[j_ch]]$display_name
                cls_renaming_list <- c(cls_renaming_list, display_name)
                
                
            }
            

          
     
            cls_num <- length(cls_list)

            
            if (debug_mode== TRUE) {
                print("##DEBUG##")
                print("disp_index is")
                print(disp_index)
                print("cls_list=")
                print(cls_list)

        
            }
            
            tree_data <- ggtree(tree440)$data
            leaves_id_from_tree1 <- tree_data[tree_data$isTip==TRUE,'label']
            leaves_id_from_tree <- as.list(leaves_id_from_tree1['label'])$label
            

  
            
            list_id_by_class<- func.make.list_id_by_class(
                cls_num, cls_renaming_list,yaml_file,title.id,leaves_id_from_tree,readfile440,
                according_list,debug_mode,id_tip_trim_flag,id_tip_trim_start,id_tip_trim_end)
  
            
            dx_rx_types1_short <- names(list_id_by_class)
            
    
            for (j in 1:length(according_list)) {
                acc <- according_list[[j]]
  

                j_ch <- as.character(j)


                display_name <- as.character(acc[[j_ch]]$display_name)


                color <- acc[[j_ch]]$color   
                
                if (length(list_id_by_class[[display_name]])>1 ) {
                    colors_scale1 <- c(colors_scale1, color)
                }

            }
            
            print("LENN")
            #print(max(nchar(list_id_by_class)))
            #print(list_id_by_class)
            #print(ids_list)

            
            non_cluster_color <- yaml_file[['visual definitions']]$'classification'[[disp_index]][[disp_indx_ch]]$'non_cluster_color'
            colors_scale1 <- c(colors_scale1,non_cluster_color)
            
            #print(yaml_file[['visual definitions']]$'classification'[[disp_index]])
            
            na_name <- yaml_file[['visual definitions']]$'classification'[[disp_index]][[disp_indx_ch]]$'na_name'
            if (is.null(na_name)) {
                na_name <-"NA"
            }
            
            cls <- function.create.cls.list(leaves_id_from_tree,dx_rx_types1_short,list_id_by_class,na_name)

            
            if (na_name %in% dx_rx_types1_short) {
                
            } else {
                dx_rx_types1_short <- c(dx_rx_types1_short,na_name)
                
                colors_scale1 <- c(colors_scale1,no_name_color)
            }
            

            if (id_tip_trim_flag == TRUE) {
                leaves_id_from_tree_num <- as.numeric(substring(leaves_id_from_tree,id_tip_trim_start))
                
            } else {
                leaves_id_from_tree_num <- leaves_id_from_tree
                
            }
            
    
            
            readfile440 <- readfile440[readfile440[[title.id]] %in% leaves_id_from_tree_num, ]
            
       
            
            co_nam <- colnames(readfile440)


            index_id <- which(co_nam==title.id)
            #print("index_id is")
            #print(index_id)

            colnames(readfile440)[index_id] <- "Sample.Reads.ID"
            
        
            readfile440<- distinct(readfile440, Sample.Reads.ID, .keep_all = TRUE)
            colnames(readfile440)[index_id] <- title.id


            dxdf440_dataf<-readfile440[, c(title.id,title.id)]

            
            leaves_id_ordered_for_df440 <- func.make.leaves_id_ordered_for_df440(leaves_id_from_tree1,dxdf440_dataf,title.id,
                                                                                id_tip_trim_flag,id_tip_prefix)
            
 
            
            cls2 <- function.create.cls.list(leaves_id_ordered_for_df440,dx_rx_types1_short,list_id_by_class,na_name)
            
   
            
            dxdf440_dataf['Mapping'] <- cls2 #cls

            dxdf440_dataf = subset(dxdf440_dataf, select = c(title.id,'Mapping') )
            
  

            
            if (laderize_flag == TRUE){
                rotate_flag_str <- paste0(rotate_flag_for_title,"_Ladderized")
            } else {
                rotate_flag_str <- rotate_flag_for_title
            }
            
            out_start <- yaml_file[['Individual general definitions']]$out_file$'optional text at beggining'
            out_end <- yaml_file[['Individual general definitions']]$out_file$'optional text at end'
            file_type <- yaml_file[['Individual general definitions']]$out_file$'file_type'
        
            if (flag_short_tips == TRUE) {
                tips_string <- "tips_normalized"
            } else {
                tips_string <- ""
            }
            
            classication2<- ""
            for (a in dx_rx_types1_short) {
                classication2 <- paste0(paste0(classication2,"_"),a)
            }
            
            
            file_name_end <- paste(out_start,classication2,tips_string,rotate_str, sep="_")
        

        
            if (add_date_to_text_flag ==TRUE) {
                str_date <- gsub(" " , "_", Sys.time()) 
                str_date <- gsub(":" , "_", str_date) 
                file_name_end <- paste(file_name_end,'_',str_date, sep = '')
            }
            
                 
            out_file_path <- paste0(path_base,individual)
        

            path_split <- str_split(tree_path,'/')
            path_split <- path_split[[1]]
        
        
            s_nam1 <- path_split[length(path_split)]

            s_nam<- substr(s_nam1,1,nchar(s_nam1)-7)

        
            out_file_path <- paste0(out_file_path,paste0(paste0('_for_',s_nam),'__'))
        
            flag_replace_name <- yaml_file[['Individual general definitions']]$out_file$'replace name'$flag
        
            if (func.check.bin.val.from.conf(flag_replace_name) == TRUE) {
                out_file_path <- yaml_file[['Individual general definitions']]$out_file$'replace name'$name
                out_file_path<- paste0(path_base,out_file_path)
                out_file_path <- paste0(out_file_path,out_index)
                out_file_path <- paste(out_file_path,'.',file_type, sep = '')
            
            }
            
            
            
            ooo <- str_split(out_file_path," ")
            
            
        
            oooo <- ooo[[1]]
            ind <-1
            for (letter in oooo) {
                if (letter %in% c(".","+","-","#")) {
                    oooo[ind] <- ''
                }
                ind <- ind+1
            }

        
            ooo1 <- ""
            for (ind in 1:length(oooo)) {
 
                part <- oooo[ind]


                ooo1 <- paste0(ooo1,part)

            }
            out_file_path <- ooo1
        
  
        
        #if (FLAG_BULK_DISPLAY== TRUE) {
        #    print("A")
        #    man_space_second_legend_multiplier<- man_space_second_legend_multiplier
          

        #    man_space_second_legend<- man_space_second_legend*3

        #}
        
        
        
            out_file_path <- paste(out_file_path,'.',file_type, sep = '')
        
        
            print("Output file path is")
            print(out_file_path)    
            
            
            
            
           
            #if (heat_flag == FALSE) {
                
             #   func.make.plot.tree.NEW(
             #       tree440 = tree440,
             #       dx_rx_types1_short = dx_rx_types1_short,
             #       #classication =classication,
             #       list_id_by_class= list_id_by_class,
             #       dxdf440_dataf= dxdf440_dataf,
             #       title.id=title.id,
             #       FDR_perc=FDR_perc,
             #       no_name= no_name,
             #       rotate_flag = rotate_flag,
             #       rotation_params1= rotation_params1,
             #       rotation_params2= rotation_params2,
             #       flag_short_tips = flag_short_tips, 
             #       tips_length = tips_length,
             #       show_boot_flag = show_boot_flag,
             #       FLAG_BULK_DISPLAY = FLAG_BULK_DISPLAY,
             #       how_many_hi = how_many_hi,
             #       high_label_list = high_label_list,
             #       high_color_list = high_color_list,
             #       high_title_list = high_title_list,
             #       lists_list_hi = lists_list_hi,
             #       simulate.p.value= simulate.p.value,
             #       width= width,
             #       height= height,
             #       colors_scale1 = colors_scale1,
             #       out_file_path = out_file_path,
             #       edge_width_multiplier = edge_width_multiplier,
             #       size_tip_text= size_tip_text,
             #       size_font_legend_title= size_font_legend_title,
             #       size_font_legend_text = size_font_legend_text,
             #       size_font_legend_box= size_font_legend_box,
             #       labels_not_in_legend= labels_not_in_legend,
             #       no_name_color="black",
             #       debug_mode= debug_mode,
             #       man_adjust_elipse= man_adjust_elipse,
             #       man_multiply_elipse= man_multiply_elipse,
             #       man_adj_second_legend= man_adj_second_legend,
             #       man_space_second_legend= -0.02,
             #       laderize_flag=FALSE, 
             #       cls_renaming_list <- cls_renaming_list,
             #       title_flag = TRUE,
             #       title_replace = title_replace,
             #       flag_classification_format=FALSE,
             #       id_tip_trim_flag, 
             #       id_tip_trim_start,
             #       id_tip_trim_end,
             #       units_out
             #       )
            #} else {
                
                func.make.plot.tree.heat.NEW(
                    tree440 = tree440,
                    dx_rx_types1_short = dx_rx_types1_short,
                    #classication =classication,
                    list_id_by_class= list_id_by_class,
                    dxdf440_dataf= dxdf440_dataf,
                    title.id=title.id,
                    FDR_perc=FDR_perc,
                    no_name= no_name,
                    rotate_flag = rotate_flag,
                    rotation_params1= rotation_params1,
                    rotation_params2= rotation_params2,
                    flag_short_tips = flag_short_tips, 
                    tips_length = tips_length,
                    show_boot_flag = show_boot_flag,
                    FLAG_BULK_DISPLAY = FLAG_BULK_DISPLAY,
                    how_many_hi = how_many_hi,
                    high_label_list = high_label_list,
                    high_color_list = high_color_list,
                    high_title_list = high_title_list,
                    lists_list_hi = lists_list_hi,
                    simulate.p.value= simulate.p.value,
                    width= width,
                    height= height,
                    colors_scale1 = colors_scale1,
                    out_file_path = out_file_path,
                    edge_width_multiplier = edge_width_multiplier,
                    size_tip_text= size_tip_text,
                    size_font_legend_title= size_font_legend_title,
                    size_font_legend_text = size_font_legend_text,
                    size_font_legend_box= size_font_legend_box,
                    labels_not_in_legend= labels_not_in_legend,
                    no_name_color="black",
                    debug_mode= debug_mode,
                    man_adjust_elipse= man_adjust_elipse,
                    man_multiply_elipse= man_multiply_elipse,
                    man_adj_second_legend= man_adj_second_legend,
                    man_space_second_legend= -0.02,
                    laderize_flag=FALSE, 
                    cls_renaming_list <- cls_renaming_list,
                    title_flag = TRUE,
                    title_replace = title_replace,
                    flag_classification_format=FALSE,
                    heat_flag= heat_flag,
                    dxdf440_for_heat = dxdf440_for_heat,
                    heat_map_title_list= heat_map_title_list,
                    man_adjust_image_of_second_legend = man_adjust_image_of_second_legend,
                    man_adj_heat_loc= man_adj_heat_loc,
                    man_boot_x_offset= man_boot_x_offset,
                    man_adj_heat_loc2= man_adj_heat_loc2,
                    man_adj_heat_loc3= man_adj_heat_loc3,
                    id_tip_trim_flag, 
                    id_tip_trim_start,
                    id_tip_trim_end, 
                    debug_print_data_tree,
                    heat_display_vec,
                    heat_display_params_list,
                    man_multiply_second_legend,
                    man_multiply_second_legend_text,
                    size_font_heat_map_text,
                    man_multiply_first_legend_text,
                    man_multiply_first_legend_title_size,
                    man_space_second_legend_multiplier,
                    man_offset_for_highlight_legend_x,
                    units_out
                    )
           # }
            
        out_list[[out_index]] <- out_file_path
        out_index <- out_index+1   
            
        #close for of display options   
        }

        
        
    # close for on trees structure     
    }
    out_file_path <- out_file_path

#close func  
}

func.make.leaves_id_ordered_for_df440 <- function(leaves_id_from_tree1,dxdf440_dataf,title.id,id_tip_trim_flag,id_tip_prefix) {
    leaves_id_ordered_for_df440 <- c()
    base.list <- dxdf440_dataf[[title.id]]

    tips_list <- as.list(leaves_id_from_tree1)

    for (id in base.list) {
        if (id_tip_trim_flag== "TRUE") {
            temp <- paste0(id_tip_prefix,id)
        } else {
            temp <- id
        }
        

        
        if (temp %in% tips_list$label) {
 
            leaves_id_ordered_for_df440 <- c(leaves_id_ordered_for_df440,temp)
        }
    }
    

    leaves_id_ordered_for_df440 <- leaves_id_ordered_for_df440
}

func.make.plot.tree.heat.NEW <- function( tree440,
                    dx_rx_types1_short,
                    list_id_by_class,
                    dxdf440_dataf,
                    title.id,
                    FDR_perc,
                    no_name,
                    rotate_flag,
                    rotation_params1,
                    rotation_params2,
                    flag_short_tips , 
                    tips_length ,
                    show_boot_flag ,
                    FLAG_BULK_DISPLAY ,
                    how_many_hi =0 ,
                    high_label_list ,
                    high_color_list ,
                    high_title_list ,
                    lists_list_hi ,
                    simulate.p.value,
                    width,
                    height,
                    colors_scale1 ,
                    out_file_path ,
                    edge_width_multiplier=1 ,
                    size_tip_text=3,
                    size_font_legend_title=30,
                    size_font_legend_text=20 ,
                    size_font_legend_box=15,
                    labels_not_in_legend,
                    no_name_color,
                    debug_mode,
                    man_adjust_elipse=0,
                    man_multiply_elipse,
                    man_adj_second_legend= 0.15,
                    man_space_second_legend= -0.02,
                    laderize_flag=FALSE, 
                    cls_renaming_list ,
                    title_flag ,
                    title_replace ,
                    flag_classification_format,
                    heat_flag= FALSE,
                    dxdf440_for_heat ,
                    heat_map_title_list= NA,                 
                    man_adjust_image_of_second_legend,
                    man_adj_heat_loc,
                    man_boot_x_offset,
                    man_adj_heat_loc2,
                    man_adj_heat_loc3,
                    id_tip_trim_flag, 
                    id_tip_trim_start,
                    id_tip_trim_end,
                    debug_print_data_tree,
                    heat_display_vec,
                    heat_display_params_list,
                    man_multiply_second_legend,
                    man_multiply_second_legend_text,
                    size_font_heat_map_text,
                    man_multiply_first_legend_text,
                    man_multiply_first_legend_title_size,
                    man_space_second_legend_multiplier,
                    man_offset_for_highlight_legend_x,
                                        units_out) {
    

    
    if (debug_mode==TRUE){
        print("In func.make.plot.tree.HEAT")
    }
    
    
    dx_rx_types1 <- dx_rx_types1_short
    
    pr440 <- ggtree(tree440) 
    d440 <- pr440$data
    cc_tipss <- func.create.cc_tipss(d440)
    cc_nodss <- func.create.cc_nodss(d440)
    cc_totss <- func.create.cc_totss(d440)
    
    nods_num <- length(cc_nodss)
    tips_num <- length(cc_tipss)
    tree_size <- nods_num +tips_num 
    
    if (debug_mode==TRUE){
        print("tree data is")
        print(d440)
    }
    
    if (debug_print_data_tree== TRUE) {
        print(tree_size)
        print(d440)
    }


    
    #####

    s <- subset(d440,isTip=="TRUE")
        
    
    

    
    list_node_by_class <- func.create.list_node_by_class(dx_rx_types1_short,
                                                     list_id_by_class,dxdf440_dataf,
                                                         title.id,tree_size,d440,cc_totss,debug_mode,
                                                        id_tip_trim_flag,id_tip_prefix)
    
    

    list_rename_by_class <-list_node_by_class
    
    
    y_off_base= -8
    
    yet_another_multiplier<- 0 

    
    
    tree_with_group <- groupOTU(tree440, list_node_by_class)

    
    
        
    subframe_of_nodes <- d440[d440$isTip=="FALSE",]
    cc_nodss90 <- func.create.cc_nodss90(subframe_of_nodes)
    cc_nodss80 <- func.create.cc_nodss80(subframe_of_nodes)
    cc_nodss70 <- func.create.cc_nodss70(subframe_of_nodes)
    
    
    
    
    #create array of number of nodes of each type for the whole tree
    df_count_FULL_tree_populations <- data.frame(idx = 1:length(dx_rx_types1_short),
                type = dx_rx_types1_short,
                count = rep(0, length(dx_rx_types1_short))
                )
    
    
    
    for (opt in dx_rx_types1){
        indx <- which(dx_rx_types1_short == opt)
        df_count_FULL_tree_populations[indx,'count'] <- length(list_id_by_class[[opt]])-1
    }
    

    

    
    tree_with_group_CPY <- groupOTU(tree440, list_rename_by_class)
    #debug_tree <- tree_with_group_CPY
    
    
    levels_groups <- levels(ggtree(tree_with_group_CPY)$data$group)
    
    
    
   
    
    tree_TRY <- ggtree(tree_with_group_CPY,aes(color=new_class, size=p_val_new),ladderize=laderize_flag,
                  ) #size=p_val
    
  test_fig <- tree_TRY
    
    
     new_colors_list<- func.create.new_colors_list(FDR_perc,tree_TRY,tree_with_group,no_name,tree_size)
    
    if (no_name %in% new_colors_list) {
        flag_no_name <- 1
    } else {
        flag_no_name <- 0
    }
    
    tree_TRY$data$new_class <- new_colors_list


    
    op_list <- c(paste0("p>",FDR_perc ), 
             paste0(paste0(0.5*FDR_perc,"<p<="),FDR_perc),
             paste0("p<=",0.5*FDR_perc ))
    
    p_list_of_pairs<- func.create.p_list_of_pairs(list_node_by_class,d440,dx_rx_types1_short,
                                            cc_nodss,tree_with_group,FDR_perc,tree, cc_tipss,
                                            tree_TRY,tree_size,no_name,simulate.p.value)
    

    
    p_PAIRS_pval_list <- func.create.p_val_list_FROM_LIST(FDR_perc,tree_TRY,p_list_of_pairs,op_list)
    
     

    
    tree_TRY$data$p_val_new <-  p_PAIRS_pval_list 

    tree_TRY$data$p_val_new <- factor(tree_TRY$data$p_val_new,                 # Relevel group factor
                         levels = op_list)
    
    ######################
    
    ##########################

    #debug_tree <- tree_TRY
    
    if (FLAG_BULK_DISPLAY== TRUE){

        
        for (index_high in 1:how_many_hi) {
            list_high_for <- lists_list_hi[[index_high]]

            if (index_high==1) {
                tree_TRY$data$'high1' <-  list_high_for 
                #tree_TRY$data$high1 <- factor(tree_TRY$data$high1,                 # Relevel group factor
                 #        levels = c(TRUE,FALSE))
            } else {
                if (index_high==2) {
                   tree_TRY$data$'high2' <-  list_high_for 
                    #tree_TRY$data$high2 <- factor(tree_TRY$data$high2,                 # Relevel group factor
                    #     levels = c(TRUE,FALSE)) 
                } else {
                    if (index_high==3) {
                   tree_TRY$data$'high3' <-  list_high_for 
                #    tree_TRY$data$high3 <- factor(tree_TRY$data$high3,                 # Relevel group factor
                #         levels = c(TRUE,FALSE)) 
                    } else {
                        stop ("Too many highlight options. Only 3 are supported")
                    }
                }
            
            
            }
    
        }

    #tree_TRY$data$high_id <-  list_high_with_id 
    #tree_TRY$data$high_id <- factor(tree_TRY$data$high_id,                 # Relevel group factor
    #                     levels = unique(list_high_with_id))
    }
    
    tree_TRY1<- tree_TRY

    
    
       #rotate tree if required
    
    if (rotate_flag %in% c("RX_first", "FRAC_first") ) {
 
        list_weight_dx.rx <- rotation_params1[['list_weight_dx.rx']]
        list_weight_frac <- rotation_params2[['list_weight_dx.rx']]
        TREE_OTU_dx.rx  <- rotation_params1[['TREE_OTU_dx.rx']]
        TREE_OTU_frac <- rotation_params2[['TREE_OTU_dx.rx']]
    }
    
    if (rotate_flag =="RX_first"){
        

        

        tree_TRY2 <- func.rotate.tree.based.on.weights (tree_TRY1,list_weight_dx.rx,list_weight_frac,
                                    TREE_OTU_dx.rx,TREE_OTU_frac,tree_size)
   
        
        
        
        } else if (rotate_flag =="FRAC_first"){
            
        
            tree_TRY2 <- func.rotate.tree.based.on.weights (tree_TRY1,list_weight_frac,list_weight_dx.rx,
                                    TREE_OTU_frac,TREE_OTU_dx.rx,tree_size)
        } else {
        
        
            tree_TRY2 <- tree_TRY1
        }
    
    
        temp <- tree_TRY
       #tree with short tips if required
        pr440_short_tips_TRY <- tree_TRY2 #tree_TRY
    
        if (flag_short_tips == TRUE) {
        for(i in cc_tipss) {
            par <- pr440_short_tips_TRY$data$parent[i]
            parX <- pr440_short_tips_TRY$data$x[par]
            pr440_short_tips_TRY$data[pr440_short_tips_TRY$data$node[i], "x"] = parX +tips_length
        }
    
        if (debug_mode== TRUE) {
            print("##DEBUG##")
            print("flag_short_tips is")
            print(flag_short_tips)
        
        }
    } 

    
    
        #make column for bootstrap values only in the tree, used in display
    zs =rep(0, tree_size)
    pr440_short_tips_TRY$data$boot_val <-  zs
    for (i in 1:tree_size) {
 
        if (i %in% cc_tipss) {
            pr440_short_tips_TRY$data$boot_val[i] = 0
        } else {
            pr440_short_tips_TRY$data$boot_val[i] =pr440_short_tips_TRY$data$label[i]
        }
    }
    
    cls_renaming_list_with_no_name <- c(cls_renaming_list,no_name)
    
    
    if (missing(how_many_hi)) {
        how_many_hi <- 0
    }
    
    
    if (flag_no_name== 1) {

        fin_color_list <- cls_renaming_list_with_no_name
    } else {

        fin_color_list <- cls_renaming_list
    }
    

    

    
     pr440_short_tips_TRY$data$new_class <- factor(pr440_short_tips_TRY$data$new_class,                 # Relevel group factor
                         levels = unique(dx_rx_types1_short))
    

    
    
    list_of_sizes <- c(1,  1.6,  2.2)*edge_width_multiplier
    


    pr440_short_tips_TRY_new <- pr440_short_tips_TRY + scale_color_manual(values=colors_scale1) + 
    scale_size_manual(values=list_of_sizes)
    

    
        
    t<- pr440_short_tips_TRY_new

    #Add bootstrap values to the tree if required

   
    #x_range_min = min(pr440_short_tips_TRY_new$data$x)
    #x_range_max = max(pr440_short_tips_TRY_new$data$x)
    
    square_size_legend <- size_font_legend_box # 30  #16
    text_legend_title_size <- size_font_legend_title # 50 #30
    text_legend_size <- size_font_legend_text #35 # 25
    
    #size_90 <- 26
    #size_80 <- 24
    #size_70 <- 22
    
    size_70 <- edge_width_multiplier
    size_80 <- size_70+1
    size_90 <- size_70+2
    
    
    if (show_boot_flag == TRUE) {
        print("man_boot_x_offset is")
        print(man_boot_x_offset)

    pr440_short_tips_TRY_new_with_boot <- pr440_short_tips_TRY_new + 
    geom_nodepoint(position=position_nudge(x = man_boot_x_offset, y = 0),aes(subset = boot_val >= 0.9 ),size = size_90,
               shape=24, fill="grey36",colour = "grey20",show.legend = FALSE,alpha=1/2) + 
    geom_nodepoint(position=position_nudge(x = +man_boot_x_offset, y = 0),aes(subset = boot_val >= 0.8 & boot_val < 0.9 ),size = size_80,
               shape=24, fill="grey36",colour = "grey20",show.legend = FALSE,alpha=1/2) + #orchid3
    geom_nodepoint(position=position_nudge(x = man_boot_x_offset, y = 0),aes(subset = boot_val >= 0.7 & boot_val < 0.8),size = size_70,
               shape=24, fill="grey36", colour = "grey20",show.legend = FALSE,alpha=1/2)  #orchid1 #"#b5e521"
        
    } else {

        pr440_short_tips_TRY_new_with_boot <- pr440_short_tips_TRY_new
    }
    
        
    # add legends to the tree
    group_display_title <- "Cell type"
    if (exists('title_flag')) {
        if (title_flag== TRUE) {
            group_display_title <- title_replace
        }
    } 
    


   
    new_base <-0
    
    pr440_short_tips_TRY_new_with_boot<- pr440_short_tips_TRY_new_with_boot#+geom_treescale(x=-20,y=0) 
    

     
    pr440_short_tips_TRY_new_with_boot_more1 <- pr440_short_tips_TRY_new_with_boot
    
    text_mu <- round(width/100)-2
    if (units_out=="cm") {
        text_mu <- text_mu*10
    }
    mar <- round(width/50)
    mar1 <- mar+5
    print("text_mu is")
    print(text_mu)
 
    #print(pr440_short_tips_TRY_new_with_boot_more1)
 
    pr440_short_tips_TRY_new_with_boot_more1 <- pr440_short_tips_TRY_new_with_boot_more1 + 
    layout_dendrogram() +
    guides(colour = guide_legend(title = group_display_title, byrow = TRUE,
                                 override.aes = list(size = square_size_legend*text_mu, alpha = 1,linetype="solid", linewidth=1),
                                breaks = levels(pr440_short_tips_TRY_new_with_boot_more1$data$group)),
       size = guide_legend(title = "p value")) + 
    theme(legend.text = element_text(size = text_legend_size*text_mu,margin = margin(t = mar, b=mar, unit = "pt")),
      legend.title = element_text(size = text_legend_title_size,face = "bold"),
     legend.spacing.y = unit(text_mu/6, 'cm')) +
  guides(fill = guide_legend(byrow = TRUE))+
    guides(colour = guide_legend(title = group_display_title,byrow = TRUE))+
    scale_y_reverse()  +
    geom_rootedge()
    
    
    
    x_range_min <-  min(pr440_short_tips_TRY_new_with_boot_more1$data$x)
    x_range_max <-  max(pr440_short_tips_TRY_new_with_boot_more1$data$x)
    
    
    
    
    x_off_base<-  x_range_min*0.4  -man_adj_second_legend  #-0.9 #-0.76

    # highlight bullk if required
    
    print("x_range_min is")
    print(x_range_min)
    print("x_range_max is")
    print(x_range_max)
    print("x_off_base is is")
    print(x_off_base)
    
 op <- -3 
        
    pr440_short_tips_TRY_new_with_boot_more1 <- pr440_short_tips_TRY_new_with_boot_more1 +
    geom_tiplab(size =size_tip_text,angle=90,hjust=1,show.legend = FALSE)
    
    how_mant_rows <-0
    how_many_boxes <- 0
    if (FLAG_BULK_DISPLAY== TRUE) {
        how_mant_rows <- how_mant_rows+ how_many_hi
        how_many_boxes <- how_many_boxes+1
    }
    if (show_boot_flag== TRUE) {
        how_mant_rows <- how_mant_rows+3
        how_many_boxes <- how_many_boxes+1
    }
    print("how_mant_rows is")
    print(how_mant_rows)
    print("how_many_boxes is")
    print(how_many_boxes)
    
    if (heat_flag== FALSE) {
            new_base_for_second_legend_non <- 0.5+ 0.5*(how_mant_rows)*0.05+ 0.5*how_many_boxes*0.1
            new_base_for_second_legend_normalized <- new_base_for_second_legend_non*x_range_min
            new_step <- 0.02*x_range_min
            new_big_step <- 0.1*x_range_min
            extra <- 1
        move <- 0
    } else {
        print("AAAAAAAAAAAAAAAAAAAAAAAA")
            new_base_for_second_legend_non <- 4
            new_base_for_second_legend_normalized <- new_base_for_second_legend_non*x_range_min
            new_step <- 0.15*x_range_min
            new_big_step <- 0.2*x_range_min
         extra <- 6
        move <- 0.9
    }

    print("new_base_for_second_legend_non is")
    print(new_base_for_second_legend_non)
  print("new_base_for_second_legend_normalized is")
    print(new_base_for_second_legend_normalized)



        if (FLAG_BULK_DISPLAY ==TRUE) {
        
    #b <- 0.28 #+ ((100-width)*0.005)#0.45 #0.45 + (height-60)*0.04
    #a <-   0.017 + ((100-height)*0.0004) #0.035 + (height-60)*(-0.002) #0.15
    #x_adj_hi <- 0.017  + ((100-height)*0.0001)#0.025
            x_adj_hi <- 0
    if (heat_flag== FALSE) {
            a <- id_tip_trim_end*0.002  #height of elipse
            b <- 0.45 #width of elipse
    } else {
        a <- id_tip_trim_end*0.045
        b <- 0.45 #width of elipse
    }   
    #id_tip_trim_end
            print("a is")
            print(a)
            print("b is")
            print(b)

    #if (x_adj_hi<=0) {
    #    x_adj_hi <-0.01
    #}
    #if (a<=0) {
    #    a<-0.01
    #}
    #if (b<=0) {
    #    b<-0.01
    #}
    
    #b<- 0.45 #0.15 #0.45 /1.5 * (x_range_min *(-1))

    #a<- 0.035 /1.5 * (x_range_min *(-1)) *man_multiply_elipse #0.01 #0.035
        
     #x_adj_hi <- 0 #0.017
        
    up_offset<- -1 #-3
    #a<- 0.5 #0.02
   print("man_adjust_elipse is")         
   print(man_adjust_elipse)
           

            #x_adj_hi
    
    y_off_base<-  -8
    x_off_base<-round(x_range_min/4) #x_range_min +0.1 #x_range_min*0.6 -man_adj_second_legend  #-0.9 #-0.76
     x_off_base2 <- round(height/4)

       print("pr440_short_tips_TRY$data$x is")
            print(pr440_short_tips_TRY$data$x)
            print("max pr440_short_tips_TRY$data$x is")
            print(max(pr440_short_tips_TRY$data$x))
            print("min pr440_short_tips_TRY$data$x is")
            print(min(pr440_short_tips_TRY$data$x))
            print("man_adjust_elipse is")
            print(man_adjust_elipse)
            
        for (index_high in 1:how_many_hi) {
            #print("index_high is")
            #print(index_high)
            if (index_high==1) {
                high_nodes_table1 <- pr440_short_tips_TRY$data[tree_TRY$data$high1 == TRUE,]
                
                            print("high_nodes_table1$x is")
            print(high_nodes_table1$x)
                print("max high_nodes_table1$x is")
                print(max(pr440_short_tips_TRY$data[,'x']))
                print("max(pr440_short_tips_TRY$data[,'x'])- high_nodes_table1$x is")
                print(max(pr440_short_tips_TRY$data[,'x'])- high_nodes_table1$x)
                print("CHECK")
                print(max(pr440_short_tips_TRY$data[,'x']))
                #orig
                #(max(pr440_short_tips_TRY$data[,'x'])-x)*(-16) +man_adjust_elipse,
                
                bas <- max(max(pr440_short_tips_TRY$data[,'x'])- high_nodes_table1$x)
                 if (heat_flag== FALSE) {
                pr440_short_tips_TRY_new_with_boot_more1 <- pr440_short_tips_TRY_new_with_boot_more1 +
            geom_ellipse(data=high_nodes_table1,
                         aes(x0 = ((max(pr440_short_tips_TRY$data[,'x'])-x )*(-1)+man_adjust_elipse), #0.01
                                           y0 = y, a = a, b = b, angle = 0),
                    fill=high_color_list[[1]],alpha=0.5,linetype="blank",show.legend=FALSE)
                     } else {
                     pr440_short_tips_TRY_new_with_boot_more1 <- pr440_short_tips_TRY_new_with_boot_more1 +
            geom_ellipse(data=high_nodes_table1,
                         aes(x0 = ((max(pr440_short_tips_TRY$data[,'x'])-x )*(-15.4)+man_adjust_elipse+0.18), #0.01
                                           y0 = y, a = a, b = b, angle = 0),
                    fill=high_color_list[[1]],alpha=0.5,linetype="blank",show.legend=FALSE)
                 }
                #pr440_short_tips_TRY_new_with_boot_more1 <- pr440_short_tips_TRY_new_with_boot_more1+
                #    geom_nodepoint(position=position_nudge(x = 0, y = 0),aes(subset = high1 == TRUE ),size = 20,
               #shape=10, fill=high_color_list[[1]],colour = high_color_list[[1]],show.legend = FALSE,alpha=1/2) 
            } else {
                if (index_high==2) {
                    high_nodes_table2 <- pr440_short_tips_TRY$data[tree_TRY$data$high2 == TRUE,]
                    pr440_short_tips_TRY_new_with_boot_more1 <- pr440_short_tips_TRY_new_with_boot_more1 +
                    geom_ellipse(data=high_nodes_table2,
                         aes(x0 = ((max(pr440_short_tips_TRY$data[,'x'])-x )*(x_range_min)), 
                                           y0 = y, a = a, b = b, angle = 0),
                        fill=high_color_list[[2]],alpha=0.5,linetype="blank",show.legend=FALSE)
                } else {
                    if (index_high==3) {
                        high_nodes_table3 <- pr440_short_tips_TRY$data[tree_TRY$data$high3 == TRUE,]
                        pr440_short_tips_TRY_new_with_boot_more1 <- pr440_short_tips_TRY_new_with_boot_more1 +
                            geom_ellipse(data=high_nodes_table3,
                             aes(x0 = ((max(pr440_short_tips_TRY$data[,'x'])-x )*(x_range_min)), 
                                           y0 = y, a = a, b = b, angle = 0),
                            fill=high_color_list[[3]],alpha=0.5,linetype="blank",show.legend=FALSE)
                    }
                }
            }
            
        }
                 


            
    p<- pr440_short_tips_TRY_new_with_boot_more1
    p1 <- p$layers[1]
    p2 <- p$layers[2]
    p7 <- p$layers[7]
    p$layers[2] <- p7 
    p$layers[7] <- p2
     
    
    p3 <- p$layers[3]
    p4 <- p$layers[4]
    p5 <- p$layers[5] 
    p6 <- p$layers[6]
    p7 <- p$layers[7]
    p8 <- p$layers[8]
            
    p$layers[3] <- p6
    p$layers[4] <- p7
    p$layers[5] <- p8
    p$layers[6] <- p3
    p$layers[7] <- p4  
    p$layers[8] <- p3
     
            #old
            #x_off_base+ multiple_high_up_offset+multiple_high_down_offset-2 + (index_high-1)*1 
             #+(how_many_hi+7-index_high)*yet_another_multiplier +man_offset_for_highlight_legend_x
            
            #x_off_base+0.035+multiple_high_up_offset+multiple_high_down_offset-0.7 + (index_high-1)*1 
             #+(how_many_hi+7-index_high)*yet_another_multiplier +man_offset_for_highlight_legend_x

     yet_another_multiplier<- 0   
    if (FLAG_BULK_DISPLAY== TRUE) {

    man_space_second_legend<- man_space_second_legend*2
        yet_another_multiplier <- 1.5
    }
            print("yet")
            print(yet_another_multiplier)
            
    stair <- 0.05
    multiple_high_up_offset <- stair*(how_many_hi-1)*(-1)
            print("multiple_high_up_offset is")
            print(multiple_high_up_offset)
            print("new_big_step is")
            print(new_big_step)
            print("new_step is")
            print(new_step)
            
            

            
        for (index_high in 1:how_many_hi) {
            
            multiple_high_down_offset <- (index_high-1)*stair
            if (index_high==1) {
                if (heat_flag== FALSE) {
                    x11 <- new_base_for_second_legend_normalized- 2*new_big_step - 2*new_step
                    x22 <- new_base_for_second_legend_normalized- 2*new_big_step - 2*new_step - index_high*(new_step)    
                    norm <- 1
                } else {
                    x11 <- op
                    x22 <- op - index_high*(new_step*10)
                    norm <- 0.8
                }
                print("ya")
                print(new_base_for_second_legend_normalized+ new_big_step + 2*new_step + new_big_step)
                #print(x_off_base2+multiple_high_up_offset+multiple_high_down_offset+man_offset_for_highlight_legend_x- 
                #    (how_many_hi+2-index_high)*yet_another_multiplier)
                p <- p +
                annotate(
            geom = "text", 
            label = high_title_list[[index_high]],  size = size_font_legend_title*man_multiply_second_legend_text, 
                x = x11,
                y = y_off_base,  fontface = "bold"
              ) 
                
                if (how_many_hi>1) {
                    p<- p+
                    geom_ellipse(aes(x0 = x22,
                         y0 = y_off_base+0.7+(width/400)+man_adjust_image_of_second_legend, 
                         a = a, b = b, angle = 0),
            fill=high_color_list[[index_high]],alpha=0.5,linetype="blank",show.legend=FALSE)
                }

            }
            print("ya ya")
            print(new_base_for_second_legend_normalized+ new_big_step + 2*new_step + new_big_step + index_high*(new_step))
            p <- p +
         annotate(
        geom = "text", 
        label = high_label_list[[index_high]],  size = size_font_legend_text* man_multiply_second_legend* norm, 
             x = x22,
             y = y_off_base-2,  
      )+
        geom_ellipse(aes(x0 = x22,
                         y0 = y_off_base+0.7+(width/400)+man_adjust_image_of_second_legend, 
                         a = a, b = b, angle = 0),
            fill=high_color_list[[index_high]],alpha=0.5,linetype="blank",show.legend=FALSE)
       
          #new_base <- 0.035+multiple_high_up_offset+multiple_high_down_offset  
         #print(index_high)   
            #print(high_color_list[[index_high]])
            
        }
    } else {

        pr440_short_tips_TRY_new_with_boot_more1 <-pr440_short_tips_TRY_new_with_boot_more1 +
        geom_tiplab(size =size_tip_text,angle=90,hjust=1,,show.legend = FALSE)
        p <- pr440_short_tips_TRY_new_with_boot_more1
    }


    
        #add second legend for bootstrap and highlight if required
    
    second_legend_title_font_size<-  round(text_legend_title_size*2,digits = 0)  #16 # 10.5
    second_text_legend_size<- round(text_legend_size*2,digits = 0) #8.6

    #boot_spacing_legend<- 0.1 #0.04+   man_space_second_legend
    
    
       boot_spacing_legend <-0.5 +    man_space_second_legend #0.04+    man_space_second_legend
    offset_x_xtra <- -2

    #second_legend_title_font_size
    base_from_hi <- (how_many_hi-1)*(-0.5)
    print("man_space_second_legend_multiplier")
    print(man_space_second_legend_multiplier)
    print("new_base_for_second_legend_normalized is")
    print(new_base_for_second_legend_normalized)
        if (show_boot_flag == TRUE) {
         
    print("new_base is")
            print(new_base)
            print("x_off_base is")
            print(x_off_base)
            print("offset_x_xtra is")
            print(offset_x_xtra)
            print("base_from_hi is")
            print(base_from_hi)
            print("yet_another_multiplier is")
            print(yet_another_multiplier)
    #new_base_for_second_legend_normalized        
            #old
            #x_off_base+ new_base -offset_x_xtra+base_from_hi -2* yet_another_multiplier
            #x_off_base+ new_base -offset_x_xtra+ boot_spacing_legend+base_from_hi*man_space_second_legend_multiplier
    p<- p+
     annotate(
    geom = "text", 
    label = "Bootstrap",  size = round(size_font_legend_title*man_multiply_second_legend_text), 
         x = new_base_for_second_legend_normalized, 
         y = y_off_base,  fontface = "bold"
      )+ 
 annotate(
    geom = "text", 
    label = ">70%",  size = round(size_font_legend_text*man_multiply_second_legend), 
     x = new_base_for_second_legend_normalized- new_big_step*extra +move, 
     y = y_off_base-2,  
  )+
 annotate(
    geom = "text", 
    label = ">80%",  size = round(size_font_legend_text*man_multiply_second_legend), 
     x = new_base_for_second_legend_normalized- (new_big_step - new_step)*extra +move , 
     y = y_off_base-2,  
  )+
 annotate(
    geom = "text", 
    label = ">90%",  size = round(size_font_legend_text*man_multiply_second_legend), 
     x = new_base_for_second_legend_normalized- (new_big_step - 2*new_step)*extra +move, 
     y = y_off_base-2,  
  ) +
 annotate(
    geom = "point", 
    shape=24,  size = size_90, x = new_base_for_second_legend_normalized- (new_big_step - 2*new_step)*extra +move, 
     y = y_off_base+1+man_adjust_image_of_second_legend*2,  
     fill="grey36",colour = "grey20",alpha=1/2,
  )+
 annotate(
    geom = "point", 
    shape=24,  size = size_80, x = new_base_for_second_legend_normalized- (new_big_step - new_step)*extra +move , 
     y = y_off_base+1+man_adjust_image_of_second_legend*2,  
     fill="grey36",colour = "grey20",alpha=1/2,
  ) +
 annotate(
    geom = "point", 
    shape=24,  size = size_70, 
     x = new_base_for_second_legend_normalized- new_big_step*extra +move,
     y = y_off_base+1+man_adjust_image_of_second_legend*2,  
     fill="grey36",colour = "grey20",alpha=1/2,
  ) 
      print("CHECK CHECK")
            print(new_base_for_second_legend_normalized- (new_big_step - 2*new_step)*extra +move)
    }
    
    
  ##############################
    
   if (heat_flag == TRUE) { 
       
    tt<- p
    for(i in cc_totss) {
        par <- tt$data$parent[i]
        parX <- tt$data$x[par]
        tt$data[tt$data$node[i], "x"] = tt$data[tt$data$node[i], "x"]*15
    }

    off_base <- 0.6
    param <-0
    for (j in 1:length(heat_map_title_list)) {
        print(j)
        print("off_base is")
        print(off_base)
        print("param is")
        print(param)
        print("man_adj_heat_loc is")
        print(man_adj_heat_loc)

        off <- off_base +param+man_adj_heat_loc
        if (j>1) {
            off <- off+man_adj_heat_loc2
        }
        if (j>2) {
            off <- off+man_adj_heat_loc3
        }

        if (j>1) {
            tt <- pr440_short_tips_TRY_heat + new_scale_fill()
        }
        
        print("off is")
        print(off)
        print("length(colnames(dxdf440_for_heat[[j]])) is")
        print(length(colnames(dxdf440_for_heat[[j]])))
        wi <- 0.8*length(colnames(dxdf440_for_heat[[j]]))/23
        print("wi is")
        print(wi)


            pr440_short_tips_TRY_heat <- gheatmap(tt ,data= dxdf440_for_heat[[j]], colnames_angle=45, offset = off, width=wi,
                                          font.size = size_tip_text*1.6, colnames_offset_x=0.1, colnames_offset_y=6.8,
                                          legend_title = heat_map_title_list[[j]]) +
        #scale_fill_viridis_c(option = "plasma", name=heat_map_title_list[[j]], direction=-1)
        scale_fill_gradient2(
    low = "beige", 
    mid = "seashell2", 
    high = "firebrick4", 
    midpoint = .02,
    name=heat_map_title_list[[j]]       
  )
      
        param <- param+wi*length(colnames(dxdf440_for_heat[[j]]))/4.5

        if(j==1) {
            pr440_short_tips_TRY_heat <- pr440_short_tips_TRY_heat+
                    theme(legend.text = element_text(size = text_legend_size*text_mu,
                                                     margin = margin(t = mar/2, b=mar/2, unit = "pt")),
          legend.title = element_text(size = text_legend_title_size,face = "bold"),
         legend.spacing.y = unit(text_mu/12, 'cm'),
          legend.spacing.x = unit(0.7, 'cm'),
             legend.spacing=unit(0.25, 'cm'),
             legend.key.width=unit(1, 'cm'),
             legend.key.height=unit(0.5,"line"))
            
        off_base <- off
        }
        else {
            if (j==2) {
                pr440_short_tips_TRY_heat <- pr440_short_tips_TRY_heat+
                scale_fill_viridis_d(option="B", name=heat_map_title_list[[j]], direction=-1, na.value ="gray")+
                scale_fill_hue(name=heat_map_title_list[[j]])
            } else {
                pr440_short_tips_TRY_heat <- pr440_short_tips_TRY_heat+
                scale_fill_viridis_d(option="A", name=heat_map_title_list[[j]], direction=1, na.value ="gray")+
                scale_fill_hue(h = c(100, 300),name=heat_map_title_list[[j]]) 
                

            }
            
        }

    }
    p <- pr440_short_tips_TRY_heat  #tree_TRY #pr440_short_tips_TRY_heat #temp
       }
    ##################################
    #######################
    p <- p
 
    
    ####

     ggsave(out_file_path, width = width, height = height, units = units_out, limitsize = FALSE)

    out <- c()
    out[['p']] <- p
    out <- out
    
}

func.make.plot.tree.NEW <- function( tree440,
                    dx_rx_types1_short,
                    list_id_by_class,
                    dxdf440_dataf,
                    title.id,
                    FDR_perc,
                    no_name,
                    rotate_flag,
                    rotation_params1,
                    rotation_params2,
                    flag_short_tips , 
                    tips_length ,
                    show_boot_flag ,
                    FLAG_BULK_DISPLAY ,
                    how_many_hi ,
                    high_label_list ,
                    high_color_list ,
                    high_title_list ,
                    lists_list_hi ,
                    simulate.p.value,
                    width,
                    height,
                    colors_scale1 ,
                    out_file_path ,
                    edge_width_multiplier=1 ,
                    size_tip_text=3,
                    size_font_legend_title=30,
                    size_font_legend_text=20 ,
                    size_font_legend_box=15,
                    labels_not_in_legend,
                    no_name_color,
                    debug_mode,
                    man_adjust_elipse=0,
                    man_multiply_elipse=1,
                    man_adj_second_legend= 0.15,
                    man_space_second_legend= -0.02,
                    laderize_flag=FALSE, 
                    cls_renaming_list ,
                    title_flag ,
                    title_replace ,
                    flag_classification_format,
                    id_tip_trim_flag, 
                    id_tip_trim_start,
                    id_tip_trim_end,
                    units_out) {
    
    
    
    if (debug_mode==TRUE){
        print("In func.make.plot.tree")
    }
    
    dx_rx_types1 <- dx_rx_types1_short
    
    pr440 <- ggtree(tree440) 
    d440 <- pr440$data
    cc_tipss <- func.create.cc_tipss(d440)
    cc_nodss <- func.create.cc_nodss(d440)
    cc_totss <- func.create.cc_totss(d440)
    
    nods_num <- length(cc_nodss)
    tips_num <- length(cc_tipss)
    tree_size <- nods_num +tips_num 
    
    if (debug_mode==TRUE){
        print("tree data is")
        print(d440)
    }
    

    
    list_node_by_class <- func.create.list_node_by_class(dx_rx_types1_short,
                                                     list_id_by_class,dxdf440_dataf,
                                                         title.id,tree_size,d440,cc_totss,debug_mode,
                                                         id_tip_trim_flag,id_tip_prefix)
    

    list_rename_by_class <-list_node_by_class
    
    
    y_off_base= -8
    

    
    tree_with_group <- groupOTU(tree440, list_node_by_class)

    

        
    subframe_of_nodes <- d440[d440$isTip=="FALSE",]
    cc_nodss90 <- func.create.cc_nodss90(subframe_of_nodes)
    cc_nodss80 <- func.create.cc_nodss80(subframe_of_nodes)
    cc_nodss70 <- func.create.cc_nodss70(subframe_of_nodes)
    
    
    #create array of number of nodes of each type for the whole tree
    df_count_FULL_tree_populations <- data.frame(idx = 1:length(dx_rx_types1_short),
                type = dx_rx_types1_short,
                count = rep(0, length(dx_rx_types1_short))
                )
    
    for (opt in dx_rx_types1){
        indx <- which(dx_rx_types1_short == opt)
        df_count_FULL_tree_populations[indx,'count'] <- length(list_id_by_class[[opt]])-1
    }
 
    
    tree_with_group_CPY <- groupOTU(tree440, list_rename_by_class)
    

    
    levels_groups <- levels(ggtree(tree_with_group_CPY)$data$group)
    
    
    
   
    
    tree_TRY <- ggtree(tree_with_group_CPY,aes(color=new_class, size=p_val_new),ladderize=laderize_flag,
                  ) #size=p_val
    
     new_colors_list<- func.create.new_colors_list(FDR_perc,tree_TRY,tree_with_group,no_name,tree_size)
    
    if (no_name %in% new_colors_list) {
        flag_no_name <- 1
    } else {
        flag_no_name <- 0
    }
    
    tree_TRY$data$new_class <- new_colors_list

    
    op_list <- c(paste0("p>",FDR_perc ), 
             paste0(paste0(0.5*FDR_perc,"<p<="),FDR_perc),
             paste0("p<=",0.5*FDR_perc ))
    
    p_list_of_pairs<- func.create.p_list_of_pairs(list_node_by_class,d440,dx_rx_types1_short,
                                            cc_nodss,tree_with_group,FDR_perc,tree, cc_tipss,
                                            tree_TRY,tree_size,no_name,simulate.p.value)

    
    p_PAIRS_pval_list <- func.create.p_val_list_FROM_LIST(FDR_perc,tree_TRY,p_list_of_pairs,op_list)
    
   

    tree_TRY$data$p_val_new <-  p_PAIRS_pval_list 

    tree_TRY$data$p_val_new <- factor(tree_TRY$data$p_val_new,                 # Relevel group factor
                         levels = op_list)

    
    if (FLAG_BULK_DISPLAY== TRUE){
           
        
        for (index_high in 1:how_many_hi) {
            list_high_for <- lists_list_hi[[index_high]]
            #nam_high <- paste('high',as.character(index_high))
            if (index_high==1) {
                tree_TRY$data$'high1' <-  list_high_for 
                #tree_TRY$data$high1 <- factor(tree_TRY$data$high1,                 # Relevel group factor
                 #        levels = c(TRUE,FALSE))
            } else {
                if (index_high==2) {
                   tree_TRY$data$'high2' <-  list_high_for 
                    #tree_TRY$data$high2 <- factor(tree_TRY$data$high2,                 # Relevel group factor
                    #     levels = c(TRUE,FALSE)) 
                } else {
                    if (index_high==3) {
                   tree_TRY$data$'high3' <-  list_high_for 
                #    tree_TRY$data$high3 <- factor(tree_TRY$data$high3,                 # Relevel group factor
                #         levels = c(TRUE,FALSE)) 
                    } else {
                        stop ("Too many highlight options. Only 3 are supported")
                    }
                }
            
            
            }
   
        }

    #tree_TRY$data$high_id <-  list_high_with_id 
    #tree_TRY$data$high_id <- factor(tree_TRY$data$high_id,                 # Relevel group factor
    #                     levels = unique(list_high_with_id))
    }
    
    tree_TRY1<- tree_TRY
    
        
  
    
    
       #rotate tree if required
    
    if (rotate_flag %in% c("RX_first", "FRAC_first") ) {
 
        list_weight_dx.rx <- rotation_params1[['list_weight_dx.rx']]
        list_weight_frac <- rotation_params2[['list_weight_dx.rx']]
        TREE_OTU_dx.rx  <- rotation_params1[['TREE_OTU_dx.rx']]
        TREE_OTU_frac <- rotation_params2[['TREE_OTU_dx.rx']]
    }
    
    if (rotate_flag =="RX_first"){
        

        

        tree_TRY2 <- func.rotate.tree.based.on.weights (tree_TRY1,list_weight_dx.rx,list_weight_frac,
                                    TREE_OTU_dx.rx,TREE_OTU_frac,tree_size)
 
        
        } else if (rotate_flag =="FRAC_first"){
        
            tree_TRY2 <- func.rotate.tree.based.on.weights (tree_TRY1,list_weight_frac,list_weight_dx.rx,
                                    TREE_OTU_frac,TREE_OTU_dx.rx,tree_size)
        } else {
        
            tree_TRY2 <- tree_TRY1
        }
    
       #tree with short tips if required
        pr440_short_tips_TRY <- tree_TRY2 #tree_TRY
    
        if (flag_short_tips == TRUE) {
        for(i in cc_tipss) {
            par <- pr440_short_tips_TRY$data$parent[i]
            parX <- pr440_short_tips_TRY$data$x[par]
            pr440_short_tips_TRY$data[pr440_short_tips_TRY$data$node[i], "x"] = parX +tips_length
        }
    
        if (debug_mode== TRUE) {
            print("##DEBUG##")
            print("flag_short_tips is")
            print(flag_short_tips)
        
        }
    } 
 
    
    
        #make column for bootstrap values only in the tree, used in display
    zs =rep(0, tree_size)
    pr440_short_tips_TRY$data$boot_val <-  zs
    for (i in 1:tree_size) {
 
        if (i %in% cc_tipss) {
            pr440_short_tips_TRY$data$boot_val[i] = 0
        } else {
            pr440_short_tips_TRY$data$boot_val[i] =pr440_short_tips_TRY$data$label[i]
        }
    }
    
    cls_renaming_list_with_no_name <- c(cls_renaming_list,no_name)
    
    
    if (flag_no_name== 1) {

        fin_color_list <- cls_renaming_list_with_no_name
    } else {

        fin_color_list <- cls_renaming_list
    }
    

    

    
     pr440_short_tips_TRY$data$new_class <- factor(pr440_short_tips_TRY$data$new_class,                 # Relevel group factor
                         levels = unique(dx_rx_types1_short))
    

    
    
    list_of_sizes <- c(1,  1.6,  2.2)*edge_width_multiplier
    

    pr440_short_tips_TRY_new <- pr440_short_tips_TRY + scale_color_manual(values=colors_scale1) + 
    scale_size_manual(values=list_of_sizes)
    
    
    
        
    t<- pr440_short_tips_TRY_new

    #Add bootstrap values to the tree if required

   
    x_range_min = min(pr440_short_tips_TRY_new$data$x)
    x_range_max = max(pr440_short_tips_TRY_new$data$x)
    
    print("x_range_min is")
    print(x_range_min)
    print("x_range_max is")
    print(x_range_max)
    
    square_size_legend <- size_font_legend_box # 30  #16
    text_legend_title_size <- size_font_legend_title # 50 #30
    text_legend_size <- size_font_legend_text #35 # 25
    
    size_90 <- 26
    size_80 <- 24
    size_70 <- 22
    
    if (show_boot_flag == TRUE) {

    pr440_short_tips_TRY_new_with_boot <- pr440_short_tips_TRY_new + 
    geom_nodepoint(position=position_nudge(x = 0.005, y = 0),aes(subset = boot_val >= 0.9 ),size = size_90,
               shape=24, fill="grey36",colour = "grey20",show.legend = FALSE,alpha=1/2) + 
    geom_nodepoint(position=position_nudge(x = 0.005, y = 0),aes(subset = boot_val >= 0.8 & boot_val < 0.9 ),size = size_80,
               shape=24, fill="grey36",colour = "grey20",show.legend = FALSE,alpha=1/2) + #orchid3
    geom_nodepoint(position=position_nudge(x = 0.005, y = 0),aes(subset = boot_val >= 0.7 & boot_val < 0.8),size = size_70,
               shape=24, fill="grey36", colour = "grey20",show.legend = FALSE,alpha=1/2)  #orchid1 #"#b5e521"
        
    } else {

        pr440_short_tips_TRY_new_with_boot <- pr440_short_tips_TRY_new
    }
    
    
    
        
    # add legends to the tree
    group_display_title <- "Cell type"
    if (exists('title_flag')) {
        if (title_flag== TRUE) {
            group_display_title <- title_replace
        }
    } 
    

    

    new_base <-0
    
    pr440_short_tips_TRY_new_with_boot<- pr440_short_tips_TRY_new_with_boot#+geom_treescale(x=-20,y=0) 

    
    pr440_short_tips_TRY_new_with_boot_more1 <- pr440_short_tips_TRY_new_with_boot
    

    pr440_short_tips_TRY_new_with_boot_more1 <- pr440_short_tips_TRY_new_with_boot_more1 + 
    layout_dendrogram() +
    guides(colour = guide_legend(title = group_display_title, 
                                 override.aes = list(size = square_size_legend*1.5, alpha = 1,linetype=12, linewidth=30),
                                breaks = levels(pr440_short_tips_TRY_new_with_boot_more1$data$group)),
       size = guide_legend(title = "p value"),
       shape = guide_legend()) + 
    theme(legend.text = element_text(size = text_legend_size*3.1,margin = margin(t = 40, b=40, unit = "pt")),
      legend.title = element_text(size = text_legend_title_size*3.1,face = "bold"),
     legend.spacing.y = unit(3, 'cm')) +
    scale_y_reverse()  +
    geom_rootedge() 
    
    
   
    

    
    x_range_min <-  min(pr440_short_tips_TRY_new_with_boot_more1$data$x)
    x_range_max <-  max(pr440_short_tips_TRY_new_with_boot_more1$data$x)
    x_off_base<-  x_range_min*0.6  -man_adj_second_legend -3 #-0.9 #-0.76
    print("x_off_base is")
    print(x_off_base)

    # highlight bullk if required
    
    
        
    pr440_short_tips_TRY_new_with_boot_more1 <- pr440_short_tips_TRY_new_with_boot_more1 +
    geom_tiplab(size =size_tip_text,angle=90,hjust=1,,show.legend = FALSE) #+
    #guides(color = guide_legend(override.aes = list(label = "U", size = 7)))
     #guides(color = guide_legend(override.aes = list(label = "\u25A0", size = 7)))
    
        if (FLAG_BULK_DISPLAY ==TRUE) {
        
    b <- 0.28 #+ ((100-width)*0.005)#0.45 #0.45 + (height-60)*0.04
    a <-   0.017 + ((100-height)*0.0004) #0.035 + (height-60)*(-0.002) #0.15
    x_adj_hi <- 0.017  + ((100-height)*0.0001)#0.025

    if (x_adj_hi<=0) {
        x_adj_hi <-0.01
    }
    if (a<=0) {
        a<-0.01
    }
    if (b<=0) {
        b<-0.01
    }
    
    b<- 0.45 #0.15 #0.45 /1.5 * (x_range_min *(-1))

    a<- 0.035 /1.5 * (x_range_min *(-1)) *man_multiply_elipse #0.01 #0.035
     x_adj_hi <- 0 #0.017
        
    up_offset<- -3
    a<- 0.02
    
    y_off_base<-  -8
    x_off_base<- x_range_min +0.1 #x_range_min*0.6 -man_adj_second_legend  #-0.9 #-0.76
            

            
            
        for (index_high in 1:how_many_hi) {
            if (index_high==1) {
                high_nodes_table1 <- pr440_short_tips_TRY$data[tree_TRY$data$high1 == TRUE,]
   
                pr440_short_tips_TRY_new_with_boot_more1 <- pr440_short_tips_TRY_new_with_boot_more1 +
            geom_ellipse(data=high_nodes_table1,
                         aes(x0 = (max(pr440_short_tips_TRY$data[,'x'])-x)*(-1)+x_adj_hi +man_adjust_elipse+0.01, 
                                           y0 = y, a = a, b = b, angle = 0),
                    fill=high_color_list[[1]],alpha=0.5,linetype="blank",show.legend=FALSE)
            } else {
                if (index_high==2) {
                    high_nodes_table2 <- pr440_short_tips_TRY$data[tree_TRY$data$high2 == TRUE,]
                    pr440_short_tips_TRY_new_with_boot_more1 <- pr440_short_tips_TRY_new_with_boot_more1 +
                    geom_ellipse(data=high_nodes_table2,
                         aes(x0 = (max(pr440_short_tips_TRY$data[,'x'])-x)*(-1)+x_adj_hi +man_adjust_elipse+0.005, 
                                           y0 = y, a = a, b = b, angle = 0),
                        fill=high_color_list[[2]],alpha=0.5,linetype="blank",show.legend=FALSE)
                } else {
                    if (index_high==3) {
                        high_nodes_table3 <- pr440_short_tips_TRY$data[tree_TRY$data$high3 == TRUE,]
                        pr440_short_tips_TRY_new_with_boot_more1 <- pr440_short_tips_TRY_new_with_boot_more1 +
                            geom_ellipse(data=high_nodes_table3,
                             aes(x0 = (max(pr440_short_tips_TRY$data[,'x'])-x)*(-1)+x_adj_hi +man_adjust_elipse+0.005, 
                                           y0 = y, a = a, b = b, angle = 0),
                            fill=high_color_list[[3]],alpha=0.5,linetype="blank",show.legend=FALSE)
                    }
                }
            }
            
        }
                 

    
            
    p<- pr440_short_tips_TRY_new_with_boot_more1
    p1 <- p$layers[1]
    p2 <- p$layers[2]
    p7 <- p$layers[7]
    p$layers[2] <- p7 
    p$layers[7] <- p2
        

            
    stair <- 0.085
    multiple_high_up_offset <- stair*(how_many_hi-1)*(-1)
            
        for (index_high in 1:how_many_hi) {
            
            multiple_high_down_offset <- (index_high-1)*stair
            if (index_high==1) {
                p <- p +
            annotate(
            geom = "text", 
            label = high_title_list[[index_high]],  
                size = round(size_font_legend_title*man_multiply_second_legend), 
                x = x_off_base+multiple_high_up_offset+multiple_high_down_offset,
                y = y_off_base,  fontface = "bold"
              )
            }
            p <- p +
         annotate(
        geom = "text", 
        label = high_label_list[[index_high]],  
             size = round(size_font_legend_text*man_multiply_second_legend_text), 
             x = x_off_base+0.035+multiple_high_up_offset+multiple_high_down_offset,
             y = y_off_base-2,  
      )+
        geom_ellipse(aes(x0 = x_off_base+0.035+multiple_high_up_offset+multiple_high_down_offset ,
                         y0 = y_off_base+1.7, a = a, b = b, angle = 0),
            fill=high_color_list[[index_high]],alpha=0.5,linetype="blank",show.legend=FALSE)
            
          new_base <- 0.035+multiple_high_up_offset+multiple_high_down_offset  
            
            
        }
    } else {
 
        pr440_short_tips_TRY_new_with_boot_more1 <-pr440_short_tips_TRY_new_with_boot_more1 +
        geom_tiplab(size =size_tip_text,angle=90,hjust=1,,show.legend = FALSE)
        p <- pr440_short_tips_TRY_new_with_boot_more1
    }
    

    
        #add second legend for bootstrap and highlight if required
    
    second_legend_title_font_size<-  round(text_legend_title_size/3,digits = 0)  #16 # 10.5
    second_text_legend_size<- round(text_legend_size/3,digits = 0) #8.6

    boot_spacing_legend<- 0.04+   man_space_second_legend
    
    
       boot_spacing_legend<- 0.04+    man_space_second_legend
    

    #second_legend_title_font_size
    
        if (show_boot_flag == TRUE) {
    p<- p+
     annotate(
    geom = "text", 
    label = "Bootstrap",  size = size_font_legend_title, x = x_off_base+ new_base +0.07, y = y_off_base,  fontface = "bold"
      )+ 
 annotate(
    geom = "text", 
    label = ">70%",  size = size_font_legend_text, x = x_off_base+ new_base +0.08+ boot_spacing_legend, y = y_off_base-2,  
  )+
 annotate(
    geom = "text", 
    label = ">80%",  size = size_font_legend_text, x =  x_off_base+ new_base +0.08+ boot_spacing_legend*2, y = y_off_base-2,  
  )+
 annotate(
    geom = "text", 
    label = ">90%",  size = size_font_legend_text, x = x_off_base+ new_base +0.08+ boot_spacing_legend*3, y = y_off_base-2,  
  ) +
 annotate(
    geom = "point", 
    shape=24,  size = size_90, x = x_off_base+ new_base +0.08+ boot_spacing_legend, y = y_off_base+1,  fill="grey36",colour = "grey20",alpha=1/2,
  )+
 annotate(
    geom = "point", 
    shape=24,  size = size_80, x = x_off_base+ new_base +0.08+ boot_spacing_legend*2, y = y_off_base+1,  fill="grey36",colour = "grey20",alpha=1/2,
  ) +
 annotate(
    geom = "point", 
    shape=24,  size = size_70, x = x_off_base+ new_base +0.08+ boot_spacing_legend*3, y = y_off_base+1,  fill="grey36",colour = "grey20",alpha=1/2,
  ) 
            
    }
 p<- p+ xlim(x_range_max+0.02,x_range_min-0.02)
    
    
     ggsave(out_file_path, width = width, height = height, units = units_out, limitsize = FALSE)
    
    out <- c()
    out[['p']] <- p
    out <- out
    
}

func.make.highlight.params.NEW <- function(yaml_file,title.id,ids_list,tree440,readfile440,hi_def,debug_mode) {
    
    

    len_hi <- length(hi_def$according)
    indexes_hi<- 1:len_hi
    
    how_many_hi <- length(indexes_hi)
    high_label_list <- c()
    for (in_hi in indexes_hi) {
        high_label_list[[in_hi]] <- hi_def$according[[in_hi]][[as.character(in_hi)]]$display_name
    }
    
    high_color_list <- c()
    for (in_hi in indexes_hi) {
        high_color_list[[in_hi]] <- hi_def$according[[in_hi]][[as.character(in_hi)]]$color
    }
    
    high_title_list <- c()
    for (in_hi in indexes_hi) {
        high_title_list[[in_hi]] <- hi_def$according[[in_hi]][[as.character(in_hi)]]$display_title
    }
    
    tab <- ggtree(tree440)$data
    

    

    
    lists_list_hi <- c()
    for (in_hi in indexes_hi) {
        temp_list <-c()
        node_num_list <- tab$node
        title_i<- paste0('title',as.character(1))
        value_i<- paste0('value',as.character(1))
        
      
        
        title_i_name <-  hi_def$according[[in_hi]][[as.character(in_hi)]][[title_i]]
        value_i_name <-  hi_def$according[[in_hi]][[as.character(in_hi)]][[value_i]]
        
  
        
        if (! title_i_name %in% colnames(readfile440)) {
            print("wrong highlight parameter")
            stop(paste(title_i_name, "is not a title in classification file"))
        }
        
        for (node_num in node_num_list) {
            
            row <- subset(tab, node==node_num)
      
            is_high <- FALSE
            if (row$isTip== TRUE) {
                id <- substr(row$label,3,7)
                
 
                row_file <- subset(readfile440, Sample.Reads.ID== id)
       
                if (nrow(row_file)==1) {
                    val <- row_file[[title_i_name]]
    
                    if (is.na(val)) {
                        
                    } else {
                        if (val== value_i_name) {
                            is_high<- TRUE
                        }
                    }

                }
            }
            temp_list <-  c(temp_list, is_high)
            
        }
        lists_list_hi[[in_hi]]<- temp_list
    }
    highlight.params.NEW <- c()
    highlight.params.NEW$how_many_hi <- how_many_hi
    highlight.params.NEW$high_label_list <- high_label_list
    highlight.params.NEW$high_color_list <- high_color_list
    highlight.params.NEW$high_title_list <- high_title_list
    highlight.params.NEW$lists_list_hi <- lists_list_hi
    
    highlight.params.NEW <- highlight.params.NEW
        
}

func.check.if.id.in.sub.class <- function(val, title_i_name, value_i_name) {


            is_in <- FALSE
    if (is.null(val)) {
        
    } else {
        
    
                  if (length(value_i_name)>1) {
                      
                    if (val %in% value_i_name) {
                        is_in <- TRUE
  
                    }
                    if (is.na(val)) {
                        if ("na" %in% value_i_name) {
                            is_in <- TRUE
                        }
                        
                        if ("NA" %in% value_i_name) {
                            is_in <- TRUE
                        }
                    }

                } else {
                    
                   
                    if (substr(value_i_name,1,1)=='(') {
                        

                        be <- strsplit(value_i_name, "-")

                        begin <- be[[1]]

   
                        end <- begin[2]
                        begin <- begin[1]
  
                        begin <-  as.numeric(substr(begin,2,nchar(begin)))
          
                        begin <- begin[1]


                        
                        end <- as.numeric(substr(end, 1, nchar(end)-1))


                        if (is.na(begin) | is.na(end)) {
                            print("Error")
                            print(paste0(value_i_name," is an elegal value inside brackets"))
                            print("Inside brackets numerical range should be placed")
                            print("If the brackets are a part of the desired value, put them inside a list using []")
                            stop(" ")
                        }
                        if (begin>= end) {
                            stop("elegal range")
                        }
 
                        if (is.na(val)) {
                           
                        } else {
                            if ((begin <= val) & (val <= end)) {
                                is_in <- TRUE
                            }
                        }

                        
                    } else {

                        
                        if (is.na(val)) {
                            
                            if (tolower(value_i_name)== "na") {
                                
                                is_in <- TRUE
                            } else {
                                
                                is_in <- FALSE
                            }
                        } else {
                            
                            if (val == value_i_name) {
                                is_in <- TRUE
                            }
                        }

                    }
                }
    }
    is_in <- is_in
    
}

func.make.list_id_by_class <- function(cls_num, cls_renaming_list,yaml_file,
                                       title.id,leaves_id_from_tree,readfile440,acc,debug_mode,
                                       id_tip_trim_flag= TRUE,id_tip_trim_start=3,id_tip_trim_end=7) {
    list_id_by_class <- c()
    for (cls in cls_renaming_list) {

        list_id_by_class[[cls]] <- ""
    }
    



    
    list_id_is_in_clss_flags <- c()
    for (id_tip in leaves_id_from_tree) {
        list_id_is_in_clss_flags[[id_tip]] <- 0
    }
    


    for (id_tip in leaves_id_from_tree) {
        

        
        if (id_tip_trim_flag == TRUE) {
            row <- readfile440[readfile440[[title.id]] == substr(id_tip,id_tip_trim_start,id_tip_trim_end),]
            
        } else {
            row <- readfile440[readfile440[[title.id]] == id_tip,]
            
        }
        
        flag_is_in <-TRUE
        if (nrow(row)==0) {
            flag_is_in <- FALSE
        }
 
        
        


        for (ac in acc) {

            flag_op <-TRUE

            is_in_class <- c()
            index_op <- 1
            while(flag_op== TRUE) {

                
                title_i <- paste0('title',index_op)
                value_i <- paste0('value',index_op)
                
                index_op <- index_op+1

                names_ac <- names(ac)


            
                title_i_name <- ac[[names_ac]][[title_i]]
                value_i_name <- ac[[names_ac]][[value_i]]
                
  
                
                if ( (paste0('title',index_op)) %in% names(ac[[names_ac]])) {
                    
                 } else {
                    
                    flag_op <-FALSE
                }
  
                if (flag_is_in== TRUE) {
                    val <- row[[title_i_name]]
                    is_in <- func.check.if.id.in.sub.class(val, title_i_name, value_i_name)
                } else {
                    val <- NA
                    is_in <- FALSE
                }
                

     
                



                is_in_class <- c(is_in_class,is_in)


            }
 
            
                            
            
            if (sum(is_in_class)== length(is_in_class)) {
                is_in_class_res<- TRUE
            } else {
                is_in_class_res<- FALSE
            }
            
            
                          
            if (is_in_class_res == TRUE) {
                

                if (list_id_is_in_clss_flags[[id_tip]]==0) {

                    cls <- ac[[names_ac]][['display_name']]
                    

                    list_id_is_in_clss_flags[[id_tip]] <-1
           

                    list_id_by_class[[cls]] <- c(list_id_by_class[[cls]], id_tip)

                } else {
                    print(paste0("Warining: double mapping for ", id_tip))
                }
            }
        }
    }
    if (length(list_id_is_in_clss_flags[list_id_is_in_clss_flags==0])>0) {
        print("Warning: the following cells are not classified")
        print(list_id_is_in_clss_flags[list_id_is_in_clss_flags==0])
    }

    

    list_id_by_class<- list_id_by_class
}

func.calc.rotation.definitions <- function(rot1,rot2,yaml_file,title.id,ids_list,tree440,readfile440, debug_mode) {
    rotate_flag_for_title <- "no"
    rotate_str <- ""
    rotation_params1 <-c()
    rotation_params2 <-c()
    rotate1_types_list_dx.rx <-NA
    list_weight_frac1 <- NA
    
    if (func.check.bin.val.from.conf(rot1)== FALSE) {
        rotate_flag <-"no"
        rotate_flag_for_title <- "no"
        rotate_str <- "no_rotate"
    } else  {
        
        rotate_flag <- 'RX_first'
        rotate_flag_for_title <- "yes"

       
        rotation_params1 <- func.make.rot.params(1,yaml_file,title.id,ids_list,tree440,readfile440) 

        rotate1_types_list_dx.rx <- rotation_params1['types_list_dx.rx']
        
        if (debug_mode== TRUE) {
            print("##DEBUG##")
            print("rotation_params1 is")
            print(rotation_params1)
        }


        st <- c()
        for (ind_s in rotate1_types_list_dx.rx) {
            st <- paste0(st,ind_s, collapse = "_")
        }

        
        if (nchar(st)>30) {
            st <- substr(st,1,30)
            st <- paste0(st,"_ETC")
        }
        

        rotate_str <- paste0('rotate_by_',st)
        

        
        if (func.check.bin.val.from.conf(rot2)== FALSE) {
            rotation_params2 <- rotation_params1

        
        } else {
            rotation_params2 <- func.make.rot.params(2,yaml_file,title.id,ids_list,tree440,readfile440)  
            list_weight_frac1 <- rotation_params2['list_weight_dx.rx']
             
            if (debug_mode== TRUE) {
                print("##DEBUG##")
                print("rotation_params2 is")
                print(rotation_params2)
                print("list_weight_frac1=")
                print(list_weight_frac1)
            }
            
        }
    }
    
    l <- c()
    l$rotate_flag_for_title <- rotate_flag_for_title
    l$rotate_str <- rotate_str
    l$rotation_params1 <- rotation_params1
    l$rotation_params2 <- rotation_params2
    l$rotate_flag <- rotate_flag
    l$rotate1_types_list_dx.rx <- rotate1_types_list_dx.rx
    l$list_weight_frac1 <- list_weight_frac1
    
    l<- l
}

#read yaml file
func.read_yaml <- function(conf_yaml_path){
    check <- file.exists(conf_yaml_path)
    if (check == FALSE) {
        print(conf_yaml_path)
        stop("Yaml file does not exist")
    }
    siz <- file.size(conf_yaml_path)
    if (siz==0L) {
        print(yaml_file)
        stop("Empty yaml file")
    }
    
    con <- file(conf_yaml_path, "r")
    yaml_file<- read_yaml(con)
    close(con)
    yaml_file <- yaml_file
}

fix.readfile440.with.missing.leaves<- function(readfile440,title.id,tree440,ids_list,no_name,
                                               id_tip_trim_flag,id_tip_trim_start,id_tip_trim_end) {
   

    readfile440n <- readfile440
    tree_data <- ggtree(tree440)$data
    leaves_id_from_tree1 <- tree_data[tree_data$isTip==TRUE,'label']
    leaves_id_from_tree <- as.list(leaves_id_from_tree1['label'])$label
    


    for (id in leaves_id_from_tree) {
  

        if (id == "root" ) {
            idd <- "root"
        } else {
            if (id_tip_trim_flag== TRUE) {
                idd <- substring(id,id_tip_trim_start,nchar(id))
            } else {
                idd <- id
            }
            
        }

        

        if (!(idd %in%  ids_list)) {

            
            place <- nrow(readfile440n) + 1
            readfile440n[place,] = rep(NA,length(colnames(readfile440n)))
            readfile440n[place,title.id] <- idd
 

        }
    } 
    

    for (i in ids_list) {
        if (i=="root"){
            next
        }
        if (id_tip_trim_flag== TRUE) {
            temp <- paste0("ID",i)
        } else {
            temp <- i
        }
        
        if (!(temp %in% leaves_id_from_tree)) {

            vec <- readfile440n[[title.id]]

            hw <- which (vec==i)

            readfile440n <- readfile440n[-c(hw), ]
        }
    }
    

    readfile440_new <- readfile440n
}

func.check.bin.val.from.conf <- function(val) {
    #function for rotation flag
    if (tolower(val) == "yes"  || val == TRUE || tolower(val) == "true") {
        out <- TRUE
    } else {
        out <- FALSE
    }

    out <- out
}

func.make.rot.params <- function(rot_index,yaml_file,title.id,ids_list,tree440,readfile) {
    #rotation function

    rotation_name <- paste0('rotation',rot_index)

    if (rot_index==1){
        rot_class_list <- yaml_file[['visual definitions']]$rotation1$'according'
    } else {
        rot_class_list <- yaml_file[['visual definitions']]$rotation2$'according'
    }
    how_many_rot <- length(rot_class_list)
    column_base<-3
    
    


    list_id_by_rot<- func.make.mapping.list(how_many_rot,rot_class_list,yaml_file,title.id,ids_list,readfile)


    
    
    

    TREE_OTU_dx.rx <- groupOTU(tree440, list_id_by_rot)

    types_list_dx.rx <- names(list_id_by_rot)
    list_weight_dx.rx <- func.set.list_weight(types_list_dx.rx)
    
    ret_rot= c()
    
    ret_rot[['list_id_by_rot']] <- list_id_by_rot
    ret_rot[['TREE_OTU_dx.rx']] <- TREE_OTU_dx.rx
    ret_rot[['types_list_dx.rx']] <- types_list_dx.rx
    ret_rot[['list_weight_dx.rx']] <- list_weight_dx.rx
    
    ret_rot <- ret_rot

}

function.create.cls.list <- function(ids_list,dx_rx_types1_short,list_id_by_class,na_name) {
    # make list of classes

    cls <- c()
    for (j in ids_list) {
 
        l<-""
        for (n in dx_rx_types1_short) {
            b<- as.character(j)

            if (b %in% list_id_by_class[[n]]) {
                l<- n
            }
    
        } 

        if (nchar(l)==0)
            {
            
            l<- na_name
        }
        cls <- c(cls,l)
    }

    

    cls <- cls
}

#make list of tips
func.create.cc_tipss <- function(d440){
    subframe_of_tips<- d440[d440$isTip=="TRUE",2]
    list_tips <- list(subframe_of_tips$node)
    cc_tips=c(list_tips)
    cc_tipss <-cc_tips[[1]]
}

#make list of intermediate nodes
func.create.cc_nodss <- function(d440){
    subframe_of_nodes<- d440[d440$isTip=="FALSE",2]
    list_nodes <- list(subframe_of_nodes$node)
    cc_nodes=c(list_nodes)
    cc_nodss <-cc_nodes[[1]]
}

#make lists of all nodes
func.create.cc_totss <- function(d440){
    list_all_tree <- list(d440$node)
    cc_tot=c(list_all_tree)
    cc_totss <-cc_tot[[1]]
}

func.create.list_node_by_class <- function(dx_rx_types1_short,list_id_by_class,
                                           dxdf440_dataf,title.id,tree_size,d440,cc_totss,debug_mode,
                                           id_tip_trim_flag,id_tip_prefix) { 
    
    classication<- "Mapping"
    if(debug_mode==TRUE) {
        print("In func.create.list_node_by_class")
    }
    

    
    #mapping nodes to classes
    convert_id_to_node <- func.create.convert_id_to_node(tree_size,d440,cc_totss)

    
    if(debug_mode==TRUE) {
        print("convert_id_to_node is")
        print(convert_id_to_node)
        print("dx_rx_types1_short is")
        print(dx_rx_types1_short)
    }

    list_node_by_class=c()
    for(i in 1:length(dx_rx_types1_short)) {

        name<-dx_rx_types1_short[i]

        
  
 
        subframe<-dxdf440_dataf[is.element(dxdf440_dataf[[classication]], dx_rx_types1_short[i]),]

        if (nrow(subframe)==0) {
            vec=NA
            vec2 =NA
        } else {
            if (id_tip_trim_flag == TRUE) {
                vec=c(paste0(id_tip_prefix,subframe[[title.id]]))
            } else {
                vec=c(paste0("",subframe[[title.id]]))
            }
            
            vec2=c()
            for (v in vec) {
 
                if (v=="IDroot") {
                    vec2 <- c(vec2,"root")
                } else {
                    vec2 <- c(vec2,v)
                }
            }
        }
        vec <- vec2

 
        
        vec1 <- vec

        
        if (is.na(vec[1])&& length(vec)==1){
            
            vec1 <- c(" ")
            list_node_by_class[[name]]<-as.numeric(vec1)
            
        } else {
        
        for( i in 1: length(vec)) {

     
            vec1[i]<-as.numeric(convert_id_to_node[vec[i]])
            }
        vec1 <- c(" ", vec1)

        list_node_by_class[[name]]<-as.numeric(sort(as.numeric(vec1)))
        }
    }
    if(debug_mode==TRUE) {
        print("list_node_by_class is")
        print(list_node_by_class)
    }
    
    list_node_by_class<- list_node_by_class
}

func.create.convert_id_to_node <- function(tree_size,d440,cc_totss){
    #converting id to node number
    convert_node_to_id <- vector(mode="list", length=tree_size)
    convert_id_to_node <- vector(mode="list", length=tree_size)
    names(convert_node_to_id) <- c(d440$node)
    names(convert_id_to_node) <- c(d440$label)
    for(i in cc_totss){

        convert_node_to_id[[i]] <-  d440$label[i]
        convert_id_to_node[[i]] <- d440$node[i]
    }
    convert_id_to_node <- convert_id_to_node
}

#list of nodes ids for nodes with bootstrap greater than 0.9
func.create.cc_nodss90 <- function(subframe_of_nodes){
    subframe_of_nodes90<- subframe_of_nodes[subframe_of_nodes$label >=0.9,]

    list_nodes90 <- list(subframe_of_nodes90$node)

    cc_nods90=c(list_nodes90)

    cc_nodss90 <-cc_nods90[[1]]

    temp <-c()
    for (el in cc_nodss90) {
        temp <- c(temp,el)
    }
    temp <- temp
}


#list of nodes ids for nodes with bootstrap greater than 0.8 and smaller than 0.9
func.create.cc_nodss80 <- function(subframe_of_nodes){
    subframe_of_nodes80<- subframe_of_nodes[subframe_of_nodes$label >=0.8 & subframe_of_nodes$label <0.9,]
    list_nodes80 <- list(subframe_of_nodes80$node)
    cc_nods80=c(list_nodes80)
    cc_nodss80 <-cc_nods80[[1]]
}


#list of nodes ids for nodes with bootstrap greater than 0.7 and smaller than 0.8
func.create.cc_nodss70 <- function(subframe_of_nodes){
    subframe_of_nodes70<- subframe_of_nodes[subframe_of_nodes$label >=0.7 & subframe_of_nodes$label <0.8,]
    list_nodes70 <- list(subframe_of_nodes70$node)
    cc_nods70=c(list_nodes70)
    cc_nodss70 <-cc_nods70[[1]]
}

func.create.new_colors_list<- function(FDR_perc,tree_TRY,tree_with_group,no_name,tree_size){
    #Define color to each node in the tree
    list_uncertain <- c()
    list_certain <- c() 
    list_equals <- c()
    new_colors_list <- rep(no_name,tree_size)
    for (nod in tree_TRY$data$node) {

        isT <- tree_TRY$data[[nod,"isTip"]]
        count_list <-c()
        if (isT == TRUE) {
            new_colors_list[nod]<- as.character.factor(tree_TRY$data[[nod,"group"]])
        } else {
            kids <- ggtree:::getSubtree(tree_with_group, nod)
            tips_list <-c()
            for (j in kids) {
                check <- tree_TRY$data[[j,"isTip"]]
                if (check == TRUE){
                    tips_list <- append(tips_list,j)
                }
            }
            for (j in tips_list) {
                color <- as.character.factor(tree_TRY$data[[j,"group"]])
 
                if (color %in% names(count_list)) {
                    count_list[color] <- count_list[color] +1
                } else {
                    count_list[color] <- 1
                }
                
            }

            identical_test <- length(which(count_list ==max(count_list)))
            if (identical_test ==1){ 
                color <- names(which(count_list ==max(count_list)))

                new_colors_list[nod]<- color

                list_certain <- append(list_certain, nod)
            } else {
                list_uncertain <- append(list_uncertain, nod)
                vec <- names(which(count_list ==max(count_list)))
                list_equals[[nod]] <- vec
            }
        }

    }
    for (nod in list_uncertain){

        flag_end <- FALSE
        flag <- FALSE
        prev_par <- nod
        color <- no_name
        color_temp <- no_name
        vec <- list_equals[[nod]]
        while (flag_end == FALSE & flag ==FALSE){
            par <- tree_TRY$data[[prev_par,"parent"]]
            if (par == prev_par){
                flag_end <- TRUE
            }
            if (par %in% list_certain){
                flag <- TRUE
                color_temp <- new_colors_list[par]
                if (color_temp %in% vec) {
                    color <- color_temp
                }
            } 
            prev_par <- par
        }
        new_colors_list[nod] <- color
        

    }
    new_colors_list <- new_colors_list
}


func.create.p_list_of_pairs <- function(list_node_by_class,d440,dx_rx_types1_short,
                                            cc_nodss,tree_with_group,FDR_perc,tree, cc_tipss,
                                            tree_TRY,tree_size,no_name,simulate.p.value) {
    
    in_sub_tree <- rep(0,2)
    Not_in_sub_tree <- rep(0,2)
    p_list_of_pairs<- rep(1,tree_size)
    for(nod in cc_nodss) {   
  
        kids <- ggtree:::getSubtree(tree_with_group, nod)  #kids contain all the nods in the subtree including nod (its root)
        sub_tree_size <- length(kids)
        tips_list <-c()
        
        in_sub_tree <- rep(0,2)
        Not_in_sub_tree <- rep(0,2)
  
        hyp <- tree_TRY$data[nod,"new_class"]
        if (hyp != no_name){
            for (j in kids) {
                check <- tree_TRY$data[[j,"isTip"]]
                if (check == TRUE){
                    tips_list <- append(tips_list,j)
                }
            }
        
            for (j in cc_tipss) {
                j_typ <- tree_TRY$data[j,"new_class"]
                if (j %in% tips_list){
                    if (j_typ == hyp) {
                        in_sub_tree[1] <- in_sub_tree[1]+1
                    } else {
                        in_sub_tree[2] <- in_sub_tree[2]+1
                    }
                } else {
                    if (j_typ == hyp) {
                        Not_in_sub_tree[1] <- Not_in_sub_tree[1]+1
                    } else {
                        Not_in_sub_tree[2] <- Not_in_sub_tree[2]+1
                    }
                }
     
            }

        Not_in_sub_tree1 <-  abs(Not_in_sub_tree)
        in_sub_tree1<- abs(in_sub_tree)
            
   
     
        xtab <- as.table(rbind(
          in_sub_tree1,
          Not_in_sub_tree1
        ))
 
        dimnames(xtab) <- list(
            In_sub_tree = c("Yes", "No"),
            Groups = c("In_type","Not_in_type")
        )
   
        conf_comp <- 1- FDR_perc
            

   

        test <- fisher.test(xtab,conf.level = conf_comp,simulate.p.value=simulate.p.value) #conf.level = 0.95
        p_val <- test$p.value
        p_list_of_pairs[nod]<- p_val  
        }
        
    }

    p_list_of_pairs <- p_list_of_pairs
}


func.create.p_val_list_FROM_LIST <-function(FDR_perc,tree_TRY,p_list_of_pairs,op_list){
    
    #Define p_val for each node in the tree

p_val_list <- c()



for (nod in tree_TRY$data$node) {
    p_val_op <- 0
    p_val_precise <- p_list_of_pairs[nod]
    
    if (FDR_perc < p_val_precise) {
        
        p_val_op <-op_list[1]
    } 
    if(FDR_perc >= p_val_precise  & p_val_precise >0.5*FDR_perc) {
        p_val_op <-op_list[2]
    }
    if(0.5*FDR_perc >= p_val_precise  ) {
        
        p_val_op <-op_list[3]
    }

    p_val_list <- append(p_val_list, p_val_op)
}
    p_val_list<- p_val_list
}

#mapping id to class
func.make.mapping.list<- function(cls_num,cls_list,yaml_file,title.id,ids_list,readfile,debug_mode =FALSE,no_nam){

    if(debug_mode==TRUE) {
        print("In func.make.mapping.list")
        #print("ids_list is is")
        #print(ids_list)
    }
    

    
    nams <- c()
    for (i in names(cls_list)) {
        value <- cls_list[[i]]$'value'
        nams <- c(nams,value)
    }
    


    class1 <- c()
    for (ty in 1:length(nams)){

        name <- nams[ty]


        vec_num =c("")
        vec = c("")
        for (ii in ids_list) {

            if (ii == "root") {
                i <- "root"
            } else {
                if (nchar(ii)>5) {
                    i <- substring(ii,3,nchar(ii))
                } else {
                    i <- ii
                }
                
            }
            

            check <- func.check.if.id.in.class.NEW(name,i,yaml_file,readfile,title.id,cls_list) 

            
            if (check == TRUE) {
                if (i=="root") {
                    vec <- c(vec, i)
                } else {
                    vec <- c(vec, paste0('ID',i))
                }
                
  
            }
        }
    
        

        
        class1[[name]]<-vec
 
    }
    

    list_id_by_class <- class1
}
    

func.check.if.id.in.class.NEW <- function(class_name,id,yaml_file,readfile,title.id,cls_list) {
    lin <- readfile[which(readfile[title.id] ==id ), ]

    
    ans <- FALSE
    for (i in names(cls_list)) {
        title <-cls_list[[i]]$'title'

        value <- cls_list[[i]]$'value'

        
        if (title %in% colnames(lin)) {
            val <- lin[[title]]
            if (length(val)>0){
                if (val== class_name) {
                    ans <- TRUE
                }
            }

        }
    }
    
    ans<- ans
}

#make list of weights for each nodes for rotation
func.set.list_weight <- function(types_list){
    list.weight <- c()
    start <-1
    jmp<- length(types_list)+1

    pow <- 1
    weight <- start
    for (type in types_list){
        list.weight[type]<- weight
        weight <- weight +((jmp^pow))
        pow <- pow+1
    }
    list.weight <- list.weight
}



func.rotate.tree.based.on.weights <- function(tree_TRY,list_weight_dx.rx,list_weight_frac,
                                              TREE_OTU_dx.rx,TREE_OTU_frac,tree_size){
    
#Rotating tree

    tree_return <- tree_TRY
    

    
    list_weights_for_nodes_dx.rx <- func.create.weight_list(tree_return,list_weight_dx.rx,TREE_OTU_dx.rx,tree_size)
    list_weights_for_nodes_frac <- func.create.weight_list(tree_return,list_weight_frac,TREE_OTU_frac,tree_size)
    

    
    if (NA %in% list_weights_for_nodes_dx.rx) {
        print("Error: rotation classes do not cover all tips")
    }
        
    tree_TRY$data$weight1 <- list_weights_for_nodes_dx.rx
    tree_TRY$data$weight2 <- list_weights_for_nodes_frac
    
    for (nod in tree_TRY$data$node){
        isT <- tree_TRY$data[[nod,"isTip"]]
        if (isT ==FALSE){
            children <- which(tree_TRY$data$parent==nod)
            children_weights <-rep(0,2)
            children_weights_SECOND <-rep(0,2)
            children_side_orig <- rep("l",2)
            children_side_new <- rep("l",2)
            if (length(children)==2){
                
                children_weights[1]<- list_weights_for_nodes_dx.rx[children[1]]
                children_weights[2]<- list_weights_for_nodes_dx.rx[children[2]]
                children_weights_SECOND[1]<- list_weights_for_nodes_frac[children[1]]
                children_weights_SECOND[2]<- list_weights_for_nodes_frac[children[2]]
                
                if (tree_TRY$data$y[children[1]]>tree_TRY$data$y[children[2]]){
                    children_side_orig[1]<- "l"
                    children_side_orig[2]<- "r"
                } else {
                    children_side_orig[1]<- "r"
                    children_side_orig[2]<- "l"
                }

                if (children_weights[1] < children_weights[2]) {
                    children_side_new[1] <- "l"
                    children_side_new[2] <- "r"    
                } else if (children_weights[1] >children_weights[2]) {
                    children_side_new[1] <- "r"
                    children_side_new[2] <- "l"   
                } else {
                    if (children_weights_SECOND[1] <children_weights_SECOND[2]){
                        children_side_new[1] <- "l"
                        children_side_new[2] <- "r"    
                    } else if (children_weights_SECOND[1] >children_weights_SECOND[2]){
                        children_side_new[1] <- "r"
                        children_side_new[2] <- "l"  
                    } else {
                        children_side_new[1] <-children_side_orig[1]
                        children_side_new[2] <-children_side_orig[2]
                    }
                }
                if (!children_side_new[1] == children_side_orig[1]){
                    tree_return<- flip(tree_return, children[1], children[2])
                }
            }
        } 
    }   
    tree_return <- tree_return
}

func.create.weight_list <- function(tree_TRY,weights_list,TREE_OTU_class,tree_size){
    #defining weitghs to tips for rotation
    

    
    g <- ggtree(TREE_OTU_class)
    

    
    weight_list <- rep(0,tree_size)
    for (nod in g$data$node){
        kids <- ggtree:::getSubtree(TREE_OTU_class, nod)
        tips_list <-c()
        weight <- 0
        s <- 0
        for (j in kids) {
            check <- tree_TRY$data[[j,"isTip"]]
            if (check == TRUE){
                tips_list <- append(tips_list,j)
            }
        }
        
        for (j in tips_list) {
            typ <- as.character.factor(g$data[[j,"group"]])
            
            weight_tip <- weights_list[typ]
            s<- s+ weight_tip
        }

        weight <- s/length(tips_list)

        weight_list[nod] <- weight
    }
    weight_list<- weight_list
}
