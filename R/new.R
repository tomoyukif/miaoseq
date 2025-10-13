
edit_result <- subset(edit_result,
                      select = -c(data_type, total_read),
                      subset = data_type == "genotype")

edit_result_pat <- apply(subset(edit_result,
                                select = -c(index_pair_id:name)), 1, paste, collapse = "_")
dup_list <- tapply(seq_along(edit_result_pat), sub("_.+", "", edit_result$name), function(i){
    dup <- !duplicated(edit_result_pat[i])
    out <- data.frame(i = i, dup = dup)
    return(out)
})
dup_list <- do.call("rbind", dup_list)
dup_list <- dup_list[order(dup_list$i), ]
edit_result$uniq <- "Dup"
edit_result$uniq[dup_list$dup] <- "Uniq"

long_edit_result <- pivot_longer(edit_result,
                                 cols = -c(index_pair_id:name, uniq),
                                 names_to = "gene",
                                 values_to = "edit")

long_edit_result$edit_eval <- sapply(long_edit_result$edit, function(x){
    x <- unlist(strsplit(x, "/"))
    if(all(is.na(x))){
        return(NA)
    }
    x1 <- gsub("[0-9]", "", x)
    x2 <- as.numeric(gsub("[a-zA-Z]", "", x))
    is_ref <- x1 %in% "ref"
    is_inframe <- x2 %% 3 %in% 0
    if(all(is_ref)){
        return("ref")

    } else if(all(!is_ref)){
        if(all(!is_inframe)){
            return("alt")

        } else if(all(is_inframe)){
            return("alt_inframe_homo")

        } else {
            return("alt_inframe_het")
        }

    } else {
        if(all(!is_inframe)){
            return("het")

        } else {
            return("het_inframe")
        }
    }
})

eval_levels <- c("ref",
                 "alt", "alt_inframe_het", "alt_inframe_homo",
                 "het", "het_inframe")
long_edit_result$edit_eval <- factor(long_edit_result$edit_eval, eval_levels)
long_edit_result <- subset(long_edit_result, name != "")
n_name <- length(unique(long_edit_result$name))

library(ggplot2)
for(i in unique(long_edit_result$plate)){
    p <- ggplot(subset(long_edit_result, plate == i)) +
        geom_tile(aes(x = sub("_.+", "", gene), y = 0, fill = edit_eval)) +
        geom_text(aes(x = (n_name + 1) / 2, y = 1.6, label = name), vjust = 1, hjust = 0.5, size = 3) +
        geom_text(aes(x = (n_name + 1) / 2, y = 1, label = uniq), vjust = 1, hjust = 0.5, size = 2.5) +
        facet_grid(rows = vars(well_row), cols = vars(well_col), switch = "y") +
        scale_fill_manual(values = c("yellow", "darkblue", "blue",
                                     "lightblue", "darkgreen", "green"),
                          breaks = eval_levels,
                          name = NULL) +
        labs(title = paste0("Plate ", i)) +
        theme(axis.title = element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7),
              axis.ticks.y = element_blank(),
              panel.grid = element_blank())
    pdf(file.path(wd, paste0("edit_viewer_plate", i, ".pdf")), width = 11, height = 6)
    print(p)
    dev.off()
}
