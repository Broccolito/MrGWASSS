merge_saige = function(autosomal_pattern = "FHS_EA_MRS_chrXXX.txt",
                       chrx_filename = "FHS_EA_MRS_plink2_run_chrX.txt",
                       autosomal_only = FALSE,
                       write_file = TRUE){

  files = as.list(
    sapply(1:22, function(x){
      gsub("XXX", x, autosomal_pattern)
    })
  )

  all_files = c(unlist(files), chrx_filename)
  sanity_check_table = tibble(
    File = all_files,
    Found = all_files %in% list.files()
  )
  print(kable(sanity_check_table))
  if(any(!sanity_check_table$Found)){
    return("Some of the GWAS results files specified do not exist...\n")
  }else{
    cat("\nAll GWAS results files are found...\n\n")
  }

  filter_autosomal_tables = function(f){
    d = as.data.frame(fread(f))
    d_filtered = d %>%
      filter(imputationInfo > 0.3) %>%
      filter(AC_Allele2 > 10) %>%
      filter(SE > 0) %>%
      filter(p.value >= 0, p.value <= 1) %>%
      filter(BETA != 0)
    cat(paste0(f, " Merged...\n"))
    return(d_filtered)
  }

  filter_chrx_table = function(chrx_filename, sample_size){
    d = as.data.frame(fread(chrx_filename))
    chrx_table = tibble(CHR = 23, POS = d$POS, SNPID = d$ID,
                        Allele1 = d$REF, Allele2 = d$ALT,
                        AF_Allele2 = d$A1_FREQ, N = sample_size,
                        BETA = d$BETA, SE = d$SE, Tstat = d$T_STAT, p.value = d$P) %>%
      mutate(AC_Allele2 = sample_size*AF_Allele2) %>%
      filter(AC_Allele2 > 10) %>%
      filter(SE > 0) %>%
      filter(p.value >= 0, p.value <= 1) %>%
      filter(BETA != 0)
    cat(paste0(chrx_filename, " Merged...\n"))
    return(chrx_table)
  }

  # Filter autosomal table
  autosomal_table = map(files, filter_autosomal_tables)
  autosomal_table = autosomal_table %>%
    reduce(full_join, by = names(autosomal_table[[1]]))

  if(!autosomal_only){
    # Filter sex chromosome table
    sample_size = max(autosomal_table$N)
    chrx_table = filter_chrx_table(chrx_filename, sample_size)
    chrx_table = full_join(head(autosomal_table,0), chrx_table,
                           by = names(chrx_table))

    # Merge the tables
    merged_table = rbind.data.frame(autosomal_table, chrx_table) %>%
      mutate(`1KG_ID` = SNPID) %>%
      mutate(NEA = Allele1) %>%
      mutate(EA = Allele2) %>%
      mutate(EAF = AF_Allele2) %>%
      mutate(P_gc = p.value) %>%
      mutate(SE_gc = SE) %>%
      select(`1KG_ID`, SNPID,
             CHR, POS, NEA, EA,AC_Allele2, EAF,
             imputationInfo, N, BETA, SE, Tstat,
             p.value, varT, varTstar, P_gc, SE_gc)
  }else{
    merged_table = autosomal_table %>%
      mutate(`1KG_ID` = SNPID) %>%
      mutate(NEA = Allele1) %>%
      mutate(EA = Allele2) %>%
      mutate(EAF = AF_Allele2) %>%
      mutate(P_gc = p.value) %>%
      mutate(SE_gc = SE) %>%
      select(`1KG_ID`, SNPID,
             CHR, POS, NEA, EA,AC_Allele2, EAF,
             imputationInfo, N, BETA, SE, Tstat,
             p.value, varT, varTstar, P_gc, SE_gc)
  }

  if(write_file){
    merged_filename = gsub("XXX", "_MERGED", autosomal_pattern)
    write.table(merged_table, file = merged_filename, sep = " ",
                quote = FALSE, row.names = FALSE)
  }

  return(merged_table)
}
