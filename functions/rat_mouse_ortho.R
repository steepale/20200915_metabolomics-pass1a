# Functions to relabel RNASeq read names to orthologs
###############################################################################################
# mmusculus_gene_ensembl: Mouse genes (GRCm38.p6)
# rnorvegicus_gene_ensembl: Rat genes (Rnor_6.0)
rat_mouse_ortho <- function(x, column = 'ENSEMBL_RAT', direction = 'rat2mouse') {
        # Ensure 'x' is a data.frame
        if ( class(x) != 'data.frame' ) {
                stop('x must be a data frame', class.= FALSE)
        }
        if ( column %!in% c('ENSEMBL_RAT','ENSEMBL_MOUSE') ){
                stop('column must be either ENSEMBL_RAT or ENSEMBL_MOUSE', class.= FALSE)
        }
        if ( direction %!in% c('rat2mouse', 'mouse2rat')){
                stop('direction must be either rat2mouse or mouse2rat')
        }
        
        # Load requirements
        library(biomaRt)
        library(purrr)
        library(dplyr)
        # Load in annotations
        mart_mm_ens = useMart('ensembl', dataset='mmusculus_gene_ensembl')
        mart_rn_ens = useMart('ensembl', dataset='rnorvegicus_gene_ensembl')
        # Create ortholog table
        if (direction == 'rat2mouse'){
                ortho_df <- getLDS(attributes=c('ensembl_gene_id',
                                                'mmusculus_homolog_orthology_confidence'),
                                   filters='ensembl_gene_id', 
                                   values = x[[column]], 
                                   mart=mart_rn_ens,
                                   attributesL=c('ensembl_gene_id'), 
                                   martL=mart_mm_ens) # Use biomart to get orthologs
                # Filter out any low confidence orthologs and any genes 
                #that are not one-to-one orthologs in both directions
                #ortho_df <- ortho_df[ortho_df$Mouse.orthology.confidence..0.low..1.high. == '1',]
                ortho_df <- ortho_df[!duplicated(ortho_df[,1]),]
                ortho_df <- ortho_df[!duplicated(ortho_df[,3]),]
                names(ortho_df) <- c('ENSEMBL_RAT','CONFIDENCE','ENSEMBL_MOUSE') 
                ortho_df <- ortho_df %>%
                        select(-CONFIDENCE)
                # Assign the symbols to a new column
                x <- left_join(x, ortho_df, by = 'ENSEMBL_RAT') %>%
                        filter(!is.na(ENSEMBL_RAT)) %>%
                        mutate(SYMBOL_MOUSE = mapIds(org.Mm.eg.db, ENSEMBL_MOUSE, 'SYMBOL', 'ENSEMBL'))
                
        }else{
                ortho_df <- getLDS(attributes=c('ensembl_gene_id',
                                                'rnorvegicus_homolog_orthology_confidence'),
                                   filters='ensembl_gene_id', 
                                   values = x[[column]], 
                                   mart=mart_mm_ens,
                                   attributesL=c('ensembl_gene_id'), 
                                   martL=mart_rn_ens) # Use biomart to get orthologs
                # Filter out any low confidence orthologs and any genes 
                #that are not one-to-one orthologs in both directions
                #ortho_df <- ortho_df[ortho_df$Rat.orthology.confidence..0.low..1.high. == '1',]
                ortho_df <- ortho_df[!duplicated(ortho_df[,1]),]
                ortho_df <- ortho_df[!duplicated(ortho_df[,3]),]
                names(ortho_df) <- c('ENSEMBL_MOUSE','CONFIDENCE','ENSEMBL_RAT') 
                ortho_df <- ortho_df %>%
                        select(-CONFIDENCE)
                # Assign the HUGO symbols to a new column
                x <- left_join(x, ortho_df, by = 'ENSEMBL_MOUSE') %>%
                        filter(!is.na(ENSEMBL_RAT)) %>%
                        mutate(SYMBOL_RAT = mapIds(org.Rn.eg.db, 
                                                   ENSEMBL_RAT, 
                                                   'SYMBOL', 'ENSEMBL'))
        }
        # Return the output
        x
}
#########################################################################################

