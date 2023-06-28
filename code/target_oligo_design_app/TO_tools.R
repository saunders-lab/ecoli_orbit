

get_replichore <- function(pos, ori = 3923882.5, ter = 1590250.5 ){
  

    # Determine the replichore of a bacterial chromosome for a certain position. Requires origin and terminus positions. Assumes E. coli like organization.
    # 
    # pos : int
    #     Genomic coordinate of interest.
    # ori : float
    #     Genomic coordinate of the origin of replication. 
    # ter : float
    #     Genomic coordinate of the replication terminus.

  
  pos = as.integer(pos)
  
  if((pos<0) | (pos>4641652)){
    stop("position must be within genome.")
  }
  
  if((pos > ori) | (pos<ter)){
    rep = 1
  } 
  
  if((pos<ori) & (pos>ter)){
    rep = 2
  }
   
    
  return(rep)
  
}


get_target_oligo <- function(left_pos, right_pos, genome, homology = 90, attB_dir = '+', attB_fwd_seq = 'ggcttgtcgacgacggcggtctccgtcgtcaggatcat',  verbose = F, html_format = F, surround = 0){

    # Given a set of parameters, get an ORBIT oligo that targets the lagging strand. 
    # Left and right positions are absolute genomic coordinates that specify the final nucleotides to keep unmodified in the genome, 
    # everything in between will be replaced by attB. In other words the left position nucleotide is the final nt before attB in the oligo.
    # The right position nt is the first nt after attB in the oligo.
    # 
    # This function determines the lagging strand by calling `get_replichore()` on the left_pos.
    # Typically attB_dir should be set to the same direction as the gene of interest, such that the integrating plasmid will insert with payload facing downstream.
    # attB_fwd_seq can be modified, and the total homology can be modified, but should be an even number since homology arms are symmetric. 
    # 
    # Verbose prints helpful statements for testing functionality.
    # 
    # Parameters
    # -----------------
    # left_pos : int
    #     Left genomic coordinate of desired attB insertion. attB is added immediately after this nt.
    # right_pos : int
    #     Right genomic coordinate of desired attB insertion. attB is added immediately before this nt.
    # genome : str
    #     Genome as a string.
    # homology : int (even)
    #     Total homology length desired for oligo. Arm length = homology / 2.
    # attB_dir : chr ('+' or '-')
    #     Desired direction of attB  based on genomic strand. Typically same direction as gene.
    # attB_fwd_seq : str
    #     Sequence of attB to insert between homology arms.
    # verbose : bool
    #     If true, prints details about genomic positions and replichore.
    # Returns
    # ---------------
    # oligo : str
    #     Targeting oligo against lagging strand, including the attB sequence in the correct orientation.


left_pos <-  as.integer(left_pos)

right_pos <-  as.integer(right_pos)

# Arm length is 1/2 total homology. Arms are symmetric
arm_len <-  as.integer(homology / 2)

# As written, should be inclusive with indexing from 1
left_arm <-  str_sub(genome,(left_pos - arm_len + 1),left_pos) #rewrite with str_sub()

right_arm <- str_sub(genome,right_pos,(right_pos - 1 + arm_len))

surround_len <- as.integer(surround/2)

left_start_lab <- (left_pos - arm_len + 1)-surround_len

right_start_lab <- (right_pos - 1 + arm_len)+surround_len

left_surround <- str_sub(genome,left_start_lab, (left_pos - arm_len + 1)-1)

right_surround <- str_sub(genome,(right_pos - 1 + arm_len)+1,(right_pos - 1 + arm_len)+surround_len )

# Generate attB reverse sequence
#as.character(seq_complement(seq_reverse(dna(sticky_ends))))

#seq_attB <-  dna(attB_fwd_seq)
attB_rev_seq = tolower(as.character(seq_complement(seq_reverse(dna(attB_fwd_seq)))))

# Replichore 1
if(get_replichore(left_pos) == 1){
  
  rep = 1

# Reverse complement replichore 1 sequences.
#left_arm_seq = Seq(left_arm)
left_arm_prime <- as.character(seq_complement(seq_reverse(dna(left_arm))))

left_surround_prime <- tolower(as.character(seq_complement(seq_reverse(dna(left_surround)))))

#right_arm_seq = Seq(right_arm)
right_arm_prime <-  as.character(seq_complement(seq_reverse(dna(right_arm))))

right_surround_prime <-  tolower(as.character(seq_complement(seq_reverse(dna(right_surround)))))



# Determine attB direction and paste fwd/rev seq accordingly
if(attB_dir == '+'){
  
  oligo <-  paste0(right_arm_prime, attB_rev_seq, left_arm_prime)
  
  if(html_format == T){
    oligo <- color_html_strings(c(right_surround_prime,right_arm_prime, attB_rev_seq, left_arm_prime, left_surround_prime),
                                c("#FFFFFF", "#56B4E9","#D55E00","#E69F00","#FFFFFF"),
                                start_lab = right_start_lab, end_lab = left_start_lab)
  }

}
if(attB_dir == '-'){
  
  oligo <- paste0(right_arm_prime,attB_fwd_seq,left_arm_prime)
  
  if(html_format == T){
    oligo <- color_html_strings(c(right_surround_prime,right_arm_prime,attB_fwd_seq,left_arm_prime, left_surround_prime),
                                c("#FFFFFF", "#56B4E9","#D55E00","#E69F00","#FFFFFF"),
                                start_lab = right_start_lab, end_lab = left_start_lab)
  }
}
}


# Replichore 2
if(get_replichore(left_pos) == 2){
  
  rep = 2

# '+' arm sequence used. Determine attB direction and paste accordingly.
if(attB_dir == '+'){
  
  oligo <-  paste0(left_arm, attB_fwd_seq, right_arm)
  
  if(html_format == T){
    oligo <- color_html_strings(c(left_surround, left_arm, attB_fwd_seq, right_arm, right_surround),
                                c("#FFFFFF","#E69F00" ,"#D55E00","#56B4E9","#FFFFFF"),
                                start_lab = left_start_lab, end_lab = right_start_lab
                                )
  }
  
}
if(attB_dir == '-'){
  
  oligo <- paste0(left_arm, attB_rev_seq, right_arm)
  
  if(html_format == T){
    oligo <- color_html_strings(c(left_surround,left_arm, attB_rev_seq, right_arm, right_surround),
                                c("#FFFFFF","#E69F00" ,"#D55E00","#56B4E9","#FFFFFF"),
                                start_lab = left_start_lab, end_lab = right_start_lab)
  }
}
}

# Verbose print statements
if(verbose){
  
print(paste0('left_arm_coord = ', left_pos - arm_len +1,' : ', left_pos))
print(paste0('right_arm_coord = ', right_pos, ' : ', right_pos -1 + arm_len))
print(paste0('Replichore = ', rep))
}

return(oligo)
  
}

color_html_strings <- function(strings, colors = c("#E69F00","#D55E00","#56B4E9"), start_lab = '', end_lab = ''){
  #strings and colors are vectors of the same length

  color_strings = c()
  
  for(i in 1:length(strings)){
    color_strings[i] <- paste0('<span style=\"word-wrap: break-word;background-color: ',colors[i],'\">',strings[i],'</span>')
  }

  paste0(c("[", start_lab, "]  ", "5'  ",color_strings,"  3'", "  [",end_lab, "]"), collapse = '')

}

