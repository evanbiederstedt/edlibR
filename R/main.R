#' @useDynLib edlibR
#' @import Rcpp
NULL

#' Align query with target using edit distance
#' 
#' @param query character string Combined with target must have no more than 256 unique values
#' @param target character string Combined with query must have no more than 256 unique values
#' @param mode character string (default="NW") Alignment method to be used. Possible values are:
#'         - 'NW' for global (default). Note that 'NW' stands for 'Needleman-Wunsch'.
#'         - 'HW' for infix. Note that 'HW' stands for 'Hybrid Wunsch'.
#'         - 'SHW' for prefix. Note that 'SHW' stands for 'Semi-Hybrid Wunsch'.
#' @param task character string (default="distance") Specifies what to calculate. The less there is to calculate,
#'         the faster it is. Possible options are (ranked from fastest to slowest):
#'         - 'distance': Find the edit distance and the end locations in the target (default).
#'         - 'locations': Find the edit distance, the end locations, and the start locations.
#'         - 'path': Find the edit distance, the start and end locations, and the alignment path.
#' @param k integer (default=-1) Max edit distance to search for --- the lower this value,
#'         the faster the calculation. Set to -1 (default) to have no limit on edit distance.
#' @param cigarFormat character string (default="extended") Specifies which format to use for writing out the CIGAR string. 
#'         The two possible values are 'standard' and 'extended' (Note: the function getNiceAlignment() only accepts 'cigarFormat="extended"'):
#'         - 'standard': Standard uses the following symbols to generate a CIGAR string: 
#'             Match: 'M', Insertion: 'I', Deletion: 'D', Mismatch: 'M'. 
#'             Note that 'M' in this setting can denote either a sequence match or mismatch.
#'         - 'extended': Extended uses the following symbols to generate a CIGAR string: 
#'             Match: '=', Insertion to target: 'I', Deletion from target: 'D', Mismatch: 'X'.
#'             e.g. CIGAR of "5=1X1=1I" means "5 matches, 1 mismatch, 1 match, 1 insertion (to target)".
#'         For more details on the CIGAR format, please check <http://samtools.github.io/hts-specs/SAMv1.pdf> and <http://drive5.com/usearch/manual/cigar.html>.
#' @param additionalEqualities List of vectors contains pairs of characters (default=NULL) Allows users to extend the definition of equality used in the alignment. 
#'         The input 'additionalEqualities' must be a list of character vectors whereby each character vector 
#'         contains a pair of character strings. (NOTE: the character vectors must contain exactly two strings, a pair.) 
#'         Each pair defines two values as equal. This can be useful e.g. when you want edlib to be case insensitive, 
#'         or if you want certain characters to act as wildcards. If NULL, there will be no additional extensions to edlib's default equality definition.
#' @return List with the following fields:
#'         - editDistance: (integer) The edit distance. This is set to -1 if it is larger than k.
#'         - alphabetLength: (integer) Length of unique characters in 'query' and 'target'
#'         - locations: (list of vectors) List of R vectors of locations, in the format list(c(start, end)). 
#'           Note: if the start or end positions are NULL, this is encoded as 'NA' to work correctly with R vectors.
#'         - cigar: (character string) CIGAR is a standard format for the alignment path.
#'         - cigarFormat: (character string) Format provided by the parameter 'cigarFormat' 
#'           in the function align() which is returned here for the function getNiceAlignment(). 
#'           (Note: the function getNiceAlignment() only accepts 'extended')
#' @examples
#' align("ACTG", "CACTRT", mode="HW", task="path")
#' align("elephant", "telephone")
#' align("ACTG", "CACTRT", mode="HW", task="path", additionalEqualities=list(c("R", "A"), c("R", "G")))
#'
#' @export 
align <- function(query, target, mode="NW", task="distance", k=-1, cigarFormat="extended", additionalEqualities=NULL) {
	if (!string_check(query)) {
		stop("The 'query' input must be a string. Please fix.")
	}

	if (!string_check(target)) {
		stop("The 'target' input must be a string. Please fix.")
	}
	if (!is.numeric(k)) {
		stop("The parameter 'k' must be an integer. Please fix.")
	}
	if (k != as.integer(k)) {
		stop("The parameter 'k' must be an integer. Please fix.")
	}
	'%ni%' <- Negate('%in%')
	if (!is.character(mode)) {
		stop("The parameter 'mode' must be a character string. Please fix.")
	}
	if (!is.character(task)) {
		stop("The parameter 'task' must be a character string. Please fix.")
	}
	if (mode %ni% c("NW", "HW", "SHW")) {
		stop("The parameter 'mode' must match one of these values: 'NW', 'HW', or 'SHW'")
	}
	if (task %ni% c("distance", "locations", "path")) {
		stop("The parameter 'task' must match one of these values: 'distance', 'locations', or 'path'")
	}
	if (cigarFormat %ni% c("extended", "standard")) {
		stop("The parameter 'cigarFormat' must match either the value 'extended' or 'standard'")
	}
	if (!is.null(additionalEqualities)) {
		if (!is.list(additionalEqualities)){
			stop("The 'additionalEqualities' provided must be a list. Please fix.")
		}
		if (length(additionalEqualities) == 0) {
			stop("The 'additionalEqualities' provided has a length 0. It must be a list 
				containing vectors of pairs of characters, e.g. list(c('X', 'Y')) or list(c('A', 'T'), c('C', 'G')). Please fix.")
		}
		if (unique(lengths(additionalEqualities)) != 2) {
			stop("The 'additionalEqualities' provided does not contain pairs of individual
			 characters. It must be a list containing vectors of pairs of characters, e.g. list(c('X', 'Y')) or list(c('A', 'T'), c('C', 'G')). Please fix.")
		}
		if (!all(sapply(additionalEqualities, is.character))){
			stop("The 'additionalEqualities' provided does not contain pairs of individual characters, i.e.,
				character vectors length of 2. Please fix. Note that the parameter 'additionalEqualities' 
				must be a list containing vectors of pairs of characters, e.g. list(c('X', 'Y')) or list(c('A', 'T'), c('C', 'G')).")			
		}
		if (!all(unique(sapply(additionalEqualities, nchar)) == 1)) {
			stop("The 'additionalEqualities' provided must contain pairs of single individual characters, i.e., 
				these cannot be strings but individual characters. Please fix. Note that the parameter 
				'additionalEqualities' must be a list containing vectors of pairs of characters, e.g. list(c('X', 'Y')) or list(c('A', 'T'), c('C', 'G')). 
				Lists of strings are not acceptable, such as list(c('abc', 'def'))")	
		}
	}

	edlibalign(query=query, target=target, input_mode=mode, input_task=task, k=k, cigar_format=cigarFormat, input_additionalEqualities=additionalEqualities)
}




#' Output alignments from align() in NICE format. This outputs the alignment from align() in a visually informative format for human inspection.
#' 
#' @param alignResult list Output of the method align() 
#'     Note: align() requires the argument task="path" for 'alignResult' to output a CIGAR for getNiceAlignment()
#'     Note: Also, align() requires the argument cigarFormat="extended" in order for getNiceAlignment() to work
#' @param query character string The exact query used for alignResult
#' @param target character string The exact target used for alignResult
#' @param gapSymbol character (default="-") Character used to represent gaps in the alignment between query and target. This must be a single character, i.e. a string of length 1.
#' @return Alignment in NICE format, which is an informative visual representation of how the query and target align to each other. 
#'     e.g., for "telephone" and "elephant", it would look like:
#'        telephone
#'         |||||.|.
#'        -elephant
#'     It is represented as an R list with the following fields:
#'       - query_aligned (character string)
#'       - matched_aligned (character string) ('|' for match, '.' for mismatch, ' ' for insertion/deletion)
#'       - target_aligned (character string)
#'     Normally you will want to print these three in order above with the function nice_print(), or another method to apply pretty-printing to R lists
#' @examples
#' query = "elephant"
#' target = "telephone"
#' result = align(query, target, task = "path")
#' nice_algn = getNiceAlignment(result, query, target)
#' 
#' @export 
getNiceAlignment <- function(alignResult, query, target, gapSymbol="-") {
	if (!is.list(alignResult)){
		stop("The object 'alignResult' must be an output of align(). Please check the input 'alignResult'.")
	}
	'%ni%' <- Negate('%in%')
	if ('locations' %ni% names(alignResult)) {
		stop("The object 'alignResult' is expected to have the field 'locations'. Please check the input 'alignResult'.")
	}
	if ('cigar' %ni% names(alignResult)) {
		stop("The object 'alignResult' is expected to have the field 'cigar'. Please check the input 'alignResult'.")
	}
	if (alignResult$cigar == "" || is.null(alignResult$cigar)) {
		stop("The object 'alignResult' contains an empty CIGAR string. 
			Users must run align() with task='path'. Please check the input 'alignResult'.")
	}
	if (alignResult$cigarFormat != "extended") {
		stop("The object 'alignResult' must be generated with the CIGAR format 'extended'. 
			The CIGAR format 'standard' is too ambiguous to  reconstruct the alignment. Please recreate 'alignResult' using align().")
	}

	## check CIGAR correct
	if (!grepl("^[XDI1-9=]*$", alignResult$cigar)) {
		stop("The CIGAR string contains invalid characters, i.e. characters 
			besides 'X', 'I', 'D', '=' or numerals besides 1-9. Please fix.")
	}

	if (!string_check(query)) {
		stop("The 'query' input must be a string. Please fix. Note that that input 'query' must be the exact query used for alignResult.")
	}

	if (!string_check(target)) {
		stop("The 'target' input must be a string. Please fix. Note that that input 'target' must be the exact target used for alignResult.")
	}

	gapsymbol_check <- function(input) {
		return(is.character(input) && length(input)==1 && nchar(input))
	}

	if (!gapsymbol_check(gapSymbol)) {
		stop("The parameter 'gapSymbol' must be a single character. Please fix.")
	}	


	target_pos = alignResult$locations[[1]][1]
	if (is.na(target_pos)) {
		target_pos = 1
	} else if (target_pos==0) {
		## R is 1-based, so this must be 1, not 0
		target_pos = 1
	}
	query_pos = 1
	target_aln = match_aln = query_aln = ""

	cigar_parsing = stringr::str_extract_all(alignResult$cigar, "(\\d+)(\\D)")

	## cigar parsing, motivated by yech1990: https://github.com/Martinsos/edlib/issues/127
	## 'num_occurrences' == Number of occurrences of the alignment operation
	## 'alignment_operation' == Cigar symbol/code that represent an alignment operation

	## check that the strings in the resulting character vector have
	## the format "number of  occurrences" + "alignment_operation"

	for (i in 1:length(cigar_parsing[[1]])) {
		num_occurrences = as.numeric(gsub("([0-9]+).*$", "\\1", cigar_parsing[[1]][i]))
		if (is.na(num_occurrences)) {
			stop("The CIGAR string is in an invalid format. Cannot find number of occurrences for the operation. 
				The expected format should be 'number of occurrences' + 'CIGAR operation', e.g. '4=' or '1D5=1X1=1X'. Please fix.")
		}
		alignment_operation = sub("[[:digit:]]", "", cigar_parsing[[1]][i])
		if (nchar(alignment_operation)>1) {
			stop("The CIGAR string is in an invalid format. Operation detected is not a single character '=' or 'X' or 'D'. 
				The expected format should be 'number of occurrences' + 'CIGAR operation', e.g. '4=' or '1D5=1X1=1X'. Please fix.")       	
		}

		if (alignment_operation == "=") {
			target_aln_addition = substr(target, target_pos, target_pos + num_occurrences - 1)
			target_aln = paste0(target_aln, target_aln_addition)
			target_pos = target_pos + num_occurrences
			query_aln_addition = substr(query, query_pos, query_pos + num_occurrences - 1)
			query_aln = paste0(query_aln, query_aln_addition)
			query_pos = query_pos + num_occurrences
			match_aln_addition = paste(rep("|", num_occurrences), collapse="")
			match_aln = paste0(match_aln, match_aln_addition)
		} else if (alignment_operation == "X") {
			target_aln_addition = substr(target, target_pos, target_pos + num_occurrences - 1)
			target_aln = paste0(target_aln, target_aln_addition)
			target_pos = target_pos + num_occurrences
			query_aln_addition = substr(query, query_pos, query_pos + num_occurrences - 1)
			query_aln = paste0(query_aln, query_aln_addition)
			query_pos = query_pos + num_occurrences
			match_aln_addition = paste(rep(".", num_occurrences), collapse="")
			match_aln = paste0(match_aln, match_aln_addition)
		} else if (alignment_operation == "D") {
			target_aln_addition = substr(target, target_pos, target_pos + num_occurrences - 1)
			target_aln = paste0(target_aln, target_aln_addition)
			target_pos = target_pos + num_occurrences   
			query_aln_addition = paste(rep(gapSymbol, num_occurrences), collapse="")               	
			query_aln = paste0(query_aln, query_aln_addition)
			## query_pos = query_pos + 0
			match_aln_addition = paste(rep(gapSymbol, num_occurrences), collapse="")
			match_aln = paste0(match_aln, match_aln_addition)        	
		} else if (alignment_operation == "I") {
			target_aln_addition = paste(rep(gapSymbol, num_occurrences), collapse="") 
			target_aln = paste0(target_aln, target_aln_addition)	
			## target_pos = target_pos + 0
			query_aln_addition = substr(query, query_pos, query_pos + num_occurrences - 1)
			query_aln = paste0(query_aln, query_aln_addition)
			query_pos = query_pos + num_occurrences
			match_aln_addition = paste(rep(gapSymbol, num_occurrences), collapse="")
			match_aln = paste0(match_aln, match_aln_addition)  
		} else {
			stop("An error occurred due to an incorrect CIGAR format. The CIGAR string contains invalid characters, i.e. characters 
			      besides 'X', 'I', 'D', '=' or numerals besides 1-9. Please fix.")
		}       
	}
	
	return(list('query_aligned'=query_aln, 'matched_aligned'=match_aln, 'target_aligned'=target_aln))
}


#' @keywords internal
string_check <- function(input) {
	return(is.character(input) && length(input)==1)
}


#' Prints the output of getNiceAlignment() in a visually informative format in order to inspect the alignment
#' 
#' @param niceAlignment list Output of the method getNiceAlignment() 
#' @return Pretty-prints the list returned by getNiceAlignment() 
#' @examples
#' query = "elephant"
#' target = "telephone"
#' result = align(query, target, task = "path")
#' nice_algn = getNiceAlignment(result, query, target)
#' nice_print(nice_algn)
#' 
#' @export 
nice_print <- function(niceAlignment){
	if (!is.list(niceAlignment)){
		stop("The object 'niceAlignment' must be an output of getNiceAlignment(). Please check the input 'niceAlignment'.")
	}
	'%ni%' <- Negate('%in%')
	if ('query_aligned' %ni% names(niceAlignment)) {
		stop("The object 'niceAlignment' is expected to have the field 'query_aligned'. Please check the input 'niceAlignment'.")
	}
	if ('matched_aligned' %ni% names(niceAlignment)) {
		stop("The object 'niceAlignment' is expected to have the field 'matched_aligned'. Please check the input 'niceAlignment'.")
	}
	if ('target_aligned' %ni% names(niceAlignment)) {
		stop("The object 'niceAlignment' is expected to have the field 'target_aligned'. Please check the input 'niceAlignment'.")
	}
	## make sure list in specific order, with spacing in names
	nice_algn = list('query:   '=niceAlignment$query_aligned, 'matched: '=niceAlignment$matched_aligned, 'target:  '=niceAlignment$target_aligned)
	for (i in 1:3) {
		print(paste0(names(nice_algn)[i], nice_algn[[i]]))
	}
}


