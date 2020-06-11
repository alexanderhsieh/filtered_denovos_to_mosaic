## Copyright Broad Institute, 2020
## This script takes as input a set of filtered and QC'd de novo SNVs and run the following steps
## (1) combine filtered de novo callsets from each sample into a cohort-wide callset
## (2) remove variants belonging to samples with an abnormally high number of de novos (outliers)
## (3) run EM-mosaic to annotate de novo variants with additional columns (p-value, 
##     likelihood ratio, and posterior odds) and produce a candidate mosaics file (variants with 
##     posterior odds above a given cutoff; default=10), as well as associated plots
## 
## TESTED: 
## Versions of other tools on this image at the time of testing:
##
## LICENSING : This script is released under the WDL source code license (BSD-3) (see LICENSE in https://github.com/broadinstitute/wdl). 
## Note however that the programs it calls may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running this script. 
## Please see the docker for detailed licensing information pertaining to the included programs.
##

###########################################################################
#WORKFLOW DEFINITION
###########################################################################
workflow filtered_denovo_to_mosaic {
  
  Array[File] dnsnvs 

  File em_mosaic_script
  File outlier_script
  
  Int postcut
  Int cohortsize # number of samples in this cohort
  Int expected_dnsnvs # rough expectation; exomes: ~1 denovo/sample; genomes: ~100 denovos/sample

  String output_prefix

  parameter_meta{
    dnsnvs: "filtered de novo SNVs file with cols {id, chr, pos, ref, alt, refdp, altdp} at minimum"
    em_mosaic_script: "EM-mosaic.R script"
    outlier_script: "filter_outlier.R script"
    postcut: "posterior odds cutoff score (default: 10) for defining mosaic vs. germline"
    cohortsize: "number of samples in this cohort; for outlier calculation"
    expected_dnsnvs: "rough expected number of de novos/sample; exomes=1, genomes=100"
    output_prefix: "filename prefix for output files"
  }

  meta{
    author: "Alex Hsieh"
    email: "ahsieh@broadinstitute.org"
  }

  # gather sample-level callsets into a cohort-level callset
  call gather_callsets {
    input:
    callsets = dnsnvs,
    outprefix = output_prefix
  }


  output {
      File cohort_callset = gather_callsets.out
    }

}


###########################################################################
#Task Definitions
###########################################################################

# Gathers sample-level callsets into a cohort-level callset
task gather_callsets {
  Array[File] callsets
  String outprefix

  String outfname = "ADfile.${outprefix}.txt"
 
  command {

    while read file; do
      cat $file >> "tmp.cohort.denovos.txt"
    done < ${write_lines(callsets)};

    grep "^id" "tmp.cohort.denovos.txt" | head -n 1 > "header.txt"

    (cat header.txt; grep -v "^id" "tmp.cohort.denovos.txt") > ${outfname}
  }

  runtime {
    docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
  }

  output {
    File out = "${outfname}"
  }
}

#Adds Outlier filter-related columns to variants file for downstream filtering
## NOTE:  THIS IS A COHORT-LEVEL FILTER
task filter_OUT {
  File infile
  File script
  Int cohort_size
  Int cutoff
  String outprefix = basename(infile, '.txt')

  command {
    Rscript ${script} ${infile} ${outprefix} ${cohort_size} ${cutoff}
  }

  runtime {
    docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
  }

  output {
    File outfile = "${outprefix}.OUT.txt"
  }
}

#Scores variants in dnSNVs file
task detect_mosaic {
  File infile
  File script
  Int postcut 
  String outprefix = basename(infile, '.txt')

  command {
    Rscript ${script} ${infile} ${outprefix} ${postcut} 
  }

  runtime {
    docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
  }

  output {
    File denovos = "${outprefix}.denovo.txt"
    File mosaics = "${outprefix}.candidates.txt"
    File plot_EM = "${outprefix}.EM.pdf"
    File plot_QQ = "${outprefix}.QQ.pdf"
    File plot_overdispersion = "${outprefix}.overdispersion.pdf"
    File plot_dp_vs_vaf = "${outprefix}.dp_vs_vaf.pdf"
    File plot_vaf_vs_post = "${outprefix}.vaf_vs_post.pdf"
  }
}