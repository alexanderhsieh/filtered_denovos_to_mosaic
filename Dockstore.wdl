version 1.0

import "https://raw.githubusercontent.com/alexanderhsieh/EM-mosaic-pipeline/master/tasks/tasks_filtered_denovo_to_mosaic.wdl" as call_mosaic

workflow filtered_denovo_to_mosaic {

	input {
		File dnsnvs
		Int postcut 
		Int sample_avg_dp
		Int cohort_size
		String output_prefix
	}

	#run EM_mosaic
	call call_mosaic.detect_mosaic {
		input:
			infile = dnsnvs,
			postcut = postcut,
			cohort_size = cohort_size,
			sample_avg_dp = sample_avg_dp,
			outprefix = output_prefix
	}

	output {
		File denovos = detect_mosaic.denovos
		File mosaics = detect_mosaic.mosaics
		Array[File] output_plots = detect_mosaic.output_plots
	}

}

