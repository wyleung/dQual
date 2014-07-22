#!/usr/bin/env rdmd

/*
 @author: Wai Yi Leung
 @copyright: 2013-2014 Wai Yi Leung
 @contact: w.y.leung@e-sensei.nl
 @license: MIT Licence
*/

/*
    Region definition: a genomic region with at least 2 bases
    We summarize the genome regions which have genomic DNA reads aligned, 
    but they are reported discordant (abnormal insert-size, mate on other chromosome)
    
    
*/


import std.getopt;
import std.stdio;
import std.parallelism;
import std.algorithm : count, max, reduce;

// public import sambamba.view;

import bio.bam.reader;
import bio.bam.pileup;
import bio.bam.read;

class DiscordantRegion {
	ulong start=0;
	ulong end=0;
	short chromosome=-1;
	ushort[] coverages;
	bool working;
	
	void addCoverage( ushort coverage ) {
		// append coverage to the coverages
		coverages ~= coverage;
	}
	
	ushort maxCoverage(){
		if(coverages.length){
			return reduce!(max)(coverages);
		} else {
			return 0;
		}
	}
	
	ushort avgCoverage(){
		/* FIXME: cast from ulong to int is not so performing? 
			We are only interested in coverages upto ushort (65,535) large?
		*/
		auto sum = reduce!((a,b) => a + b)(0, coverages);
		return cast(ushort)(sum / coverages.length);
	}
	
	void resetCoverage() {
		coverages = [];
	}
	
	ulong length() {
		return end-start;
		
	}
	
	void startReporting( ) {
		working = true;
	}
}


uint n_threads;

bool is_concordant( BamRead read ) {
	switch ( read.flag() ) {
		case 99:
		case 47:
		case 83:
		case 163:
			return true;
		default:
		return false;
	}
	return false;
}

bool has_minimum_quality( BamRead read, int minQuality ) {
	if(read.mapping_quality == 0) return false;
	return minQuality >= read.mapping_quality;
}


// auto makePileup(R)(R reads, int window_width, int cutoff_min, int cutoff_max ) {
// 	auto pileup = makePileup(reads, true);
// 	
// 	int lastChrom = -1;
// 	ulong lastPos = 0;
// 	ulong lastPosStart = 0;
// 
// 	DiscordantRegion bed = new DiscordantRegion();
// 
// 	foreach (column; pileup) {
// 		auto coverage = 0;
// 
// 		foreach (read; column.reads) {
// 			if( is_concordant(read) ) coverage++;
// 		}
// 		auto discordant_coverage = column.coverage-coverage;
// 		// should be dynamic range
// 		if ( cutoff_min <= discordant_coverage ) {
// 			// check whether to output or not
// 			// quick patch
// 			if ( lastChrom == -1 ) {
// 				lastChrom = column.ref_id;
// 				bed.chromosome = column.ref_id;
// 			}
// 			
// 			if ( bed.chromosome == column.ref_id && (column.position - window_width) < bed.end && bed.end < column.position ) {
// 				bed.end = column.position;
// 				bed.addCoverage( discordant_coverage );
// 			}
// 			else {
// 				if ( bed.working ) {
// 					writefln("%s\t%d\t%d\tDP=%d\tDPdis=%d\tCLEN=%d", bam.reference_sequences[bed.chromosome].name, bed.start+1, bed.end, coverage, bed.maxCoverage(), bed.length());
// 				}
// 				bed.startReporting();
// 				bed.chromosome = column.ref_id;
// 				bed.start = column.position;
// 				bed.end = column.position;
// 				bed.resetCoverage();
// 			}
// 		}
// 	}
// }

void printUsage() {
    stderr.writeln("Usage: yamsvp-bedregion [options] <input.bam | input.sam>");
    stderr.writeln();
    stderr.writeln("Options: -T, --threads=NTHREADS");
    stderr.writeln("                    maximum number of threads to use");
    stderr.writeln("         -w, --window=WINDOWWIDTH");
    stderr.writeln("                    width of window (bp) to scan [200]");
    stderr.writeln("         -l, --cutoff_min");
    stderr.writeln("                    low threshold for reads covering a breakpoint [4]");
    stderr.writeln("         -h, --cutoff_max");
    stderr.writeln("                    high threshold for reads covering a breakpoint [90]");
    stderr.writeln("         -v, --verbose");
    stderr.writeln("                    Turn on verbose mode");
}

int main(string[] args) {
	int n_threads = 4;
	string bamFile;
	bool verbose;
	int window_width = 200;
	int cutoff_max = 90;
	int cutoff_min = 4;

    getopt(
        args,
        std.getopt.config.caseSensitive,
        "threads|T",  &n_threads,    // numeric
        "window|w", &window_width,
        "cutoff_max|h", &cutoff_max,
        "cutoff_min|l", &cutoff_min,
        "verbose|v", &verbose,   // flag
        );

    if (args.length < 2) {
        printUsage();
        return 0;
    }
    
    auto task_pool = new TaskPool(n_threads);
    scope(exit) task_pool.finish();
    
    auto bam = new BamReader(args[1], task_pool);
//	 foreach (refseq; bam.reference_sequences) {
//	 	writeln( refseq.name );
//	 }

	foreach (refseq; bam.reference_sequences) {
		// auto reads = bam[ refseq ].reads;
		auto reads = bam[ refseq.name ][1 .. refseq.length];
//		findDiscordantRegions( reads, window_width, cutoff_min, cutoff_max );

//		 auto reads = bam["1"][41_300 .. 102_800];
		auto pileup = makePileup(reads, true);

		int lastChrom = -1;
		ulong lastPos = 0;
		ulong lastPosStart = 0;
		
		DiscordantRegion bed = new DiscordantRegion();
		
		foreach (column; pileup) {
			auto coverage = 0;
		
			foreach (read; column.reads) {
				// writefln("%30s\t%s\t%.2d\t%s\t%2s/%2s\t%2s/%2s\t%10s \t%s\tcc_%s_cc\t%s \t%s \t%s %s", 
				// 	read.name, 
				// 	read.current_base,
				// 	read.current_base_quality,
				// 	read.cigar_operation,
				// 	read.cigar_operation_offset + 1, read.cigar_operation.length,
				// 	read.query_offset + 1, read.sequence.length,
				// 	read.cigarString(),
				// 	read.flag(),
				// 	is_concordant(read),
				// 	read.is_paired(),
				// 	read.proper_pair(),
				// 	read.cigar_before, read.cigar_after);
				if( is_concordant(read) && has_minimum_quality(read, 30) ) coverage++;
			}
			auto discordant_coverage = column.coverage-coverage;
			// should be dynamic range
			if ( cutoff_min <= discordant_coverage ) {
				// check whether to output or not
				// quick patch
				if ( lastChrom == -1 ) {
					lastChrom = column.ref_id;
					bed.chromosome = cast(short) column.ref_id;
				}

				if ( bed.chromosome == column.ref_id && (column.position - window_width) < bed.end && bed.end < column.position ) {
					bed.end = column.position;
					bed.addCoverage( cast(ushort) discordant_coverage );
				}
				else {
					if ( bed.working ) {
						writefln("%s\t%d\t%d\tDP=%d;DPdis=%d;CLEN=%d", bam.reference_sequences[bed.chromosome].name, bed.start, bed.end+1, coverage, bed.maxCoverage(), bed.length());
					}
					bed.startReporting();
					bed.chromosome = cast(ushort) column.ref_id;
					bed.start = column.position;
					bed.end = column.position;
					bed.resetCoverage();
				}
			}
		}
	}
	return 0;
}
