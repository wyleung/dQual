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


import std.algorithm;
import std.getopt;
import std.math;
import std.parallelism;
import std.stdio;
import std.regex;
import std.range;

// public import sambamba.view;

import bio.bam.reader;
import bio.bam.pileup;
import bio.bam.read;
import bio.bam.writer;

import sambamba.utils.common.filtering;

import yamsvc.utils;
import yamsvc.datatypes;

void printBedRegion(string reference_name, ref DiscordantRegion bed ) {
//    writefln("%s\t%d\t%d\tDP=%d;DPdis=%d;CLEN=%d;ORIENTATION=%d;MATEPOS=%s;MATECHR=%s", 
//        reference_name, 
//        bed.start, 
//        bed.end+1, 
//        bed.avgCoverage(), 
//        cast(ulong) bed.cov_disc / cast(ulong) bed.size, 
//        bed.size,
//        bed.orientation,
//        bed.mate_positions,
//        bed.mate_chromosome
//    );
    writefln("%s\t%d\t%d\tDP=%d;DPdis=%d;CLEN=%d;ORIENTATION=%d", 
        reference_name, 
        bed.start, 
        bed.end+1, 
        bed.avgCoverage(), 
        cast(ulong) bed.cov_disc / cast(ulong) bed.size, 
        bed.size,
        bed.orientation
    );
}

void makeRegion(BamReader bam, 
                ReferenceSequenceInfo reference_sequence, 
                ushort cutoff_min, 
                ushort cutoff_max, 
                ushort window_width,
                uint estimated_insertsize, 
                uint estimated_insertsize_stdev,
                uint max_sd_degree, 
                TaskPool tp) {
    auto reads = bam[ cast(string) reference_sequence.name() ][1 .. reference_sequence.length];
    auto pileup = makePileup(reads, true);

    int lastChrom = -1;

//    string bam_outputfilename = "/tmp/" ~ cast(string)(reference_sequence.name()) ~ ".bam";

//    DiscordantRegion bed = new DiscordantRegion();
//    BamRegion[] regions;
//    ulong[] reads_seen;

//    auto name_regex = regex(`[.*]:(.*)$`, "g");

	// 
    DiscordantRegion[ubyte] regions;
    regions[0] = new DiscordantRegion(0);
    regions[1] = new DiscordantRegion(1);
    regions[2] = new DiscordantRegion(2);
    regions[3] = new DiscordantRegion(3);


    foreach (column; pileup) {
        auto total_coverage = 0;
        auto normal_coverage = 0;
        auto discordant_coverage = 0;
        
        foreach (BamRead read; column.reads) {
            total_coverage++;
            
            ubyte orientation = 0;
            
            if (read.is_first_of_pair && !read.is_reverse_strand) {
            	// F1
            	orientation = 0;
            } else if(read.is_first_of_pair && read.is_reverse_strand) {
            	// R1
            	orientation = 1;
        	} else if (read.is_second_of_pair && !read.is_reverse_strand) {
        		// F2
        		orientation = 2;
        	} else if (read.is_second_of_pair && read.is_reverse_strand) {
        		// R2
        		orientation = 3;
        	}

        	// work on this bedregion
            DiscordantRegion bed = regions[ orientation ];
//            writeln(bed);
            // check whether to output or not
            // quick patch for initializing the lastChrom
            if ( lastChrom == -1 ) {
                lastChrom = column.ref_id;
                bed.chromosome = cast(short) column.ref_id;
            }

            // We are working on the same chromosome as before
            bool same_chrom = bed.chromosome == column.ref_id;
            bool in_window = (min(abs(column.position - window_width), column.position) <= bed.end);
            
            if ( same_chrom && in_window ) {
	            if( is_concordant(read) && has_minimum_quality(read, 20) ) {
	                bed.addCoverage( 1 );
	            }
	            if ( isDiscordant(read, estimated_insertsize, estimated_insertsize_stdev ) ) {
	                bed.addCoverageDiscordant( 1 );
	                bed.end = column.position + read.basesCovered;
	                
	                bed.add_mate_pos( read.mate_position );
	                bed.add_matechromosome( read.mate_ref_id );
	                
//	                writeln( read.template_length ); 
	            }
            }
            else {
            	version(Debug){
//            	if ( same_chrom ) {
//            		writeln("Same chromosome");
//        		}
//            	if ( in_window ) {
//            		writeln("In window");
//        		} else {
//        			writeln((min(abs(column.position - window_width), column.position) <= bed.end));
//        			writeln((min(abs(column.position - window_width), column.position)));
//        			writeln(bed.end <= column.position);
//            		writefln("%d %d", column.position, bed.end);
//        		}
            	}
            	auto bp_perbase = cast(ulong)(bed.cov_disc+1) / cast(ulong)(bed.size+1);
            	if ( bp_perbase  >= cutoff_min &&
            		 bp_perbase <= cutoff_max ) {

//writefln("Dbases: %d, clen: %d, DP/base: %d", bed.cov_disc, bed.size, bp_perbase);

	                if ( bed.working ) {
	                    printBedRegion( bam.reference_sequences[bed.chromosome].name, bed );
	                }
                }
                bed.startReporting( cast(ushort) column.ref_id, read.position, read.position+read.basesCovered );
            }

//
//            if( is_concordant(read) && has_minimum_quality(read, 20) ) {
//                normal_coverage++;
//                bed.addCoverage( 1 );
//            }
//            if ( isDiscordant(read, estimated_insertsize, estimated_insertsize_stdev ) ) {
//                discordant_coverage++;
//                bed.addCoverageDiscordant( 1 );             
//            }
//            
            
        }

//        // should be dynamic range
//        if ( discordant_coverage >= cutoff_min && 
//            discordant_coverage <= cutoff_max ) {
//            // check whether to output or not
//            // quick patch for initializing the lastChrom
//            if ( lastChrom == -1 ) {
//                lastChrom = column.ref_id;
//                bed.chromosome = cast(short) column.ref_id;
//            }
//
//            // We are working on the same chromosome as before
//            bool same_chrom = bed.chromosome == column.ref_id;
//            bool in_window = (abs(column.position - window_width) < bed.end) && (bed.end < column.position);
//            
//            if ( same_chrom && in_window ) {
//                bed.end = column.position;
//                bed.addCoverage( cast(ushort) normal_coverage );
//                bed.addCoverageDiscordant( cast(ushort) discordant_coverage );
//            }
//            else {
//                if ( bed.working ) {
//                    // report only if we were working on a bed track.
//                    // and if we have logical results
//                    printBedRegion( bam.reference_sequences[bed.chromosome].name, bed );
////                    regions ~= BamRegion(cast(uint) bed.chromosome,
////                                 cast(uint) bed.start, cast(uint) bed.end);
//                    
//                }
//
//                // start new bed track for the new region
//                bed.startReporting( cast(ushort) column.ref_id, column.position, column.position );
//                
//            }
//        }
    }
    
    /* Continu with calling the discordant reads, now we know the regions to scan for */
    
//    std.algorithm.sort(regions);
/*
    Filter read_filter = new NullFilter();
    read_filter = new AndFilter(read_filter, new ValidAlignmentFilter());

    Filter only_mapped = new AndFilter(
        new NotFilter(new FlagFilter!"is_unmapped"()), 
        new NotFilter(new FlagFilter!"mate_is_unmapped"())
    );

    Filter quality = new IntegerFieldFilter!">="("mapping_quality", 20);

    Filter flags = new NotFilter(
        new OrFilter( new AlignmentFlagFilter!"=="(99), 
        new OrFilter( new AlignmentFlagFilter!"=="(147), 
        new OrFilter( new AlignmentFlagFilter!"=="(83), 
        new AlignmentFlagFilter!"=="(163) ) ) )
    );
    
    Filter template_size = new OrFilter( 
        new TemplateSizeFilter!">="( (estimated_insertsize + ( max_sd_degree * estimated_insertsize_stdev )) ),
        new TemplateSizeFilter!"<="( (estimated_insertsize - ( max_sd_degree * estimated_insertsize_stdev )) )
        );
        
    
    read_filter = new AndFilter(read_filter, only_mapped);
    read_filter = new AndFilter(read_filter, quality);
    read_filter = new AndFilter(read_filter, new OrFilter( flags, template_size )); 
 
 
    int processAlignments(P)(P processor) {
    
        void runProcessor(SB, R, F)(SB bam, R reads, F filter) {
            if (processor.is_serial)
                bam.assumeSequentialProcessing();
            if (cast(NullFilter) filter)
                processor.process(reads, bam);
            else
                processor.process(reads.filtered(filter), bam);
        }
        
        // adapted from sambamba
        auto _reads = bam.getReadsOverlapping(regions);
        runProcessor(bam, _reads, read_filter);
        
        foreach( BamRead read; _reads ){
            CallRegion cr = new CallRegion();
        }
        
        
        return 0;
    }

    processAlignments(new BamSerializer(bam_outputfilename, cast(int) 9, tp));
*/
}

void makeRegion(string bamfile, 
                ReferenceSequenceInfo reference_sequence, 
                ushort cutoff_min, 
                ushort cutoff_max, 
                ushort window_width,
                uint estimated_insertsize, 
                uint estimated_insertsize_stdev,
                uint max_sd_degree) {
    auto tp = new TaskPool(4);
    scope(exit) tp.finish();
    auto bam = new BamReader(bamfile, tp);
    auto input_buf_size = 64_000_000;
    bam.setBufferSize(input_buf_size);
    makeRegion(bam, 
                reference_sequence, 
                cutoff_min, 
                cutoff_max, 
                window_width,
                estimated_insertsize, 
                estimated_insertsize_stdev,
                max_sd_degree, tp);
}

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
    stderr.writeln("         -i, --insertsize");
    stderr.writeln("                    insertsize of the input bam [450]");
    stderr.writeln("         -s, --insertsd");
    stderr.writeln("                    stdev in the insertsize [15]");
    stderr.writeln("         -v, --verbose");
    stderr.writeln("                    Turn on verbose mode");
}

int main(string[] args) {

    int n_threads = std.parallelism.totalCPUs;
    ushort window_width = 200;
    ushort cutoff_max = 90;
    ushort cutoff_min = 4;
    uint estimated_insertsize = 450;
    uint estimated_insertsize_stdev = 15;
    uint max_sd_degree = 2;
    bool verbose;

    getopt(
        args,
        std.getopt.config.caseSensitive,
        "threads|T",  &n_threads,    // numeric
        "window|w", &window_width,   // numeric
        "cutoff_max|h", &cutoff_max, // numeric
        "cutoff_min|l", &cutoff_min, // numeric
        "insertsize|i", &estimated_insertsize, // numeric
        "insertsd|s", &estimated_insertsize_stdev, // numeric
        "verbose|v", &verbose,       // flag
        );

    if (args.length < 2) {
        printUsage();
        return 0;
    }

    auto task_pool = new TaskPool(n_threads);
    scope(exit) task_pool.finish();
    
    string _bamfilename = args[1];
    auto src = new BamReader(_bamfilename);
    
    foreach (refseq; src.reference_sequences) {
//      writefln("Starting thread for %s", refseq);
        auto t = task!makeRegion(_bamfilename, 
                                refseq, 
                                cutoff_min, 
                                cutoff_max, 
                                window_width,
                                estimated_insertsize, 
                                estimated_insertsize_stdev, 
                                max_sd_degree);
        task_pool.put(t);
        }

    return 0;
}
