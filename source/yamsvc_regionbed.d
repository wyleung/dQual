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
import std.math : abs;


// public import sambamba.view;

import bio.bam.reader;
import bio.bam.pileup;
import bio.bam.read;
import bio.bam.writer;

import yamsvc.utils;

void printBedRegion(string reference_name, ref DiscordantRegion bed ) {
    if ( bed.length() && bed.maxCoverageDiscordant() ) {
        writefln("%s\t%d\t%d\tDP=%d;DPdis=%d;CLEN=%d", 
            reference_name, 
            bed.start, 
            bed.end+1, 
            bed.avgCoverage(), 
            bed.maxCoverageDiscordant(), 
            bed.length()
        );
    }
}

void makeRegion(BamReader bam, 
				ReferenceSequenceInfo reference_sequence, 
				ushort cutoff_min, 
				ushort cutoff_max, 
				ushort window_width,
				uint estimated_insertsize, 
				uint estimated_insertsize_stdev) {
    auto reads = bam[ cast(string) reference_sequence.name() ][1 .. reference_sequence.length];
    auto pileup = makePileup(reads, true);

    int lastChrom = -1;
        
    DiscordantRegion bed = new DiscordantRegion();
        
    foreach (column; pileup) {
        auto total_coverage = 0;
        auto normal_coverage = 0;
        auto discordant_coverage = 0;
        
        foreach (BamRead read; column.reads) {
        	total_coverage++;
        	
            if( is_concordant(read) && has_minimum_quality(read, 20) ) {
                normal_coverage++;
            }
            if ( isDiscordant(read, estimated_insertsize, estimated_insertsize_stdev ) ) {
                discordant_coverage++;
            }
        }

        // should be dynamic range
        if ( discordant_coverage >= cutoff_min && 
            discordant_coverage <= cutoff_max ) {
            // check whether to output or not
            // quick patch for initializing the lastChrom
            if ( lastChrom == -1 ) {
                lastChrom = column.ref_id;
                bed.chromosome = cast(short) column.ref_id;
            }

            // We are working on the same chromosome as before
            bool same_chrom = bed.chromosome == column.ref_id;
            bool in_window = (abs(column.position - window_width) < bed.end) && (bed.end < column.position);
            
            if ( same_chrom && in_window ) {
                bed.end = column.position;
                bed.addCoverage( cast(ushort) normal_coverage );
                bed.addCoverageDiscordant( cast(ushort) discordant_coverage );
            }
            else {
                if ( bed.working ) {
                    // report only if we were working on a bed track.
                    // and if we have logical results
                    printBedRegion( bam.reference_sequences[bed.chromosome].name, bed );
                    
                }

                // start new bed track for the new region
                bed.startReporting( cast(ushort) column.ref_id, column.position, column.position );
                
            }
        }
    }
}

void makeRegion(string bamfile, 
				ReferenceSequenceInfo reference_sequence, 
				ushort cutoff_min, 
				ushort cutoff_max, 
				ushort window_width,
				uint estimated_insertsize, 
				uint estimated_insertsize_stdev) {
    auto tp = new TaskPool(2);
    scope(exit) tp.finish();
    auto bam = new BamReader(bamfile, tp);
    makeRegion(bam, 
    			reference_sequence, 
				cutoff_min, 
				cutoff_max, 
				window_width,
				estimated_insertsize, 
				estimated_insertsize_stdev);
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
    
    // free some threads for the bamreader, which also uses a Taskpool for reading operations11
    n_threads = n_threads > 5 ? n_threads - 2 : n_threads;
    // this a dangerous operation, what if we still overuse the amount of 'threads' on the host sytem?


    auto task_pool = new TaskPool(n_threads);
    scope(exit) task_pool.finish();
    
    string bamfile = args[1];
    auto src = new BamReader(bamfile);
    
    foreach (refseq; src.reference_sequences) {
//      writefln("Starting thread for %s", refseq);
        auto t = task!makeRegion(bamfile, 
        						refseq, 
    							cutoff_min, 
    							cutoff_max, 
    							window_width,
    							estimated_insertsize, 
    							estimated_insertsize_stdev);
        task_pool.put(t);
        }

    return 0;
}
