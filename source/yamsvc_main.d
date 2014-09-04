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


import core.atomic;

import std.algorithm;
import std.getopt;
import std.math;
import std.parallelism;
import std.stdio;

import bio.bam.reader;
import bio.bam.pileup;
import bio.bam.read;
import bio.bam.writer;

import sambamba.utils.common.filtering;

import bam.SortedBamReader;
import yamsvc.datatypes;
import yamsvc.utils;

void printBedRegion(string reference_name, ref DiscordantRegion bed ) {
    if ( bed.size && bed.maxCoverageDiscordant() ) {
        writefln("%s\t%d\t%d\tDP=%d;DPdis=%d;CLEN=%d", 
            reference_name, 
            bed.start, 
            bed.end+1, 
            bed.avgCoverage(), 
            bed.maxCoverageDiscordant(), 
            bed.size
        );
    }
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

    string bam_outputfilename = "/tmp/" ~ cast(string)(reference_sequence.name()) ~ ".bam";
        
    DiscordantRegion bed = new DiscordantRegion(0);
    BamRegion[] regions;

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
                    regions ~= BamRegion(cast(uint) bed.chromosome,
                                 cast(uint) bed.start, cast(uint) bed.end);
                    
                }

                // start new bed track for the new region
                bed.startReporting( cast(ushort) column.ref_id, column.position, column.position );
                
            }
        }
    }
    
    /* Continu with calling the discordant reads, now we know the regions to scan for */
    std.algorithm.sort(regions);

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
            /* Do indept analysis */
            CallRegion cr = new CallRegion();
            
        }
        return 0;
    }

    processAlignments(new BamSerializer(bam_outputfilename, cast(int) 9, tp));
    auto src = new BamReader(bam_outputfilename, tp);
    // Fix index, we need this for random access
    if( !src.has_index ) {
        src.createIndex();
    }

}

void extractReadInfo( BamReader bam, 
                ReferenceSequenceInfo reference_sequence, 
                ushort cutoff_min, 
                ushort cutoff_max, 
                ushort window_width,
                uint estimated_insertsize, 
                uint estimated_insertsize_stdev,
                uint max_sd_degree, 
                TaskPool tp) {
    auto reads = bam[ cast(string) reference_sequence.name() ][0 .. reference_sequence.length];
    
    ReadInfo[] contigreads;
    BamRead read;

    // test how mush storage would be neede to store whole genome
    while (!reads.empty) {
        reads.popFront();
//        read = reads.front;
        read = reads.front.dup; // duplicate only for low memory usage ?
        if ( !isDiscordant(read, estimated_insertsize, estimated_insertsize_stdev ) ) {
            version(extraVerbose) {
                if (has_minimum_quality(read, 20)) {
                    writeln("1. Has MAPQ: %d", read.mapping_quality);
                } else {
                    writefln("----> 1. Failing MAPQ: %d", read.mapping_quality);
                    }
                if (read.is_paired) {
                    writeln("2. Read is paired");
                    }
                if (!read.mate_is_unmapped) {
                    writeln("3. Mate is mapped");
                    }
                if (!is_concordant(read)) {
                    writeln("4. Is not concordant");
                }
                if (has_abnormal_isize(read, estimated_insertsize, estimated_insertsize_stdev)) {
                    writefln("5. Has abnormal isize: %d", read.template_length );
                }
            }
            continue;
        }
        bool found_read = false;
        foreach(int j, ReadInfo cread; contigreads) {
        	if(cread.name == read.name) {
//              writefln("%s %s. %d - %d  = %d", cread.name, read.name, cread.end, read.position, read.position - cread.end );
//              writefln("%s", read.toSam );
//              writefln( "Found matching pair in %s name: %s", reference_sequence.name(), read.name );
//              writefln("=================START REM (%d)============", contigreads.length);
//              ReadInfo r = contigreads[ i ];
//              writefln("Readname mem: %s", r.name);
//              writefln("Readname for: %s", cread.name);
//              writefln("Readname cur: %s", read.name);
                found_read = true;
                contigreads = std.algorithm.remove( contigreads, j );
//              writefln("=================END REM (%d)============", contigreads.length);
        		break;
    		}
    	}
        if (!found_read) {
            contigreads ~= ReadInfo(
                    read.ref_id,
                    read.position,
                    read.basesCovered,
                    read.flag, 
                    read.name
            );
        } else {
//        	double progress = cast(double) read.position / cast(double) reference_sequence.length;
//        	writefln("%d reads in %s memory, %d/%d %0.5f", contigreads.length, reference_sequence.name(), read.position, reference_sequence.length, progress);
    	}
    }
    writefln( "Unmatched reads: %d in %s, size in mem: %d", contigreads.length, reference_sequence.name(), contigreads.sizeof );
}

void scanContig(string bamfile, 
                ReferenceSequenceInfo reference_sequence, 
                ushort cutoff_min, 
                ushort cutoff_max, 
                ushort window_width,
                uint estimated_insertsize, 
                uint estimated_insertsize_stdev,
                uint max_sd_degree,
                TaskPool task_pool) {

    auto tp = new TaskPool(4);
    scope(exit) tp.finish();

    auto bam = new BamReader(bamfile, tp);
    auto input_buf_size = 16_000_000;
    bam.setBufferSize(input_buf_size);
    
    writefln("Starting analysis of %s", reference_sequence.name);
    writefln("Read bam %s for %s", bamfile, reference_sequence.name);

    extractReadInfo(bam, 
                reference_sequence, 
                cutoff_min, 
                cutoff_max, 
                window_width,
                estimated_insertsize, 
                estimated_insertsize_stdev,
                max_sd_degree, task_pool);

    writefln("Finished extraction %s", reference_sequence.name);
    
//    makeRegion(bam, 
//                reference_sequence, 
//                cutoff_min, 
//                cutoff_max, 
//                window_width,
//                estimated_insertsize, 
//                estimated_insertsize_stdev,
//                max_sd_degree, tp);
}

void printUsage() {
    int n_threads = std.parallelism.totalCPUs;
    
    stderr.writeln("Usage: yamsvp-bedregion [options] <input.bam | input.sam>");
    stderr.writeln();
    stderr.writeln("Options: -T, --threads=NTHREADS");
    stderr.writefln("                    maximum number of threads to use [%d]", n_threads);
    
    stderr.writeln("         -w, --window=WINDOWWIDTH");
    stderr.writeln("                    width of window (bp) to scan [200]");

    stderr.writeln("         -l, --cutoff_min");
    stderr.writeln("                    low threshold for reads covering a breakpoint [2]");
    stderr.writeln("         -h, --cutoff_max");
    stderr.writeln("                    high threshold for reads covering a breakpoint [90]");

    stderr.writeln("         -o, --bam_output=BAMOUTPUT");
    stderr.writeln("                    Write the resulting bam to ...");
    stderr.writeln("         -l, --compression_level");
    stderr.writeln("                    compression level 1-9 (low-high) of bam file [9]");

    stderr.writeln("         -i, --insertsize");
    stderr.writeln("                    insertsize of the input bam [450]");
    stderr.writeln("         -s, --insertsd");
    stderr.writeln("                    stdev in the insertsize [15]");
    stderr.writeln("         -r, --reads_for_histogram");
    stderr.writeln("                    reads used to determine insert-size, a minimum of 1_000_000 is recomended [1_000_000]");



    stderr.writeln("         -d, --max_sd_degree");
    stderr.writeln("                    set degree of sd allowed to mark as concordant read in the insertsize [2]");

    stderr.writeln("         -v, --verbose");
    stderr.writeln("                    Turn on verbose mode");
}

int main(string[] args) {

    int n_threads = std.parallelism.totalCPUs;
    ushort window_width = 200;

    ushort cutoff_max = 90;
    ushort cutoff_min = 2;

    uint estimated_insertsize = 450;
    uint estimated_insertsize_stdev = 15;
    uint max_sd_degree = 2;
    ubyte compression_level = 9;
    
    uint reads_for_histogram = 1_000_000;
    
    string bam_outputfilename;
    string bed_filename;

    bool verbose;

    getopt(
        args,
        std.getopt.config.caseSensitive,
        "threads|T",  &n_threads,    // numeric

        "window|w", &window_width,   // numeric

        "cutoff_max|h", &cutoff_max, // numeric
        "cutoff_min|l", &cutoff_min, // numeric

        "bamoutput|o", &bam_outputfilename,
        "compression_level|l", &compression_level, // numeric

        "insertsize|i", &estimated_insertsize, // numeric
        "insertsd|s", &estimated_insertsize_stdev, // numeric
        "reads_for_histogram|r", &reads_for_histogram, // numeric
        "max_sd_degree|d", &max_sd_degree, // numeric

        "verbose|v", &verbose,       // flag
        );

    if (args.length < 2) {
        printUsage();
        return 0;
    }

    auto task_pool = new TaskPool(n_threads);
    scope(exit) task_pool.finish();
    
    string _bamfilename = args[1];
    auto bam = new BamReader(_bamfilename, task_pool);
    auto src = new SortedBamReader( bam );
    
    // Fix index, we need this for random access
    if( !src.has_index ) {
        writeln("No bai index found, creating one now....");
        src.createIndex();
    }
    
    
    Histogram h = Histogram();
    auto _reads = src.reads;
    auto limit = reads_for_histogram;
    int counter = 0;
    foreach( BamRead read; _reads) {
    	++counter;
        if( abs(read.template_length) < 1) continue;
        if( read.mate_ref_id != read.ref_id ) continue;
        if( read.is_duplicate ) continue;
        if( read.is_secondary_alignment) continue;
        if( read.is_unmapped ) continue;
        if( !read.is_paired ) continue;

        limit--;
        h.add( std.math.abs(read.template_length) );
        if (limit == 0){
            break;
        }
    }
    
    auto stat = h._recompute();
    
    writefln("Histogram width: %d", h.width);
    writefln("Histogram records: %d", h.size);
    writefln("Histogram median: %s", h.median);
    writefln("Histogram q25: %s", h.q25);
    writefln("Histogram q50: %s", h.q50);
    writefln("Histogram q75: %s", h.q75);
    writefln("Histogram mean: %s", h.mean);
    writefln("Histogram sd: %s", h.sd);
    
    foreach (refseq; src.reference_sequences) {
//        writefln("Starting thread for %s", refseq.name);
        auto t = task!scanContig(_bamfilename, 
                                refseq, 
                                cutoff_min, 
                                cutoff_max, 
                                window_width,
                                estimated_insertsize, 
                                estimated_insertsize_stdev,
                                max_sd_degree,
                                task_pool);
        task_pool.put(t);
    }
    return 0;
}
